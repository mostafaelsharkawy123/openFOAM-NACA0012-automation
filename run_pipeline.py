#!/usr/bin/env python3
"""
Generic OpenFOAM Automated Pipeline (v2 - generalized)
--------------------------------------------------------
Converts an ANSYS Fluent mesh (.msh), lets you choose Steady-State or
Transient and Laminar or Turbulent (RAS), auto-generates ALL required
boundary conditions and case dictionaries, picks the right solver
(icoFoam / pimpleFoam / simpleFoam), runs it (serial or parallel), and
post-processes: force coefficients (Cd/Cl), VTK export, summary report.

Requirements: OpenFOAM environment must be sourced (icoFoam, fluentMeshToFoam,
etc. must be in PATH) before running this script. Tested against the
openfoam.com branch (e.g. v2512) — command-line option syntax (e.g.
'transformPoints -scale', no '-force' on reconstructPar) follows that branch.
"""

import os
import re
import sys
import time
import shutil
import glob
import subprocess

# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def get_input(prompt, default_val, val_type=float, min_val=None):
    """Prompt user for a value, falling back to default on bad/empty input."""
    while True:
        user_val = input(f"{prompt} [Default: {default_val}]: ").strip()
        if not user_val:
            result = val_type(default_val)
        else:
            try:
                result = val_type(user_val)
            except ValueError:
                print(f"  -> Invalid format. Please enter a valid {val_type.__name__}.")
                continue
        if min_val is not None and result < min_val:
            print(f"  -> Value must be >= {min_val}. Try again.")
            continue
        return result


def get_choice(prompt, options, default_index=0):
    """Prompt user to pick one of a list of string options. Returns the chosen string."""
    print(f"\n{prompt}")
    for idx, opt in enumerate(options):
        marker = " (default)" if idx == default_index else ""
        print(f"  [{idx + 1}] {opt}{marker}")
    while True:
        raw = input(f"Choice [1-{len(options)}]: ").strip()
        if not raw:
            return options[default_index]
        try:
            choice = int(raw)
        except ValueError:
            print("  -> Please enter a number.")
            continue
        if 1 <= choice <= len(options):
            return options[choice - 1]
        print(f"  -> Please enter a number between 1 and {len(options)}.")


def confirm(prompt, default_yes=True):
    suffix = "[Y/n]" if default_yes else "[y/N]"
    raw = input(f"{prompt} {suffix}: ")
    # Strip anything that isn't a plain letter (stray \r, null bytes, or other
    # terminal artifacts seen on some Windows/PowerShell setups can otherwise
    # survive a plain .strip() and make a clear "Y" fail the startswith check).
    ans = re.sub(r'[^a-zA-Z]', '', raw).lower()
    if not ans:
        return default_yes
    result = ans.startswith("y")
    if raw.strip() != ans and raw.strip().lower() != ans:
        print(f"  [Note] Raw input was {raw!r}; interpreted as "
              f"'{'yes' if result else 'no'}'.")
    return result


def check_tool(name):
    """Verify a required OpenFOAM/system executable is available."""
    if shutil.which(name) is None:
        print(f"\n[ERROR] Required executable '{name}' was not found in PATH.")
        print("        Make sure you have sourced your OpenFOAM environment, e.g.:")
        print("        source /opt/openfoam<version>/etc/bashrc")
        sys.exit(1)


def run_cmd(cmd, shell=False, description=""):
    """Run a subprocess command with unified, friendly error handling."""
    if description:
        print(f"\n{description}")
    try:
        subprocess.run(cmd, shell=shell, check=True)
    except FileNotFoundError:
        exe = cmd if shell else cmd[0]
        print(f"[ERROR] Command not found: {exe}")
        print("        Check that your OpenFOAM environment is sourced correctly.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed (exit code {e.returncode}): {cmd}")
        sys.exit(1)


def run_checkmesh_and_parse(case_dir):
    """
    Run checkMesh, stream its output to the console, and parse a few useful
    numbers out of it: the minimum cell volume (used to suggest a stable
    deltaT) and the overall domain bounding box (used to recover the real
    span/thickness of a 2D-extruded mesh for force-coefficient
    normalization). Returns a dict; any entry may be missing if it
    couldn't be parsed from this OpenFOAM version's output format.
    """
    print("\nRunning checkMesh...")
    try:
        result = subprocess.run(
            ["checkMesh"], cwd=case_dir, check=True,
            capture_output=True, text=True,
        )
        output = result.stdout
        print(output)
    except FileNotFoundError:
        print("[ERROR] Command not found: checkMesh")
        print("        Check that your OpenFOAM environment is sourced correctly.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(e.stdout or "")
        print(e.stderr or "")
        print(f"[ERROR] checkMesh failed (exit code {e.returncode}).")
        sys.exit(1)

    info = {}
    FLOAT = r'(-?\d+\.?\d*(?:[eE][+-]?\d+)?)'

    match = re.search(rf'Min volume\s*=\s*{FLOAT}', output)
    if match:
        info["min_volume"] = float(match.group(1))

    match = re.search(rf'Max volume\s*=\s*{FLOAT}', output)
    if match:
        info["max_volume"] = float(match.group(1))

    match = re.search(rf'Minimum face area\s*=\s*{FLOAT}\.\s*Maximum face area\s*=\s*{FLOAT}', output)
    if match:
        info["min_face_area"] = float(match.group(1))
        info["max_face_area"] = float(match.group(2))

    match = re.search(rf'Mesh non-orthogonality Max:\s*{FLOAT}\s*average:\s*{FLOAT}', output)
    if match:
        info["max_non_orthogonality"] = float(match.group(1))
        info["avg_non_orthogonality"] = float(match.group(2))

    match = re.search(rf'Max aspect ratio\s*=\s*{FLOAT}', output)
    if match:
        info["max_aspect_ratio"] = float(match.group(1))

    match = re.search(rf'Max skewness\s*=\s*{FLOAT}', output)
    if match:
        info["max_skewness"] = float(match.group(1))

    bbox_match = re.search(r'Overall domain bounding box \(([^)]+)\) \(([^)]+)\)', output)
    if bbox_match:
        try:
            bbox_min = [float(x) for x in bbox_match.group(1).split()]
            bbox_max = [float(x) for x in bbox_match.group(2).split()]
            if len(bbox_min) == 3 and len(bbox_max) == 3:
                info["z_span"] = abs(bbox_max[2] - bbox_min[2])
        except ValueError:
            pass

    return info


def recommend_solver_settings(quality):
    """Map checkMesh-derived mesh-quality metrics to a set of solver
    settings, instead of assuming one fixed mesh quality every time. This
    lets the pipeline adapt automatically to whatever mesh it's given —
    a clean, well-behaved mesh gets standard settings and fast convergence;
    a rougher one automatically gets more correction and gentler relaxation
    instead of silently diverging.

    Returns a dict:
        non_orth_scheme        -- snGrad/laplacian correction scheme
        n_non_orth_correctors  -- SIMPLE non-orthogonal corrector passes
        relax_p, relax_U, relax_turb -- under-relaxation factors
        use_simplec            -- whether to use SIMPLEC (consistent=yes),
                                   OpenFOAM's closest standard equivalent
                                   to a tightly-coupled pressure-velocity
                                   scheme; used as the default since it is
                                   generally at least as robust as plain
                                   SIMPLE and often converges faster.
    """
    max_nonortho = quality.get("max_non_orthogonality", 0.0)
    max_ar = quality.get("max_aspect_ratio", 1.0)
    max_skew = quality.get("max_skewness", 0.0)

    if max_nonortho < 50:
        non_orth_scheme = "orthogonal"
        n_correctors = 0
    elif max_nonortho < 70:
        non_orth_scheme = "limited 0.5"
        n_correctors = 1
    else:
        non_orth_scheme = "limited 0.333"
        n_correctors = 2

    # Overall mesh-quality tier drives how gently we relax.
    if max_ar < 100 and max_skew < 2.0:
        tier = "good"
    elif max_ar < 300 and max_skew < 4.0:
        tier = "moderate"
    else:
        tier = "rough"

    use_simplec = True

    # SIMPLEC's consistency correction is what allows near-unity pressure
    # relaxation; plain SIMPLE needs far gentler values (p ~0.3) to stay
    # stable. Values below follow common OpenFOAM tutorial practice for
    # each algorithm, tightened progressively for rougher meshes.
    if use_simplec:
        relax_by_tier = {
            "good":     (1.0, 0.9, 0.8),
            "moderate": (0.7, 0.7, 0.6),
            "rough":    (0.4, 0.5, 0.4),
        }
    else:
        relax_by_tier = {
            "good":     (0.3, 0.7, 0.7),
            "moderate": (0.2, 0.5, 0.5),
            "rough":    (0.15, 0.3, 0.3),
        }
    relax_p, relax_U, relax_turb = relax_by_tier[tier]

    return {
        "non_orth_scheme": non_orth_scheme,
        "n_non_orth_correctors": n_correctors,
        "relax_p": relax_p,
        "relax_U": relax_U,
        "relax_turb": relax_turb,
        "use_simplec": use_simplec,
        "quality_tier": tier,
    }



# --------------------------------------------------------------------------
# Mesh handling
# --------------------------------------------------------------------------

def find_msh_file(case_dir):
    msh_files = sorted(f for f in os.listdir(case_dir) if f.lower().endswith(".msh"))
    if not msh_files:
        return None
    if len(msh_files) == 1:
        return msh_files[0]

    print("\nMultiple .msh files found:")
    for idx, f in enumerate(msh_files):
        print(f" [{idx + 1}] {f}")
    choice = get_input("Select mesh file number", 1, int, min_val=1)
    if choice > len(msh_files):
        print("  -> Out of range, defaulting to the first file.")
        choice = 1
    return msh_files[choice - 1]


def parse_boundary_patches(case_dir):
    """
    Robustly parse constant/polyMesh/boundary: for every top-level
    'name { ... }' block (except the FoamFile header), extract the
    patch name and its 'type' entry regardless of ordering inside the block.
    """
    boundary_path = os.path.join(case_dir, "constant", "polyMesh", "boundary")
    if not os.path.exists(boundary_path):
        print(f"[ERROR] Boundary file not found at {boundary_path}")
        print("        Mesh conversion likely failed or produced an unexpected layout.")
        sys.exit(1)

    with open(boundary_path, "r") as f:
        content = f.read()

    patches = []
    for match in re.finditer(r'(\w+)\s*\{([^{}]*)\}', content, re.DOTALL):
        name, body = match.group(1), match.group(2)
        if name == "FoamFile":
            continue
        type_match = re.search(r'type\s+(\w+)\s*;', body)
        patch_type = type_match.group(1) if type_match else "patch"
        patches.append((name, patch_type))

    if not patches:
        print("[ERROR] No boundary patches could be parsed. Check the mesh conversion output.")
        sys.exit(1)

    return patches


def detect_required_rotation(mesh_file):
    """Read the airfoil's actual surface coordinates straight out of the raw
    Fluent .msh file (before conversion) and determine what -rotate-z angle
    (if any) is needed so the chord ends up along +X with the leading edge
    upstream (smaller X) and the trailing edge downstream (larger X). Only
    meaningful for external-aero (single wall + farfield) cases; returns
    None (skip rotation) if no such zone can be identified, which is
    harmless for other case types. Handles all four cases (already correct,
    needs 90, -90, or 180) instead of assuming one fixed mesh orientation."""
    try:
        with open(mesh_file, "r", errors="ignore") as f:
            content = f.read()

        zone_decls = re.findall(r'\(45 \((\d+) (\S+) (\S+)\)\(\)\)', content)
        target_zone_id = None
        for zid, ztype, zname in zone_decls:
            if "foil" in zname.lower():
                target_zone_id = zid
                break
        if target_zone_id is None:
            for zid, ztype, zname in zone_decls:
                name_lower = zname.lower()
                if ztype == "wall" and not any(
                    k in name_lower for k in ("farfield", "freestream", "outer", "far")
                ):
                    target_zone_id = zid
                    break
        if target_zone_id is None:
            return None

        zone_hex = format(int(target_zone_id), "x")
        face_match = re.search(
            rf'\(13 \({zone_hex} ([0-9a-fA-F]+) ([0-9a-fA-F]+) \d+(?: \d+)?\)\((.*?)\)\)',
            content, re.DOTALL,
        )
        if not face_match:
            return None
        face_lines = [l.strip() for l in face_match.group(3).strip().split("\n") if l.strip()]

        node_ids = set()
        for line in face_lines:
            for tok in line.split()[:2]:
                try:
                    node_ids.add(int(tok, 16))
                except ValueError:
                    pass
        if not node_ids:
            return None

        node_headers = re.findall(
            r'\(10 \(([0-9a-fA-F]+) ([0-9a-fA-F]+) ([0-9a-fA-F]+) 1 2\)\((.*?)\)\)',
            content, re.DOTALL,
        )
        coords = {}
        for _zh, start_hex, _end_hex, body in node_headers:
            nums = re.findall(r'-?\d+\.?\d*(?:[eE][+-]?\d+)?', body)
            start_dec = int(start_hex, 16)
            for i in range(0, len(nums) - 1, 2):
                coords[start_dec + i // 2] = (float(nums[i]), float(nums[i + 1]))

        pts = [coords[n] for n in node_ids if n in coords]
        if len(pts) < 10:
            return None

        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        x_span, y_span = max(xs) - min(xs), max(ys) - min(ys)
        tol = 0.01 * max(x_span, y_span)

        def thickness_near(axis_is_x, target_val):
            band = [p for p in pts if abs((p[0] if axis_is_x else p[1]) - target_val) < tol]
            other = [(p[1] if axis_is_x else p[0]) for p in band]
            return (max(other) - min(other)) if other else 0.0

        if x_span >= y_span:
            t_at_min = thickness_near(True, min(xs))
            t_at_max = thickness_near(True, max(xs))
            return 0 if t_at_min >= t_at_max else 180
        else:
            t_at_min = thickness_near(False, min(ys))
            t_at_max = thickness_near(False, max(ys))
            return -90 if t_at_min >= t_at_max else 90
    except Exception:
        return None


def classify_patches(patches):
    """Return a dict: patch_name -> ('inlet'|'outlet'|'wall'|'empty'|'symmetry'|'other')"""
    classification = {}
    for name, ptype in patches:
        name_lower = name.lower()
        if ptype == "empty":
            classification[name] = "empty"
        elif ptype == "symmetry" or ptype == "symmetryPlane" or "sym" in name_lower:
            classification[name] = "symmetry"
        elif "inlet" in name_lower:
            classification[name] = "inlet"
        elif "outlet" in name_lower:
            classification[name] = "outlet"
        elif ptype == "wall" or "wall" in name_lower:
            classification[name] = "wall"
        elif "front" in name_lower or "back" in name_lower:
            # Fallback heuristic for pseudo-2D cases where fluentMeshToFoam
            # did not already tag the patch as 'empty'.
            classification[name] = "empty"
        else:
            classification[name] = "wall"  # safe default: no-slip wall
    return classification


# --------------------------------------------------------------------------
# Field / dictionary generation: U, p
# --------------------------------------------------------------------------

def generate_u_dict(case_dir, classification, velocity):
    u_file = os.path.join(case_dir, "0", "U")
    os.makedirs(os.path.dirname(u_file), exist_ok=True)

    patch_bcs = []
    for name, kind in classification.items():
        if kind == "inlet":
            bc = (f"    {name}\n    {{\n        type            fixedValue;\n"
                  f"        value           uniform ({velocity} 0 0);\n    }}")
        elif kind == "outlet":
            bc = f"    {name}\n    {{\n        type            zeroGradient;\n    }}"
        elif kind == "empty":
            bc = f"    {name}\n    {{\n        type            empty;\n    }}"
        elif kind == "symmetry":
            bc = f"    {name}\n    {{\n        type            symmetry;\n    }}"
        else:  # wall
            bc = f"    {name}\n    {{\n        type            noSlip;\n    }}"
        patch_bcs.append(bc)

    boundary_block = "\n".join(patch_bcs)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({velocity} 0 0);

boundaryField
{{
{boundary_block}
}}
"""
    with open(u_file, "w") as f:
        f.write(content)


def generate_p_dict(case_dir, classification):
    p_file = os.path.join(case_dir, "0", "p")
    os.makedirs(os.path.dirname(p_file), exist_ok=True)

    patch_bcs = []
    for name, kind in classification.items():
        if kind == "outlet":
            bc = (f"    {name}\n    {{\n        type            fixedValue;\n"
                  f"        value           uniform 0;\n    }}")
        elif kind == "empty":
            bc = f"    {name}\n    {{\n        type            empty;\n    }}"
        elif kind == "symmetry":
            bc = f"    {name}\n    {{\n        type            symmetry;\n    }}"
        else:  # inlet or wall -> zeroGradient
            bc = f"    {name}\n    {{\n        type            zeroGradient;\n    }}"
        patch_bcs.append(bc)

    boundary_block = "\n".join(patch_bcs)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{{
{boundary_block}
}}
"""
    with open(p_file, "w") as f:
        f.write(content)


# --------------------------------------------------------------------------
# Turbulence: momentumTransport + k / epsilon|omega / nut
# --------------------------------------------------------------------------

TURB_DIMENSIONS = {
    "k": "[0 2 -2 0 0 0 0]",
    "epsilon": "[0 2 -3 0 0 0 0]",
    "omega": "[0 0 -1 0 0 0 0]",
    "nut": "[0 2 -1 0 0 0 0]",
}

TURB_WALL_FUNCTION = {
    "k": "kqRWallFunction",
    "epsilon": "epsilonWallFunction",
    "omega": "omegaWallFunction",
    # nutkWallFunction only works for y+ in the ~30-300 log-law range.
    # nutUSpaldingWallFunction blends continuously from the viscous
    # sublayer through the log-law region, so it stays valid whether a
    # given mesh's first cell lands at y+~1 or y+~200 — this lets the
    # pipeline adapt to whatever near-wall resolution the mesh actually
    # has, instead of requiring the user to tune the mesh to match one
    # fixed wall-function assumption.
    "nut": "nutUSpaldingWallFunction",
}


def generate_turbulence_properties(case_dir, turbulent, ras_model):
    """Generates constant/turbulenceProperties for OpenFOAM v2512 compatibility."""
    tp_file = os.path.join(case_dir, "constant", "turbulenceProperties")
    os.makedirs(os.path.dirname(tp_file), exist_ok=True)
    sim_type = "RAS" if turbulent else "laminar"
    model_str = ras_model if (turbulent and ras_model) else "laminar"
    
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}}

simulationType  {sim_type};

RAS
{{
    RASModel        {model_str};
    turbulence      on;
    printCoeffs     on;
}}
"""
    with open(tp_file, "w") as f:
        f.write(content)


def estimate_turbulence_quantities(velocity, turb_intensity_pct, length_scale, ras_model):
    """Standard engineering estimate of inlet k and epsilon/omega from
    turbulence intensity and a turbulent length scale."""
    Cmu = 0.09
    intensity = turb_intensity_pct / 100.0
    k = 1.5 * (velocity * intensity) ** 2
    if ras_model == "kEpsilon":
        second = (Cmu ** 0.75) * (k ** 1.5) / length_scale
    else:  # kOmegaSST
        second = (k ** 0.5) / ((Cmu ** 0.25) * length_scale)
    return k, second


def generate_turbulence_field(case_dir, field_name, internal_value, classification):
    """Generic writer for k, epsilon, or omega field files."""
    f_file = os.path.join(case_dir, "0", field_name)
    os.makedirs(os.path.dirname(f_file), exist_ok=True)
    wall_type = TURB_WALL_FUNCTION[field_name]

    patch_bcs = []
    for name, kind in classification.items():
        if kind == "inlet":
            bc = (f"    {name}\n    {{\n        type            fixedValue;\n"
                  f"        value           uniform {internal_value};\n    }}")
        elif kind == "outlet":
            bc = (f"    {name}\n    {{\n        type            inletOutlet;\n"
                  f"        inletValue      uniform {internal_value};\n"
                  f"        value           uniform {internal_value};\n    }}")
        elif kind == "empty":
            bc = f"    {name}\n    {{\n        type            empty;\n    }}"
        elif kind == "symmetry":
            bc = f"    {name}\n    {{\n        type            symmetry;\n    }}"
        else:  # wall
            bc = (f"    {name}\n    {{\n        type            {wall_type};\n"
                  f"        value           uniform {internal_value};\n    }}")
        patch_bcs.append(bc)

    boundary_block = "\n".join(patch_bcs)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      {field_name};
}}

dimensions      {TURB_DIMENSIONS[field_name]};
internalField   uniform {internal_value};

boundaryField
{{
{boundary_block}
}}
"""
    with open(f_file, "w") as f:
        f.write(content)


def generate_nut_field(case_dir, classification):
    """nut is derived, not solved directly: 'calculated' everywhere except
    walls, which use a wall function."""
    f_file = os.path.join(case_dir, "0", "nut")
    os.makedirs(os.path.dirname(f_file), exist_ok=True)

    patch_bcs = []
    for name, kind in classification.items():
        if kind == "empty":
            bc = f"    {name}\n    {{\n        type            empty;\n    }}"
        elif kind == "symmetry":
            bc = f"    {name}\n    {{\n        type            symmetry;\n    }}"
        elif kind == "wall":
            bc = (f"    {name}\n    {{\n        type            nutUSpaldingWallFunction;\n"
                  f"        value           uniform 0;\n    }}")
        else:  # inlet / outlet
            bc = (f"    {name}\n    {{\n        type            calculated;\n"
                  f"        value           uniform 0;\n    }}")
        patch_bcs.append(bc)

    boundary_block = "\n".join(patch_bcs)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      {TURB_DIMENSIONS['nut']};
internalField   uniform 0;

boundaryField
{{
{boundary_block}
}}
"""
    with open(f_file, "w") as f:
        f.write(content)


# --------------------------------------------------------------------------
# transportProperties
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# transportProperties & physicalProperties
# --------------------------------------------------------------------------

def generate_transport_properties(case_dir, nu):
    """Generates transportProperties and physicalProperties for OpenFOAM v2512 compatibility."""
    filenames = ["transportProperties", "physicalProperties"]
    
    for filename in filenames:
        tp_file = os.path.join(case_dir, "constant", filename)
        os.makedirs(os.path.dirname(tp_file), exist_ok=True)
        content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      {filename};
}}

transportModel  Newtonian;
viscosityModel  Newtonian;

nu              [0 2 -1 0 0 0 0] {nu};
"""
        with open(tp_file, "w") as f:
            f.write(content)


# --------------------------------------------------------------------------
# controlDict (+ optional forceCoeffs functions block)
# --------------------------------------------------------------------------

def generate_control_dict(case_dir, application="icoFoam", steady=False, delta_t=0.00005,
                           end_time=10.0, write_interval=200, functions_block="",
                           adjustable_dt=False, max_co=0.8, max_delta_t=1.0,
                           start_from="startTime", purge_write=0):
    cd_file = os.path.join(case_dir, "system", "controlDict")
    os.makedirs(os.path.dirname(cd_file), exist_ok=True)

    if steady:
        dt_line = "deltaT          1;"
        adjust_block = ""
    else:
        dt_line = f"deltaT          {delta_t};"
        if adjustable_dt:
            adjust_block = (f"adjustTimeStep  yes;\n"
                             f"maxCo           {max_co};\n"
                             f"maxDeltaT       {max_delta_t};\n")
        else:
            adjust_block = ""

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}

application     {application};
startFrom       {start_from};
startTime       0;
stopAt          endTime;
endTime         {end_time};
{dt_line}
writeControl    timeStep;
writeInterval   {write_interval};
purgeWrite      {purge_write};
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
{adjust_block}{functions_block}"""
    with open(cd_file, "w") as f:
        f.write(content)


def build_force_coeffs_block(wall_patches, velocity, lRef, Aref, liftDir, dragDir):
    """Return a controlDict 'functions{}' block that computes force
    coefficients (Cd, Cl) on the given wall patches at every write.
    Aref must be computed by the caller from the mesh's REAL span (e.g.
    from checkMesh's bounding box) — this function does not assume any
    particular mesh's dimensions."""
    if not wall_patches:
        return ""
    patches_str = " ".join(wall_patches)
    return f"""
functions
{{
    forceCoeffs1
    {{
        type            forceCoeffs;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   1;
        patches         ({patches_str});
        rho             rhoInf;
        rhoInf          1;
        liftDir         ({liftDir});
        dragDir         ({dragDir});
        CofR            (0 0 0);
        pitchAxis       (0 0 1);
        magUInf         {velocity};
        lRef            {lRef};
        Aref            {Aref};
    }}
}}
"""


def build_yplus_block(wall_patches):
    """Return a controlDict 'functions{}' block that reports y+ (min/max/avg
    per wall patch) at every write, so wall-function assumptions can be
    checked against the near-wall mesh resolution. Auto-enabled whenever the
    case is turbulent and has wall patches — no user input needed."""
    if not wall_patches:
        return ""
    patches_str = " ".join(wall_patches)
    return f"""
functions
{{
    yPlus1
    {{
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        patches         ({patches_str});
    }}
}}
"""


def read_latest_yplus_summary(case_dir):
    """Return the header + last data row from the most recent
    postProcessing/yPlus1/<time>/yPlus.dat, or None if not found."""
    pattern = os.path.join(case_dir, "postProcessing", "yPlus1", "*", "yPlus.dat")
    files = glob.glob(pattern)
    if not files:
        return None

    def time_key(path):
        try:
            return float(os.path.basename(os.path.dirname(path)))
        except ValueError:
            return -1.0

    latest_file = max(files, key=time_key)
    with open(latest_file, "r", errors="ignore") as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]
    header_lines = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if not line.startswith("#")]
    return {
        "file": latest_file,
        "header": header_lines,
        "last_row": data_lines[-1] if data_lines else None,
    }


# --------------------------------------------------------------------------
# fvSchemes / fvSolution (depend on solver: icoFoam / pimpleFoam / simpleFoam)
# --------------------------------------------------------------------------

def generate_fv_schemes(case_dir, solver, turbulent, ras_model=None, non_orth_scheme="limited 0.5"):
    fv_file = os.path.join(case_dir, "system", "fvSchemes")
    os.makedirs(os.path.dirname(fv_file), exist_ok=True)

    steady = (solver == "simpleFoam")
    ddt_line = "default         steadyState;" if steady else "default         Euler;"

    # div(phi,U): bounded/upwind-biased for SIMPLE (steady) and for turbulent
    # PIMPLE runs (stability); plain linear for the classic laminar icoFoam case.
    if steady:
        div_u_line = "div(phi,U)      bounded Gauss linearUpwind grad(U);"
    elif turbulent:
        div_u_line = "div(phi,U)      bounded Gauss linearUpwind grad(U);"
    else:
        div_u_line = "div(phi,U)      Gauss linear;"

    div_lines = [div_u_line]
    if solver != "icoFoam":
        # icoFoam's momentum equation has no nuEff term (no turbulence support);
        # pimpleFoam/simpleFoam always include it (nut = 0 when laminar).
        div_lines.append("div((nuEff*dev2(T(grad(U)))))  Gauss linear;")
    if turbulent:
        second_field = "epsilon" if ras_model == "kEpsilon" else "omega"
        turb_div = "bounded Gauss limitedLinear 1" if steady else "bounded Gauss upwind"
        div_lines.append(f"div(phi,k)      {turb_div};")
        div_lines.append(f"div(phi,{second_field})  {turb_div};")
    div_block = "\n    ".join(div_lines)

    wall_dist_block = ""
    if turbulent and ras_model == "kOmegaSST":
        wall_dist_block = "\nwallDist\n{\n    method meshWave;\n}\n"

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

ddtSchemes
{{
    {ddt_line}
}}

gradSchemes
{{
    default         Gauss linear;
    grad(p)         Gauss linear;
}}

divSchemes
{{
    default         none;
    {div_block}
}}

laplacianSchemes
{{
    default         Gauss linear {non_orth_scheme};
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         {non_orth_scheme};
}}
{wall_dist_block}"""
    with open(fv_file, "w") as f:
        f.write(content)


def generate_fv_solution(case_dir, solver, turbulent, ras_model=None, n_outer_correctors=2,
                          n_non_orth_correctors=2, relax_p=0.15, relax_U=0.3, relax_turb=0.3,
                          use_simplec=False):
    sol_file = os.path.join(case_dir, "system", "fvSolution")
    os.makedirs(os.path.dirname(sol_file), exist_ok=True)

    steady = (solver == "simpleFoam")
    second_field = None
    if turbulent:
        second_field = "epsilon" if ras_model == "kEpsilon" else "omega"

    # ---- solvers{} block ----
    solvers = []
    if steady:
        solvers.append("""    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.1;
        smoother        GaussSeidel;
    }""")
        solvers.append("""    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.1;
    }""")
        if turbulent:
            solvers.append(f"""    "(k|{second_field})"
    {{
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.1;
    }}""")
    else:
        solvers.append("""    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.05;
    }

    pFinal
    {
        $p;
        relTol          0;
    }""")
        solvers.append("""    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }""")
        if solver == "pimpleFoam":
            solvers.append("""    UFinal
    {
        $U;
        relTol          0;
    }""")
        if turbulent:
            solvers.append(f"""    "(k|{second_field})"
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }}

    "(k|{second_field})Final"
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0;
    }}""")
    solvers_block = "\n\n".join(solvers)

    # ---- algorithm control block ----
    if steady:
        rc_lines = ["        p               1e-6;", "        U               1e-6;"]
        if turbulent:
            rc_lines.append(f'        "(k|{second_field})"  1e-6;')
        rc_block = "\n".join(rc_lines)

        eq_lines = [f"        U               {relax_U};"]
        if turbulent:
            second_field = "omega" if "omega" in str(ras_model).lower() else "epsilon"
            eq_lines.append(f"        \"(k|{second_field})\" {relax_turb};")
        eq_block = "\n".join(eq_lines)

        consistent_word = "yes" if use_simplec else "no"

        algo_block = f"""SIMPLE
{{
    nNonOrthogonalCorrectors {n_non_orth_correctors};
    consistent          {consistent_word};

    residualControl
    {{
{rc_block}
    }}
}}

relaxationFactors
{{
    fields
    {{
        p               {relax_p};
    }}
    equations
    {{
{eq_block}
    }}
}}
"""
    elif solver == "pimpleFoam":
        algo_block = f"""PIMPLE
{{
    nOuterCorrectors {n_outer_correctors};
    nCorrectors      1;
    nNonOrthogonalCorrectors 1;
    pRefCell         0;
    pRefValue        0;
}}
"""
    else:  # icoFoam
        algo_block = """PISO
{
    nCorrectors     2;
    nNonOrthogonalCorrectors 1;
    pRefCell        0;
    pRefValue       0;
}
"""

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}}

solvers
{{
{solvers_block}
}}

{algo_block}"""
    with open(sol_file, "w") as f:
        f.write(content)


def generate_decompose_par_dict(case_dir, n_cores):
    dp_file = os.path.join(case_dir, "system", "decomposeParDict")
    os.makedirs(os.path.dirname(dp_file), exist_ok=True)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}}

numberOfSubdomains {n_cores};
method          scotch;
"""
    with open(dp_file, "w") as f:
        f.write(content)


# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------

def clean_previous_run(case_dir):
    """Remove processor directories and any numeric time-step folders
    (except '0') left over from a previous run, so stale results don't
    linger after re-running with different endTime/writeInterval."""
    for p_dir in glob.glob(os.path.join(case_dir, "processor*")):
        shutil.rmtree(p_dir)

    for entry in os.listdir(case_dir):
        full_path = os.path.join(case_dir, entry)
        if os.path.isdir(full_path) and entry != "0":
            try:
                float(entry)
            except ValueError:
                continue
            shutil.rmtree(full_path)

    pp_dir = os.path.join(case_dir, "postProcessing")
    if os.path.isdir(pp_dir):
        shutil.rmtree(pp_dir)


# --------------------------------------------------------------------------
# Post-processing: log parsing, summary report
# --------------------------------------------------------------------------

def parse_solver_log(log_path):
    """Extract final residuals (for whichever fields were actually solved),
    continuity errors and solver-reported timing from the solver log."""
    summary = {"residuals": {}}
    if not os.path.exists(log_path):
        return summary

    with open(log_path, "r", errors="ignore") as f:
        content = f.read()

    FLOAT = r'(-?\d+\.?\d*(?:[eE][+-]?\d+)?)'

    def last_match(pattern, cast=str):
        matches = re.findall(pattern, content)
        return cast(matches[-1]) if matches else None

    summary["final_time"] = last_match(rf'^Time = {FLOAT}', float)

    field_pattern = re.compile(
        rf'Solving for (\w+), Initial residual = {FLOAT}, Final residual = {FLOAT}'
    )
    for field, init_res, final_res in field_pattern.findall(content):
        summary["residuals"][field] = (float(init_res), float(final_res))

    cont_matches = re.findall(
        rf'time step continuity errors : sum local = {FLOAT}, global = {FLOAT}, cumulative = {FLOAT}',
        content,
    )
    if cont_matches:
        local_e, global_e, cum_e = cont_matches[-1]
        summary["continuity_local"] = float(local_e)
        summary["continuity_global"] = float(global_e)
        summary["continuity_cumulative"] = float(cum_e)

    exec_matches = re.findall(rf'ExecutionTime = {FLOAT} s\s+ClockTime = (\d+) s', content)
    if exec_matches:
        exec_t, clock_t = exec_matches[-1]
        summary["solver_execution_time"] = float(exec_t)
        summary["solver_clock_time"] = int(clock_t)

    return summary


def write_summary_report(case_dir, params, log_summary, wall_clock_seconds, outputs, yplus_summary=None):
    """Write a plain-text summary report of the whole run to the case directory."""
    report_path = os.path.join(case_dir, "simulation_report.txt")
    lines = []
    lines.append("=" * 60)
    lines.append("        OpenFOAM Simulation Report")
    lines.append("=" * 60)
    lines.append(f"Case directory : {case_dir}")
    lines.append(f"Generated at   : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("--- Case Setup ---")
    lines.append(f"Mesh file            : {params.get('mesh_file')}")
    lines.append(f"Mesh scale factor    : {params.get('scale')}")
    lines.append(f"Boundary patches     : {params.get('n_patches')}")
    lines.append(f"Simulation type      : {params.get('simulation_type')}")
    lines.append(f"Flow regime          : {params.get('flow_regime')}")
    if params.get("ras_model"):
        lines.append(f"Turbulence model     : {params.get('ras_model')}")
    lines.append(f"Solver               : {params.get('application')}")
    lines.append(f"Cores used           : {params.get('n_cores')}")
    lines.append("")
    lines.append("--- Parameters ---")
    lines.append(f"Inlet velocity U     : {params.get('velocity')} m/s")
    lines.append(f"Kinematic viscosity  : {params.get('nu')} m^2/s")
    if params.get("simulation_type") == "Transient":
        lines.append(f"Time step deltaT     : {params.get('delta_t')} s (initial, if adjustable)")
        lines.append(f"End time             : {params.get('end_time')} s")
    else:
        lines.append(f"Iterations (endTime) : {int(params.get('end_time'))}")
    lines.append("")
    lines.append("--- Run Results ---")
    lines.append(f"Wall-clock runtime (measured)   : {wall_clock_seconds:.1f} s")
    if "solver_execution_time" in log_summary:
        lines.append(f"Solver ExecutionTime            : {log_summary['solver_execution_time']} s")
        lines.append(f"Solver ClockTime                : {log_summary['solver_clock_time']} s")
    if "final_time" in log_summary:
        lines.append(f"Final simulation time/iteration : {log_summary['final_time']}")
    lines.append("")
    lines.append("Final residuals (last time step / iteration):")
    for field, (init_r, final_r) in log_summary.get("residuals", {}).items():
        lines.append(f"  {field:<10}: initial = {init_r:.3e}, final = {final_r:.3e}")
    if "continuity_global" in log_summary:
        lines.append("")
        lines.append("Continuity errors (last time step):")
        lines.append(f"  local = {log_summary['continuity_local']:.3e}, "
                      f"global = {log_summary['continuity_global']:.3e}, "
                      f"cumulative = {log_summary['continuity_cumulative']:.3e}")
    if yplus_summary:
        lines.append("")
        lines.append("--- Wall y+ Check ---")
        for h in yplus_summary["header"]:
            lines.append(f"  {h}")
        if yplus_summary["last_row"]:
            lines.append(f"  {yplus_summary['last_row']}")
        lines.append("  Target: y+ ~30-300 for standard wall functions, y+ < 1 for low-Re treatment.")
    lines.append("")
    lines.append("--- Post-processing Outputs ---")
    for label, path in outputs.items():
        lines.append(f"  {label:<20}: {path}")
    lines.append("=" * 60)

    with open(report_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return report_path


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("  Generic OpenFOAM Automated Pipeline")
    print("=" * 65)

    case_dir = os.getcwd()
    check_tool("fluentMeshToFoam")

    clean_previous_run(case_dir)

    # ---- Case type selection (drives everything downstream) ----
    simulation_type = get_choice(
        "Simulation type:", ["Steady-State", "Transient"], default_index=1
    )
    flow_regime = get_choice(
        "Flow regime:", ["Laminar", "Turbulent (RAS)"], default_index=0
    )
    turbulent = flow_regime.startswith("Turbulent")
    steady = (simulation_type == "Steady-State")

    ras_model = None
    if turbulent:
        ras_choice = get_choice(
            "RAS turbulence model:", ["k-epsilon", "k-omega SST"], default_index=0
        )
        ras_model = "kEpsilon" if ras_choice == "k-epsilon" else "kOmegaSST"

    if steady:
        application = "simpleFoam"
    elif turbulent:
        application = "pimpleFoam"
    else:
        application = "icoFoam"

    print(f"\n[Selected setup] {simulation_type} | {flow_regime}"
          + (f" ({ras_model})" if ras_model else "") + f" -> solver: {application}")

    mesh_file = find_msh_file(case_dir)
    if not mesh_file:
        mesh_file = input("No .msh file detected. Enter filename manually: ").strip()
        if not os.path.exists(mesh_file):
            print(f"[ERROR] Specified mesh file '{mesh_file}' does not exist.")
            sys.exit(1)

    print(f"\n[Selected Mesh]: {mesh_file}")

    # fluentMeshToFoam constructs a Foam::Time object internally, which requires
    # system/controlDict to already exist (its actual values are irrelevant at
    # this stage — they get overwritten with the real ones after user input below).
    generate_control_dict(case_dir, application=application, steady=steady)
    generate_fv_schemes(case_dir, application, turbulent, ras_model)
    generate_fv_solution(case_dir, application, turbulent, ras_model)

    print("\n[1/8] Converting ANSYS Fluent mesh using fluentMeshToFoam...")
    run_cmd(["fluentMeshToFoam", mesh_file])
    print("Mesh conversion completed successfully.")

    # Auto-detect whether the airfoil chord needs rotating onto +X with the
    # leading edge upstream, instead of assuming one fixed mesh orientation.
    required_rotation = detect_required_rotation(mesh_file)
    if required_rotation is None:
        print("\n[Note] Could not auto-detect an airfoil chord orientation from this mesh "
              "(may not be an airfoil case) — skipping automatic rotation.")
    elif required_rotation == 0:
        print("\n[Info] Airfoil chord already correctly aligned along +X — no rotation needed.")
    else:
        print(f"\n[Info] Rotating mesh {required_rotation} deg around Z so the chord aligns "
              f"with +X, leading edge upstream (detected from the actual airfoil geometry).")
        check_tool("transformPoints")
        run_cmd(
            ["transformPoints", "-rotate-z", str(required_rotation)],
            description=f"Rotating mesh {required_rotation} degrees..."
        )

    # ANSYS meshes are frequently exported in millimeters; OpenFOAM expects meters.
    scale: float = get_input(
        "\nMesh length-unit scale factor (e.g. 0.001 if ANSYS mesh was in mm, 1 for no change)",
        1.0, float, min_val=1e-12,
    )
    if scale != 1.0:
        check_tool("transformPoints")
        run_cmd(
            ["transformPoints", "-scale", f"({scale} {scale} {scale})"],
            description="Scaling mesh points...",
        )

    min_cell_volume = None
    if confirm("\nRun checkMesh to validate the converted mesh?", default_yes=True):
        check_tool("checkMesh")
        min_cell_volume = run_checkmesh_and_parse(case_dir).get("min_volume")
        if not confirm("\ncheckMesh finished. Continue with this mesh?", default_yes=True):
            print("Aborted by user after checkMesh review.")
            sys.exit(0)

    print("\n[3/8] Inspecting constant/polyMesh/boundary patches...")
    patches = parse_boundary_patches(case_dir)
    classification = classify_patches(patches)
    print(f"Found {len(patches)} boundary patches:")
    for name, ptype in patches:
        print(f" - {name:<20} (mesh type: {ptype:<10}) -> classified as: {classification[name]}")

    if not confirm("\nProceed with this boundary-condition classification?", default_yes=True):
        print("Aborted by user. Rename patches in ANSYS (or edit the script's")
        print("classify_patches logic) to fix misclassified boundaries, then re-run.")
        sys.exit(0)

    print("\n[4/8] Input Flow Parameters:")
    velocity = get_input("1. Inlet Velocity U (m/s)", 1.5, float)
    nu = get_input("2. Kinematic Viscosity nu (m^2/s)", 0.001, float, min_val=0.0)

    k_val = second_val = None
    if turbulent:
        print("\nTurbulence inlet conditions:")
        turb_intensity = get_input("  Turbulence intensity I (%)", 5.0, float, min_val=0.0)
        length_scale = get_input(
            "  Turbulent length scale (m) [~0.07 x inlet hydraulic diameter]", 0.01, float, min_val=1e-9
        )
        k_val, second_val = estimate_turbulence_quantities(velocity, turb_intensity, length_scale, ras_model)
        second_name = "epsilon" if ras_model == "kEpsilon" else "omega"
        print(f"  -> Estimated inlet k = {k_val:.4g} m^2/s^2, {second_name} = {second_val:.4g}")

    print("\n[5/8] Execution Parameters:")
    adjustable_dt = False
    max_co = 0.8
    max_delta_t = 1.0

    if steady:
        end_time = get_input("3. Number of Iterations", 1000, int, min_val=1)
        write_interval = get_input("4. Write results every N iterations", 100, int, min_val=1)
        delta_t = 1.0  # fixed for SIMPLE
    else:
        default_delta_t = 0.00005
        if min_cell_volume is not None and min_cell_volume > 0 and velocity > 0:
            target_co = 0.3  # conservative starting point
            estimated_dx = min_cell_volume ** (1.0 / 3.0)
            suggested_delta_t = target_co * estimated_dx / velocity
            print(f"\n[Suggestion] Smallest cell size (from checkMesh) ~ {estimated_dx:.3e} m")
            print(f"             For Co ~ {target_co}, suggested deltaT ~ {suggested_delta_t:.3e} s")
            print("             (Rough 3D-volume-based estimate — a starting point, not a guarantee.)")
            default_delta_t = suggested_delta_t
        else:
            print("\n[Note] No mesh size info available — using a generic default deltaT.")

        if application == "pimpleFoam":
            adjustable_dt = confirm(
                "\nEnable adjustable time-stepping (recommended - avoids manual deltaT tuning)?",
                default_yes=True,
            )
            if adjustable_dt:
                max_co = get_input("  Max Courant number", 0.8, float, min_val=0.0)
                max_delta_t = get_input("  Max allowed deltaT (s)", 1.0, float, min_val=0.0)

        delta_t = get_input("3. (Initial) Time Step deltaT (s)", default_delta_t, float, min_val=0.0)
        end_time = get_input("4. Simulation End Time (s)", 10.0, float, min_val=0.0)
        write_interval = get_input("5. Write results every N time steps", 200, int, min_val=1)

    n_cores = get_input("Number of CPU Cores (1 for Serial)", 1, int, min_val=1)

    if not steady:
        print("\n[Reminder] If a run without adjustable time-stepping crashes with")
        print("'Floating point exception' or diverges, the most common cause is deltaT")
        print("being too large relative to your smallest cell size (Courant number > ~1),")
        print("or an incorrect mesh scale factor.")

    # Wall patches are kept for the post-run forceCoeffs suggestion below;
    # we no longer ask about this upfront — it can be computed after the run
    # in one command if/when actually needed.
    wall_patches = [name for name, kind in classification.items() if kind == "wall"]
    functions_block = ""

    if turbulent and wall_patches:
        functions_block = build_yplus_block(wall_patches)
        print(f"\n[Auto-enabled] y+ monitoring on wall patches {wall_patches}")
        print("               (target: y+ ~30-300 for standard wall functions, y+ < 1 for low-Re")
        print("               near-wall resolution — checked automatically after the run).")

    export_vtk = confirm("\nExport results to VTK for ParaView after solving?", default_yes=True)

    print("\nGenerating OpenFOAM field and dictionary files...")
    generate_u_dict(case_dir, classification, velocity)
    generate_p_dict(case_dir, classification)
    generate_transport_properties(case_dir, nu)
    generate_momentum_transport(case_dir, turbulent, ras_model)
    if turbulent:
        second_name = "epsilon" if ras_model == "kEpsilon" else "omega"
        generate_turbulence_field(case_dir, "k", k_val, classification)
        generate_turbulence_field(case_dir, second_name, second_val, classification)
        generate_nut_field(case_dir, classification)

    generate_control_dict(
        case_dir, application=application, steady=steady, delta_t=delta_t, end_time=end_time,
        write_interval=write_interval, functions_block=functions_block,
        adjustable_dt=adjustable_dt, max_co=max_co, max_delta_t=max_delta_t,
    )
    generate_fv_schemes(case_dir, application, turbulent, ras_model)
    generate_fv_solution(case_dir, application, turbulent, ras_model)

    start_time = time.time()

    if n_cores > 1:
        check_tool("decomposePar")
        check_tool("mpirun")
        check_tool("reconstructPar")
        check_tool(application)

        generate_decompose_par_dict(case_dir, n_cores)
        run_cmd(
            ["decomposePar", "-force"],
            description=f"[6/8] Decomposing domain for {n_cores} processors...",
        )

        print(f"\n[7/8] Running {application} in parallel on {n_cores} cores (log: log.{application})...")
        run_cmd(
            f"mpirun --allow-run-as-root -np {n_cores} {application} -parallel | tee log.{application}",
            shell=True,
        )
        run_cmd(["reconstructPar"], description="Reconstructing parallel domain fields...")
        print(f"\nParallel {application} run completed successfully.")
    else:
        check_tool(application)
        print(f"\n[7/8] Running {application} in serial (log: log.{application})...")
        run_cmd(f"{application} | tee log.{application}", shell=True)
        print(f"\nSerial {application} run completed successfully.")

    wall_clock_seconds = time.time() - start_time

    # ---- Post-processing: VTK export ----
    outputs = {}
    if export_vtk:
        check_tool("foamToVTK")
        run_cmd(["foamToVTK"], description="Exporting results to VTK...")
        outputs["VTK export"] = os.path.join(case_dir, "VTK")

    if wall_patches:
        patches_str = " ".join(wall_patches)
        suggested_cmd = (
            f'postProcess -func "forceCoeffs(patches=({patches_str}), rhoInf=1, '
            f'liftDir=(0 1 0), dragDir=(1 0 0), CofR=(0 0 0), magUInf={velocity}, '
            f'lRef=<REFERENCE_LENGTH>, Aref=<REFERENCE_AREA>)"'
        )
        outputs["Forces (compute anytime with)"] = suggested_cmd

    yplus_summary = None
    if turbulent and wall_patches:
        yplus_summary = read_latest_yplus_summary(case_dir)
        if yplus_summary:
            print("\n[y+ check] Latest values "
                  f"(full history in {yplus_summary['file']}):")
            for h in yplus_summary["header"]:
                print(f"  {h}")
            if yplus_summary["last_row"]:
                print(f"  {yplus_summary['last_row']}")
            print("  Target: y+ ~30-300 for standard wall functions, y+ < 1 for low-Re treatment.")
            outputs["y+ report"] = yplus_summary["file"]

    # ---- Final summary report ----
    print("\n[8/8] Generating summary report...")
    log_summary = parse_solver_log(os.path.join(case_dir, f"log.{application}"))
    params = {
        "mesh_file": mesh_file,
        "scale": scale,
        "n_patches": len(patches),
        "simulation_type": simulation_type,
        "flow_regime": flow_regime,
        "ras_model": ras_model,
        "application": application,
        "n_cores": n_cores,
        "velocity": velocity,
        "nu": nu,
        "delta_t": delta_t,
        "end_time": end_time,
    }
    report_path = write_summary_report(case_dir, params, log_summary, wall_clock_seconds, outputs, yplus_summary)
    print(f"\nSummary report written to: {report_path}")
    print("\nPipeline finished successfully.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Exiting.")
        sys.exit(130)