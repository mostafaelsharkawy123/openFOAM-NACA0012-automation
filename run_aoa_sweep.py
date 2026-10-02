#!/usr/bin/env python3
"""
2D Airfoil Angle-of-Attack Sweep (run_aoa_sweep.py)
-----------------------------------------------------
Builds on run_pipeline.py's case-generation helpers. Converts a SINGLE
ANSYS mesh (airfoil at AoA = 0 inside a circular farfield domain) ONCE,
then sweeps a range of angles of attack by rotating the freestream
velocity vector (and the lift/drag directions) instead of the geometry.
"""

import os
import re
import sys
import csv
import math
import glob
import shutil
import signal
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_pipeline as base  # reuse mesh/solver/dictionary helpers


# --------------------------------------------------------------------------
# Airfoil-specific patch classification (farfield, not inlet/outlet)
# --------------------------------------------------------------------------

def classify_airfoil_patches(patches):
    """Return dict patch_name -> 'farfield' | 'wall' | 'empty' | 'symmetry'."""
    classification = {}
    for name, ptype in patches:
        name_lower = name.lower()
        if ptype == "empty":
            classification[name] = "empty"
        elif ptype == "symmetry" or ptype == "symmetryPlane" or "sym" in name_lower:
            classification[name] = "symmetry"
        elif "front" in name_lower or "back" in name_lower:
            classification[name] = "empty"
        elif any(k in name_lower for k in ("farfield", "outer", "far", "freestream")):
            classification[name] = "farfield"
        elif ptype == "wall" or any(k in name_lower for k in ("airfoil", "foil", "wall")):
            classification[name] = "wall"
        else:
            classification[name] = "farfield"
    return classification


def fix_farfield_patch_type(case_dir, classification):
    boundary_path = os.path.join(case_dir, "constant", "polyMesh", "boundary")
    with open(boundary_path, "r") as f:
        content = f.read()

    matches = list(re.finditer(r'(\w+)\s*\{([^{}]*)\}', content, re.DOTALL))
    fixed_any = False
    for match in reversed(matches):
        name = match.group(1)
        if name == "FoamFile" or classification.get(name) != "farfield":
            continue
        block_text = match.group(0)
        type_match = re.search(r'type\s+(\w+)\s*;', block_text)
        if type_match and type_match.group(1) != "patch":
            new_block = re.sub(r'type\s+\w+\s*;', 'type            patch;', block_text, count=1)
            content = content[:match.start()] + new_block + content[match.end():]
            fixed_any = True
            print(f"[Fix] Patch '{name}': geometric type changed "
                  f"'{type_match.group(1)}' -> 'patch' (so freestream behaves correctly, "
                  f"not as a solid wall for wall-distance/turbulence purposes).")

    if fixed_any:
        with open(boundary_path, "w") as f:
            f.write(content)
    return fixed_any


def detect_required_rotation(mesh_file):
    """Read the airfoil's actual surface coordinates straight out of the raw
    Fluent .msh file (before conversion) and determine what -rotate-z angle
    (if any) is needed so the chord ends up along +X with the leading edge
    upstream (smaller X) and the trailing edge downstream (larger X) —
    matching this script's AoA=0 -> flow along +X convention. Handles all
    four cases (already correct, needs 90, -90, or 180) instead of assuming
    one fixed mesh orientation. Identifies the airfoil zone directly from
    the raw zone names (this runs before conversion, so no classification
    is available yet). Returns the angle in degrees, or None if it
    couldn't be determined (caller should skip rotating and warn)."""
    try:
        with open(mesh_file, "r", errors="ignore") as f:
            content = f.read()

        zone_decls = re.findall(r'\(45 \((\d+) (\S+) (\S+)\)\(\)\)', content)
        target_zone_id = None
        # Prefer a name that clearly says "airfoil"/"foil"; fall back to any
        # wall-type zone that doesn't look like the farfield boundary.
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
            # Chord already on X. The leading edge (rounder nose) shows
            # thickness growing faster near its tip than the trailing edge does.
            t_at_min = thickness_near(True, min(xs))
            t_at_max = thickness_near(True, max(xs))
            return 0 if t_at_min >= t_at_max else 180
        else:
            t_at_min = thickness_near(False, min(ys))
            t_at_max = thickness_near(False, max(ys))
            # rotate-z -90: new_x = old_y  (old y=min -> smallest new_x)
            # rotate-z +90: new_x = -old_y (old y=max -> smallest new_x)
            return -90 if t_at_min >= t_at_max else 90
    except Exception:
        return None


# --------------------------------------------------------------------------
# Freestream-based field writers (U, p, k, epsilon/omega, nut)
# --------------------------------------------------------------------------

def _boundary_block(classification, farfield_line, wall_line, symmetry_line="symmetry"):
    parts = []
    for name, kind in classification.items():
        if kind == "farfield":
            parts.append(f"    {name}\n    {{\n{farfield_line}\n    }}")
        elif kind == "empty":
            parts.append(f"    {name}\n    {{\n        type            empty;\n    }}")
        elif kind == "symmetry":
            parts.append(f"    {name}\n    {{\n        type            {symmetry_line};\n    }}")
        else:  # wall
            parts.append(f"    {name}\n    {{\n{wall_line}\n    }}")
    return "\n".join(parts)


def generate_freestream_u(case_dir, classification, ux, uy):
    u_file = os.path.join(case_dir, "0", "U")
    os.makedirs(os.path.dirname(u_file), exist_ok=True)
    farfield_line = (f"        type            freestream;\n"
                      f"        freestreamValue uniform ({ux} {uy} 0);")
    wall_line = "        type            noSlip;"
    boundary_block = _boundary_block(classification, farfield_line, wall_line)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({ux} {uy} 0);

boundaryField
{{
{boundary_block}
}}
"""
    with open(u_file, "w") as f:
        f.write(content)


def generate_freestream_p(case_dir, classification):
    p_file = os.path.join(case_dir, "0", "p")
    os.makedirs(os.path.dirname(p_file), exist_ok=True)
    farfield_line = (f"        type            freestreamPressure;\n"
                      f"        freestreamValue uniform 0;")
    wall_line = "        type            zeroGradient;"
    boundary_block = _boundary_block(classification, farfield_line, wall_line)
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


def generate_freestream_turbulence_field(case_dir, field_name, internal_value, classification):
    f_file = os.path.join(case_dir, "0", field_name)
    os.makedirs(os.path.dirname(f_file), exist_ok=True)
    wall_type = base.TURB_WALL_FUNCTION[field_name]
    farfield_line = (f"        type            freestream;\n"
                      f"        freestreamValue uniform {internal_value};")
    wall_line = (f"        type            {wall_type};\n"
                 f"        value           uniform {internal_value};")
    boundary_block = _boundary_block(classification, farfield_line, wall_line)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      {field_name};
}}

dimensions      {base.TURB_DIMENSIONS[field_name]};
internalField   uniform {internal_value};

boundaryField
{{
{boundary_block}
}}
"""
    with open(f_file, "w") as f:
        f.write(content)


def generate_freestream_nut(case_dir, classification):
    f_file = os.path.join(case_dir, "0", "nut")
    os.makedirs(os.path.dirname(f_file), exist_ok=True)
    farfield_line = "        type            calculated;\n        value           uniform 0;"
    wall_line = "        type            nutUSpaldingWallFunction;\n        value           uniform 0;"
    boundary_block = _boundary_block(classification, farfield_line, wall_line)
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      {base.TURB_DIMENSIONS['nut']};
internalField   uniform 0;

boundaryField
{{
{boundary_block}
}}
"""
    with open(f_file, "w") as f:
        f.write(content)


def write_all_fields(case_dir, classification, ux, uy, k_val, second_name, second_val, turbulent):
    generate_freestream_u(case_dir, classification, ux, uy)
    generate_freestream_p(case_dir, classification)
    if turbulent:
        generate_freestream_turbulence_field(case_dir, "k", k_val, classification)
        generate_freestream_turbulence_field(case_dir, second_name, second_val, classification)
        generate_freestream_nut(case_dir, classification)


# --------------------------------------------------------------------------
# Warm start: patch only the boundaryField of already-solved field files
# --------------------------------------------------------------------------

def get_latest_time_dir(case_dir):
    times = []
    for entry in os.listdir(case_dir):
        full = os.path.join(case_dir, entry)
        if os.path.isdir(full):
            try:
                times.append((float(entry), entry))
            except ValueError:
                continue
    if not times:
        return None
    times.sort()
    return os.path.join(case_dir, times[-1][1])


def rewrite_boundary_field(field_file, new_boundary_block):
    with open(field_file, "r") as f:
        content = f.read()

    match = re.search(r'boundaryField\s*\{', content)
    if not match:
        raise RuntimeError(f"Could not find boundaryField section in {field_file}")

    head = content[:match.start()]
    new_content = (f"{head}boundaryField\n{{\n{new_boundary_block}\n}}\n"
                   f"\n// ************************************************************************* //\n")
    with open(field_file, "w") as f:
        f.write(new_content)


def warm_start_boundary_conditions(case_dir, classification, ux, uy, second_name, turbulent,
                                    k_val=None, second_val=None):
    latest_dir = get_latest_time_dir(case_dir)
    if latest_dir is None:
        raise RuntimeError("No previous time directory found for warm start.")

    farfield_u = (f"        type            freestream;\n"
                  f"        freestreamValue uniform ({ux} {uy} 0);")
    wall_u = "        type            noSlip;"
    u_path = os.path.join(latest_dir, "U")
    if os.path.exists(u_path):
        rewrite_boundary_field(u_path, _boundary_block(classification, farfield_u, wall_u))

    farfield_p = "        type            freestreamPressure;\n        freestreamValue uniform 0;"
    wall_p = "        type            zeroGradient;"
    p_path = os.path.join(latest_dir, "p")
    if os.path.exists(p_path):
        rewrite_boundary_field(p_path, _boundary_block(classification, farfield_p, wall_p))

    if turbulent:
        field_values = {"k": k_val, second_name: second_val}
        for field_name, value in field_values.items():
            f_path = os.path.join(latest_dir, field_name)
            if os.path.exists(f_path) and value is not None:
                wall_type = base.TURB_WALL_FUNCTION[field_name]
                farfield_line = (f"        type            freestream;\n"
                                  f"        freestreamValue uniform {value};")
                wall_line = f"        type            {wall_type};\n        value           uniform 0;"
                rewrite_boundary_field(f_path, _boundary_block(classification, farfield_line, wall_line))

    return latest_dir


# --------------------------------------------------------------------------
# Aerodynamic angle helpers
# --------------------------------------------------------------------------

def velocity_components(u_inf, aoa_deg):
    a = math.radians(aoa_deg)
    return u_inf * math.cos(a), u_inf * math.sin(a)


def lift_drag_dirs(aoa_deg):
    a = math.radians(aoa_deg)
    drag_dir = f"{math.cos(a):.6f} {math.sin(a):.6f} 0"
    lift_dir = f"{-math.sin(a):.6f} {math.cos(a):.6f} 0"
    return lift_dir, drag_dir


# --------------------------------------------------------------------------
# Convergence + force-coefficient extraction
# --------------------------------------------------------------------------

def check_simple_converged(log_path):
    if not os.path.exists(log_path):
        return False
    with open(log_path, "r", errors="ignore") as f:
        content = f.read()
    return "SIMPLE solution converged" in content


def find_force_coeffs_file(case_dir):
    pattern = os.path.join(case_dir, "postProcessing", "forceCoeffs1", "*", "*.dat")
    files = glob.glob(pattern)
    if not files:
        return None

    def time_key(path):
        try:
            return float(os.path.basename(os.path.dirname(path)))
        except ValueError:
            return -1.0

    return max(files, key=time_key)


def read_force_coeffs(case_dir, average_last_fraction=1.0):
    dat_file = find_force_coeffs_file(case_dir)
    if dat_file is None:
        return None

    header_names = None
    rows = []
    with open(dat_file, "r", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                tokens = line.lstrip("#").split()
                if "Cd" in tokens and "Cl" in tokens:
                    header_names = tokens
                continue
            rows.append(line.split())

    if header_names is None or not rows:
        return None

    try:
        cd_idx = header_names.index("Cd")
        cl_idx = header_names.index("Cl")
        time_idx = header_names.index("Time") if "Time" in header_names else 0
    except ValueError:
        return None

    parsed = []
    for row in rows:
        if len(row) <= max(cd_idx, cl_idx, time_idx):
            continue
        try:
            parsed.append((float(row[time_idx]), float(row[cd_idx]), float(row[cl_idx])))
        except ValueError:
            continue

    if not parsed:
        return None

    n_avg = max(1, int(len(parsed) * average_last_fraction))
    tail = parsed[-n_avg:]
    cd_vals = [r[1] for r in tail]
    cl_vals = [r[2] for r in tail]
    cd_mean = sum(cd_vals) / len(cd_vals)
    cl_mean = sum(cl_vals) / len(cl_vals)
    cd_std = (sum((v - cd_mean) ** 2 for v in cd_vals) / len(cd_vals)) ** 0.5
    cl_std = (sum((v - cl_mean) ** 2 for v in cl_vals) / len(cl_vals)) ** 0.5
    return {"Cd": cd_mean, "Cl": cl_mean, "Cd_std": cd_std, "Cl_std": cl_std, "n_samples": len(tail)}


# --------------------------------------------------------------------------
# Solver execution helpers
# --------------------------------------------------------------------------

def check_force_converged(case_dir, window=100, tol_rel=5e-4, cl_abs_floor=1e-4):
    """Force-based convergence: True if Cd and Cl barely changed over the last
    `window` iterations of the CURRENT angle's forceCoeffs data. Cd must change
    by less than tol_rel (relative); Cl must change by less than
    max(tol_rel*|Cl|, cl_abs_floor) so Cl ~ 0 (e.g. AoA=0) does not blow up
    a purely relative test. Returns False if there are fewer than `window` rows."""
    dat_file = find_force_coeffs_file(case_dir)
    if dat_file is None:
        return False
    header_names, rows = None, []
    with open(dat_file, "r", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                tokens = line.lstrip("#").split()
                if "Cd" in tokens and "Cl" in tokens:
                    header_names = tokens
                continue
            rows.append(line.split())
    if header_names is None or len(rows) < window:
        return False
    try:
        cd_idx, cl_idx = header_names.index("Cd"), header_names.index("Cl")
        tail = rows[-window:]
        cd_first, cd_last = float(tail[0][cd_idx]), float(tail[-1][cd_idx])
        cl_first, cl_last = float(tail[0][cl_idx]), float(tail[-1][cl_idx])
    except (ValueError, IndexError):
        return False
    cd_ok = abs(cd_last - cd_first) <= tol_rel * abs(cd_last)
    cl_ok = abs(cl_last - cl_first) <= max(tol_rel * abs(cl_last), cl_abs_floor)
    return cd_ok and cl_ok


def _run_timed(cmd, timeout_seconds, description=""):
    if description:
        print(f"\n{description}")
    if not timeout_seconds:
        base.run_cmd(cmd, shell=True)
        return True

    print(f"(wall-clock limit for this run: {timeout_seconds / 60:.1f} min)")
    proc = subprocess.Popen(cmd, shell=True, cwd=os.getcwd(), start_new_session=True)
    try:
        proc.wait(timeout=timeout_seconds)
        return proc.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"[!] Exceeded the {timeout_seconds / 60:.1f}-minute limit — stopping this run and "
              f"using whatever results were written so far.")
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            proc.wait(timeout=10)
        except Exception:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except Exception:
                pass
        return False


def run_solver(case_dir, application, n_cores, log_name, timeout_seconds=None):
    """Executes simpleFoam in serial or parallel mode."""
    if n_cores > 1:
        base.check_tool("decomposePar")
        base.check_tool("mpirun")
        base.check_tool("reconstructPar")
        base.check_tool(application)
        base.generate_decompose_par_dict(case_dir, n_cores)
        base.run_cmd(["decomposePar", "-force"], description=f"Decomposing for {n_cores} processors...")
        _run_timed(
            f"mpirun --allow-run-as-root -np {n_cores} {application} -parallel | tee {log_name}",
            timeout_seconds,
        )
        base.run_cmd(["reconstructPar"], description="Reconstructing parallel domain fields...")
    else:
        base.check_tool(application)
        _run_timed(f"{application} | tee {log_name}", timeout_seconds)


def _write_csv(csv_path, rows):
    fieldnames = ["AoA_deg", "Cl", "Cd", "Cl_over_Cd", "Cl_std", "Cd_std", "mode", "converged"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _nice_ticks(lo, hi, n=6):
    """Round tick positions covering [lo, hi] (used by the dependency-free SVG plotter)."""
    if hi <= lo:
        hi = lo + 1.0
    raw = (hi - lo) / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(raw))
    step = 10 * mag
    for m in (1, 2, 2.5, 5, 10):
        if m * mag >= raw:
            step = m * mag
            break
    first = math.floor(lo / step) * step
    last = math.ceil(hi / step) * step
    count = int(round((last - first) / step))
    return [round(first + i * step, 12) for i in range(count + 1)]


def _svg_panel(x0, y0, w, h, xs, ys, xlabel, ylabel, title, color, bad_idx):
    ml, mr, mt, mb = 66, 18, 36, 48
    px, py, pw, ph = x0 + ml, y0 + mt, w - ml - mr, h - mt - mb
    xt, yt = _nice_ticks(min(xs), max(xs)), _nice_ticks(min(ys), max(ys))
    xmin, xmax, ymin, ymax = xt[0], xt[-1], yt[0], yt[-1]
    sx = lambda v: px + (v - xmin) / (xmax - xmin) * pw
    sy = lambda v: py + ph - (v - ymin) / (ymax - ymin) * ph
    o = [f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" fill="white" stroke="#444"/>']
    for t in xt:
        X = sx(t)
        o.append(f'<line x1="{X:.1f}" y1="{py}" x2="{X:.1f}" y2="{py + ph}" stroke="#ddd"/>')
        o.append(f'<text x="{X:.1f}" y="{py + ph + 16}" font-size="11" text-anchor="middle">{t:g}</text>')
    for t in yt:
        Y = sy(t)
        o.append(f'<line x1="{px}" y1="{Y:.1f}" x2="{px + pw}" y2="{Y:.1f}" stroke="#ddd"/>')
        o.append(f'<text x="{px - 6}" y="{Y + 4:.1f}" font-size="11" text-anchor="end">{t:g}</text>')
    pts = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(xs, ys))
    o.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"/>')
    for i, (x, y) in enumerate(zip(xs, ys)):
        o.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.5" fill="{color}"/>')
        if i in bad_idx:
            o.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="8" fill="none" stroke="red" stroke-width="1.5"/>')
    o.append(f'<text x="{px + pw / 2}" y="{y0 + 22}" font-size="14" text-anchor="middle" font-weight="bold">{title}</text>')
    o.append(f'<text x="{px + pw / 2}" y="{y0 + h - 8}" font-size="12" text-anchor="middle">{xlabel}</text>')
    o.append(f'<text x="{x0 + 16}" y="{py + ph / 2}" font-size="12" text-anchor="middle" '
             f'transform="rotate(-90 {x0 + 16} {py + ph / 2})">{ylabel}</text>')
    return "\n".join(o)


def _write_svg_plots(rows, path):
    """Dependency-free fallback: 4-panel SVG (opens in any browser)."""
    a = [r[0] for r in rows]
    cl = [r[1] for r in rows]
    cd = [r[2] for r in rows]
    ld = [c / d if d else 0.0 for c, d in zip(cl, cd)]
    bad = {i for i, r in enumerate(rows) if not r[3]}
    W, H = 980, 720
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
             f'font-family="Arial, Helvetica, sans-serif">',
             f'<rect width="{W}" height="{H}" fill="white"/>',
             f'<text x="{W / 2}" y="26" font-size="17" text-anchor="middle" font-weight="bold">AoA sweep results</text>',
             _svg_panel(10, 40, 470, 330, a, cl, "Angle of attack (deg)", "Cl", "Lift coefficient vs AoA", "#1f77b4", bad),
             _svg_panel(500, 40, 470, 330, a, cd, "Angle of attack (deg)", "Cd", "Drag coefficient vs AoA", "#ff7f0e", bad),
             _svg_panel(10, 380, 470, 330, cd, cl, "Cd", "Cl", "Drag polar", "#2ca02c", bad),
             _svg_panel(500, 380, 470, 330, a, ld, "Angle of attack (deg)", "Cl / Cd", "Lift-to-drag ratio vs AoA", "#d62728", bad),
             "</svg>"]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))


def _maybe_plot_polar(csv_path, results_dir):
    """Called automatically at the end of every sweep. Draws Cl-AoA, Cd-AoA, drag polar and
    Cl/Cd-AoA from the results CSV. Uses matplotlib (PNG) if installed, otherwise falls back
    to a dependency-free SVG, so a plot is always produced. Never raises."""
    try:
        rows = []
        with open(csv_path, "r") as f:
            for row in csv.DictReader(f):
                if row["Cl"] == "" or row["Cd"] == "":
                    continue
                rows.append((float(row["AoA_deg"]), float(row["Cl"]), float(row["Cd"]),
                             str(row.get("converged", "True")).strip().lower() == "true"))
        if not rows:
            return
        rows.sort()
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            svg_path = os.path.join(results_dir, "sweep_plots.svg")
            _write_svg_plots(rows, svg_path)
            print("[Note] matplotlib not available - wrote dependency-free SVG plots instead "
                  "(open in a browser):")
            print(f"       {svg_path}")
            return

        a = [r[0] for r in rows]
        cl = [r[1] for r in rows]
        cd = [r[2] for r in rows]
        ld = [c / d if d else float("nan") for c, d in zip(cl, cd)]
        bad = [i for i, r in enumerate(rows) if not r[3]]
        fig, ax = plt.subplots(2, 2, figsize=(11, 8))

        def draw(axis, x, y, xl, yl, title, color):
            axis.plot(x, y, marker="o", color=color)
            if bad:
                axis.plot([x[i] for i in bad], [y[i] for i in bad], "o", mfc="none", mec="red",
                          ms=11, label="not converged")
                axis.legend()
            axis.set_xlabel(xl)
            axis.set_ylabel(yl)
            axis.set_title(title)
            axis.grid(True, alpha=0.3)

        draw(ax[0][0], a, cl, "Angle of attack (deg)", "Cl", "Lift coefficient vs AoA", "tab:blue")
        draw(ax[0][1], a, cd, "Angle of attack (deg)", "Cd", "Drag coefficient vs AoA", "tab:orange")
        draw(ax[1][0], cd, cl, "Cd", "Cl", "Drag polar", "tab:green")
        draw(ax[1][1], a, ld, "Angle of attack (deg)", "Cl / Cd", "Lift-to-drag ratio vs AoA", "tab:red")
        fig.suptitle("AoA sweep results", fontsize=13)
        fig.tight_layout()
        plot_path = os.path.join(results_dir, "sweep_plots.png")
        fig.savefig(plot_path, dpi=150)
        print(f"Plots saved to: {plot_path}")
    except Exception as exc:  # plotting must never break a finished sweep
        print(f"[Note] could not draw plots ({exc}); the CSV is still complete.")

# --------------------------------------------------------------------------
# Main sweep
# --------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("  2D Airfoil Angle-of-Attack Sweep")
    print("=" * 65)

    case_dir = os.getcwd()
    base.check_tool("fluentMeshToFoam")

    # Must clear out EVERYTHING from any previous run — including
    # postProcessing/ — before starting. Otherwise a stale forceCoeffs
    # directory from a prior run's last (highest-time) angle can be mistaken
    # for "the latest result" of this run's very first angle, silently
    # reporting old data instead of a fresh computation.
    base.clean_previous_run(case_dir)

    mesh_file = base.find_msh_file(case_dir)
    if not mesh_file:
        mesh_file = input("No .msh file detected. Enter filename manually: ").strip()
        if not os.path.exists(mesh_file):
            print(f"[ERROR] Specified mesh file '{mesh_file}' does not exist.")
            sys.exit(1)
    print(f"\n[Selected Mesh]: {mesh_file} (airfoil at AoA = 0, single mesh for the whole sweep)")

    turbulent = base.confirm("\nUse RAS turbulence (recommended for airfoil cases)?", default_yes=True)
    ras_model = None
    if turbulent:
        ras_choice = base.get_choice("RAS turbulence model:", ["k-omega SST", "k-epsilon"], default_index=0)
        ras_model = "kOmegaSST" if ras_choice == "k-omega SST" else "kEpsilon"
    second_name = "epsilon" if ras_model == "kEpsilon" else "omega"

    base.generate_control_dict(case_dir, application="simpleFoam", steady=True, end_time=1000)
    base.generate_fv_schemes(case_dir, "simpleFoam", turbulent, ras_model)
    base.generate_fv_solution(case_dir, "simpleFoam", turbulent, ras_model)
    base.generate_turbulence_properties(case_dir, turbulent, ras_model)

    print("\n[1/3] Converting ANSYS Fluent mesh using fluentMeshToFoam...")
    base.run_cmd(["fluentMeshToFoam", mesh_file])
    print("Mesh conversion completed successfully.")

    # Auto-detect whether the airfoil chord needs rotating onto +X with the
    # leading edge upstream, instead of assuming one fixed mesh orientation —
    # different meshes (or the same mesh rebuilt differently in ANSYS) can
    # come in already correct, or needing 90, -90, or 180 degrees.
    required_rotation = detect_required_rotation(mesh_file)
    if required_rotation is None:
        print("\n[Warning] Could not auto-detect the airfoil's chord orientation from the "
              "mesh — skipping rotation. If the flow ends up hitting the wrong face of the "
              "airfoil, check its orientation manually.")
    elif required_rotation == 0:
        print("\n[Info] Airfoil chord is already correctly aligned along +X "
              "(leading edge upstream) — no rotation needed.")
    else:
        print(f"\n[Info] Rotating mesh {required_rotation} deg around Z so the chord aligns "
              f"with +X, leading edge upstream (detected from the actual airfoil geometry).")
        base.check_tool("transformPoints")
        base.run_cmd(["transformPoints", "-rotate-z", str(required_rotation)],
                      description=f"Rotating mesh {required_rotation} degrees...")

    scale = base.get_input(
        "\nMesh length-unit scale factor (e.g. 0.001 if ANSYS mesh was in mm, 1 for no change)",
        1.0, float, min_val=0.0,
    )
    if scale != 1.0:
        base.check_tool("transformPoints")
        base.run_cmd(["transformPoints", "-scale", f"({scale} {scale} {scale})"],
                      description="Scaling mesh points...")

    # checkMesh is mandatory (not skippable) — its output drives the solver
    # settings chosen below, so every mesh gets settings matched to its own
    # actual quality instead of one fixed assumption.
    base.check_tool("checkMesh")
    checkmesh_info = base.run_checkmesh_and_parse(case_dir)
    detected_z_span = checkmesh_info.get("z_span")

    settings = base.recommend_solver_settings(checkmesh_info)
    print("\n[Auto-tuned solver settings from checkMesh]")
    print(f"  max non-orthogonality : {checkmesh_info.get('max_non_orthogonality', 'n/a')}")
    print(f"  max aspect ratio      : {checkmesh_info.get('max_aspect_ratio', 'n/a')}")
    print(f"  max skewness          : {checkmesh_info.get('max_skewness', 'n/a')}")
    print(f"  -> quality tier       : {settings['quality_tier']}")
    print(f"  -> non-orth scheme    : {settings['non_orth_scheme']} "
          f"({settings['n_non_orth_correctors']} correctors)")
    print(f"  -> relaxation (p/U/turb): {settings['relax_p']} / {settings['relax_U']} / "
          f"{settings['relax_turb']}")
    print(f"  -> algorithm          : {'SIMPLEC' if settings['use_simplec'] else 'SIMPLE'}")

    print("\n[2/3] Inspecting constant/polyMesh/boundary patches...")
    patches = base.parse_boundary_patches(case_dir)
    classification = classify_airfoil_patches(patches)
    print(f"Found {len(patches)} boundary patches:")
    for name, ptype in patches:
        print(f" - {name:<20} (mesh type: {ptype:<10}) -> classified as: {classification[name]}")
    if not base.confirm("\nProceed with this classification (farfield / wall / empty)?", default_yes=True):
        print("Aborted by user. Rename patches in ANSYS so one clearly reads as the")
        print("farfield boundary and one as the airfoil wall, then re-run.")
        sys.exit(0)

    if fix_farfield_patch_type(case_dir, classification):
        print("(If you'd rather this be correct at the source, set that zone's type to")
        print(" pressure-far-field in ANSYS and re-export — this fix is applied either way.)")

    print("\n[3/3] Flow parameters:")
    u_inf = base.get_input("Freestream velocity U_inf (m/s)", 20.0, float, min_val=0.0)
    nu = base.get_input("Kinematic viscosity nu (m^2/s)", 1.5e-5, float, min_val=0.0)

    base.generate_transport_properties(case_dir, nu)

    chord = base.get_input("Reference chord length lRef (m)", 1.0, float, min_val=1e-6)
    if detected_z_span is not None and detected_z_span > 0:
        aref = chord * detected_z_span
        print(f"[Info] Aref = chord x real mesh Z-span = {chord:g} x {detected_z_span:.6g} "
              f"= {aref:.6g} m^2 (read from checkMesh, not assumed).")
    else:
        aref = chord * 1.0
        print("[Note] Could not read the mesh's real Z-span from checkMesh — assuming a unit "
              "span (Aref = chord x 1). If checkMesh was skipped, re-run with it enabled for "
              "correctly-scaled force coefficients.")

    k_val = second_val = None
    if turbulent:
        turb_intensity = base.get_input("Freestream turbulence intensity I (%)", 0.5, float, min_val=0.0)
        length_scale = base.get_input("Turbulent length scale (m)", 0.01 * chord, float, min_val=1e-9)
        k_val, second_val = base.estimate_turbulence_quantities(u_inf, turb_intensity, length_scale, ras_model)
        print(f"  -> Estimated freestream k = {k_val:.4g} m^2/s^2, {second_name} = {second_val:.4g}")

    aoa_start = base.get_input("AoA sweep start (deg)", -4.0, float)
    aoa_end = base.get_input("AoA sweep end (deg)", 16.0, float)
    aoa_step = base.get_input("AoA sweep step (deg)", 1.0, float, min_val=0.01)

    n_iters_steady = base.get_input(
        "Steady (SIMPLE) iterations per angle", 1500, int, min_val=1
    )

    run_timeout_min = base.get_input(
        "Safety time limit per angle, in minutes (0 = no limit) — protects your "
        "overall time budget if a mesh/angle turns out slower than expected", 8.0, float, min_val=0.0
    )
    run_timeout_sec = run_timeout_min * 60 if run_timeout_min > 0 else None

    n_cores = base.get_input("Number of CPU cores per run (1 for serial)", 1, int, min_val=1)

    angles = []
    a = aoa_start
    while a <= aoa_end + 1e-9:
        angles.append(round(a, 6))
        a += aoa_step

    print(f"\nSweep plan: {len(angles)} angles from {aoa_start} to {aoa_end} deg, step {aoa_step} deg.")
    if not base.confirm("Proceed with the sweep?", default_yes=True):
        print("Aborted by user before running.")
        sys.exit(0)

    results_dir = os.path.join(case_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, "aoa_sweep_results.csv")
    csv_rows = []

    # endTime in controlDict is an ABSOLUTE time, not "N more iterations". If a
    # non-converged angle runs all the way to the ceiling, its latest time dir
    # IS that ceiling value — so the next warm-started angle must be given a
    # NEW, higher endTime (current_time + n_iters_steady), or it starts already
    # "at" endTime and the solver does zero iterations, silently freezing every
    # angle after the first one that fails to converge.
    current_time = 0.0

    for i, aoa in enumerate(angles):
        print("\n" + "-" * 65)
        print(f" Angle {i + 1}/{len(angles)}: AoA = {aoa} deg")
        print("-" * 65)

        ux, uy = velocity_components(u_inf, aoa)
        lift_dir, drag_dir = lift_drag_dirs(aoa)
        wall_patches = [name for name, kind in classification.items() if kind == "wall"]
        force_block_raw = base.build_force_coeffs_block(
            wall_patches, u_inf, chord, aref, lift_dir, drag_dir,
        )
        yplus_block_raw = base.build_yplus_block(wall_patches) if turbulent else ""
        # Both helpers each wrap their own top-level 'functions{}' — only one
        # such block is valid per file, so merge their inner entries here
        # rather than concatenating the two wrapped blocks directly.
        inner_entries = []
        for raw in (force_block_raw, yplus_block_raw):
            if raw.strip():
                inner = raw.strip()
                inner = inner[inner.index("{") + 1: inner.rindex("}")].strip()
                inner_entries.append(inner)
        functions_block = (
            "\nfunctions\n{\n" + "\n\n".join(inner_entries) + "\n}\n"
            if inner_entries else ""
        )

        is_first = (i == 0)
        if is_first:
            write_all_fields(case_dir, classification, ux, uy, k_val, second_name, second_val, turbulent)
            start_from = "startTime"
        else:
            warm_start_boundary_conditions(
                case_dir, classification, ux, uy, second_name, turbulent,
                k_val, second_val
            )
            start_from = "latestTime"

        target_end_time = current_time + n_iters_steady

        base.generate_control_dict(
            case_dir,
            application="simpleFoam",
            steady=True,
            end_time=target_end_time,
            write_interval=max(1, n_iters_steady // 10),
            functions_block=functions_block,
            start_from=start_from,
            purge_write=2,
        )
        base.generate_fv_schemes(case_dir, "simpleFoam", turbulent, ras_model,
                                  non_orth_scheme=settings["non_orth_scheme"])
        base.generate_fv_solution(case_dir, "simpleFoam", turbulent, ras_model,
                                   n_non_orth_correctors=settings["n_non_orth_correctors"],
                                   relax_p=settings["relax_p"], relax_U=settings["relax_U"],
                                   relax_turb=settings["relax_turb"],
                                   use_simplec=settings["use_simplec"])

        log_name = f"log.simpleFoam.aoa_{aoa}"
        run_solver(case_dir, "simpleFoam", n_cores, log_name, timeout_seconds=run_timeout_sec)

        # Record the ACTUAL time reached (converged early, or the ceiling) so
        # the next angle's target_end_time is computed from reality, not from
        # what we asked for — these can differ (early convergence stops short).
        latest_dir = get_latest_time_dir(case_dir)
        if latest_dir is not None:
            try:
                current_time = float(os.path.basename(latest_dir))
            except ValueError:
                current_time = target_end_time
        else:
            current_time = target_end_time

        simple_converged = check_simple_converged(os.path.join(case_dir, log_name))
        force_converged = check_force_converged(case_dir)
        converged = simple_converged or force_converged
        print(f"[Convergence] residual-based: {simple_converged} | "
              f"force-based (Cd/Cl stable over last 100 its): {force_converged}")
        mode = "steady"
        force_data = read_force_coeffs(case_dir, average_last_fraction=0.1)

        angle_dir = os.path.join(results_dir, f"aoa_{aoa}")
        os.makedirs(angle_dir, exist_ok=True)
        log_src = os.path.join(case_dir, log_name)
        if os.path.exists(log_src):
            shutil.copy(log_src, angle_dir)
        fc_file = find_force_coeffs_file(case_dir)
        if fc_file:
            shutil.copy(fc_file, angle_dir)

        row = {
            "AoA_deg": aoa,
            "Cl": force_data["Cl"] if force_data else "",
            "Cd": force_data["Cd"] if force_data else "",
            "Cl_over_Cd": (force_data["Cl"] / force_data["Cd"])
                          if force_data and force_data["Cd"] not in (0, None) else "",
            "Cl_std": force_data["Cl_std"] if force_data else "",
            "Cd_std": force_data["Cd_std"] if force_data else "",
            "mode": mode,
            "converged": converged,
        }
        csv_rows.append(row)
        _write_csv(csv_path, csv_rows)

        if force_data:
            print(f"[Result] AoA={aoa} deg | mode={mode} | Cl={force_data['Cl']:.4f} "
                  f"Cd={force_data['Cd']:.5f}")
        else:
            print(f"[Result] AoA={aoa} deg | mode={mode} | FAILED to extract force coefficients.")

    print(f"\nSweep complete. Results CSV: {csv_path}")
    _maybe_plot_polar(csv_path, results_dir)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Exiting.")
        sys.exit(130)