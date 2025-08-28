# main.py — FactSage XML → CSV analyzer for Pyodide
# Preferred entrypoint for JS:
#   analyze_xml_b64_to_file(b64data: str, outfile: str = "analysis.csv") -> str (returns outfile path)
# Back-compat (not recommended for very large outputs):
#   analyze_xml_b64(b64data: str) -> str (returns CSV text as a big string)

import re
import base64
import io
from xml.etree import ElementTree as ET

# ============== Formatting & helpers ==============

SUB_MAP = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

def subscript_formula_digits(text: str) -> str:
    """Convert digits in chemical formulas to unicode subscripts (cosmetic only)."""
    return re.sub(r'(?<=[A-Za-z\)])\d+', lambda m: m.group(0).translate(SUB_MAP), text or "")

def clean_phase_state(raw_state: str) -> str:
    """Remove FTxxx- prefixes from phase names/states."""
    return re.sub(r'^FT[a-zA-Z]+-', '', raw_state or "")

def fnum(x) -> float:
    """Safe float parse. Accepts '', None, '1.23D+02', etc."""
    if x is None:
        return 0.0
    s = str(x).strip().replace('D', 'E')
    if s == "":
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0

def fmt3(x: float) -> str:
    return "0" if abs(x) == 0 else f"{x:.3f}"

def fmt4(x: float) -> str:
    return "0" if abs(x) == 0 else f"{x:.4f}"

G_MIN = 0.01  # grams threshold to list species/solids


# ============== XML → TXT (FactSage-like layout) ==============

def convert_factsage_xml_to_txt(xml_path: str, out_path: str):
    """
    Convert a FactSage XML export to a compact TXT layout we can re-parse easily.
    Keeps:
      - Page number, T (°C), P (atm)
      - Page-scope Reactants (moles & grams)
      - Solutions (PHASE blocks): species rows with g, W, n, X, a (filtered by G_MIN)
      - PURE SOLIDS with g & a (filtered by G_MIN)
    """
    root = ET.parse(xml_path).getroot()

    # Species definition: map solution phase -> species IDs/names
    spec_def = root.find("./header/species_definition")

    solution_name_map = {}
    phase_species_ids = {}
    species_in_solutions = {}

    if spec_def is not None:
        for sol in spec_def.findall("solution"):
            pid = sol.attrib.get("phase_id")
            solution_name_map[pid] = clean_phase_state(sol.attrib.get("state", ""))
            ids = []
            for sp in sol.findall("species"):
                sid = sp.attrib["id"]
                species_in_solutions[sid] = {"name": sp.attrib.get("name", ""), "phase_id": pid}
                ids.append(sid)
            phase_species_ids[pid] = ids

    # Pure solid species mapping (OPTIONAL)
    solid_species_map = {
        sp.attrib["id"]: sp.attrib
        for sp in root.iter("species")
        if sp.attrib.get("phase") == "s"
    }

    # Reactants meta (name, MW)
    reactants_meta = {
        r.attrib["id"]: {"name": r.attrib.get("name", ""), "mw": fnum(r.attrib.get("mw"))}
        for r in root.findall("./header/reactant")
    }

    formula_pretty = subscript_formula_digits(root.attrib.get("formula", ""))

    out_lines = []

    for page in root.findall("./page"):
        pid = page.attrib.get("id", "?")
        desc = page.attrib.get("description", "").strip()

        # Temperature
        t_raw = (page.attrib.get("T") or "").strip()
        if t_raw:
            tK = fnum(t_raw)
            tC = tK - 273.15
        else:
            m = re.search(r'(-?\d+(?:\.\d+)?)\s*C', desc)
            if m:
                tC = fnum(m.group(1))
                tK = tC + 273.15
            else:
                tC, tK = 0.0, 273.15

        # Pressure
        p_raw = (page.attrib.get("P") or "").strip()
        pbar = fnum(p_raw) if p_raw else 1.01325
        patm = pbar / 1.01325

        out_lines.append("=" * 70)
        out_lines.append(f"Page {pid} - {fmt3(tC)} C")
        out_lines.append(formula_pretty)
        out_lines.append(f"T = {fmt3(tC)} C | T(K) = {fmt3(tK)} | P = {fmt3(patm)} atm")
        out_lines.append("")
        out_lines.append("Reactants (page scope):")

        # Page-scope reactants instances
        for rid, meta in reactants_meta.items():
            rnode = page.find(f"./reactant[@id='{rid}']")
            if rnode is None:
                continue
            n = fnum(rnode.attrib.get("n"))
            g = n * meta["mw"]
            out_lines.append(f"  {meta['name']:<4} n = {fmt3(n)} mol | m = {fmt3(g)} g")

        out_lines.append("")

        # Results per species (species totals across solutions)
        results = {r.attrib["id"]: r.attrib for r in page.findall("result")}

        # Solution totals to list PHASE sections in descending g
        phase_mass = {s.attrib["id"]: fnum(s.attrib.get("g")) for s in page.findall("solution")}
        ordered = sorted([ph for ph, g in phase_mass.items() if g >= G_MIN],
                         key=lambda x: phase_mass[x], reverse=True)

        for ph_id in ordered:
            pname = solution_name_map.get(ph_id, f"Phase-{ph_id}")
            out_lines.append(f"PHASE: {pname}")
            out_lines.append(f"{'Compound':<16}{'Mass(g)':>10}{'W(%)':>10}{'Mol':>8}{'X':>8}{'Activity':>12}")

            total_g = total_n = total_W = total_X = total_a = 0.0
            rows = []

            for sid in phase_species_ids.get(ph_id, []):
                sp = species_in_solutions.get(sid)
                if not sp:
                    continue
                r = results.get(sid, {})
                g = fnum(r.get("g"))
                W = fnum(r.get("W"))
                n = fnum(r.get("n"))
                X = fnum(r.get("X"))
                a = fnum(r.get("a"))

                if g < G_MIN:
                    continue

                rows.append({"name": sp["name"], "g": g, "W": W, "n": n, "X": X, "a": a})
                total_g += g
                total_n += n
                total_W += W
                total_X += X
                total_a += a

            rows.sort(key=lambda rr: rr["g"], reverse=True)
            for row in rows:
                out_lines.append(
                    f"  {row['name']:<14}{fmt3(row['g']):>10}"
                    f"{fmt3(row['W'] * 100):>10}{fmt3(row['n']):>8}"
                    f"{fmt3(row['X']):>8}{fmt4(row['a']):>12}"
                )

            out_lines.append(
                f"  {'TOTAL':<14}{fmt3(total_g):>10}{fmt3(100.00):>10}"
                f"{fmt3(total_n):>8}{fmt3(1.00):>8}{fmt4(total_a):>12}"
            )
            out_lines.append("")

        # PURE SOLIDS
        solid_rows = []
        for sid, sp in solid_species_map.items():
            r = results.get(sid)
            if not r:
                continue
            g = fnum(r.get("g"))
            a = fnum(r.get("a"))
            if g >= G_MIN:
                solid_rows.append({"name": sp.get("name", ""), "g": g, "a": a})

        if solid_rows:
            out_lines.append("PURE SOLIDS:")
            out_lines.append(f"{'Compound':<24}{'g':>12}{'Activity':>14}")
            for row in sorted(solid_rows, key=lambda x: x["g"], reverse=True):
                out_lines.append(f"  {row['name']:<22}{fmt3(row['g']):>12}{fmt4(row['a']):>14}")
            out_lines.append("")

    # Stream writing to avoid building a single huge string in memory
    with open(out_path, "w", encoding="utf-8") as f:
        for line in out_lines:
            f.write(line)
            f.write("\n")


# ============== TXT Parsing (from our layout) ==============

def parse_txt_pages(txt_path: str):
    with open(txt_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    pages = []
    i = 0
    L = len(lines)

    def next_is_page_sep(k: int) -> bool:
        if k + 1 < L:
            return (lines[k].startswith("=" * 10)) and ("Page" in lines[k + 1])
        return False

    while i < L:
        line = lines[i]
        m = re.match(r'^Page\s+(\d+)\s*-\s*([0-9.]+)\s*C', line.strip())
        if not m:
            i += 1
            continue

        page_num = int(m.group(1))
        temp_C = float(m.group(2))
        fe_mass_g = None

        # Read Fe reactant mass if present
        j = i + 1
        while j < L and not lines[j].startswith("PHASE:"):
            if lines[j].strip().startswith("Reactants (page scope):"):
                j += 1
                while j < L and lines[j].strip() and not lines[j].startswith("PHASE:"):
                    mfe = re.search(r'^\s*Fe\b.*?m\s*=\s*([0-9.]+)\s*g', lines[j])
                    if mfe:
                        try:
                            fe_mass_g = float(mfe.group(1))
                        except Exception:
                            pass
                    j += 1
                break
            j += 1

        phases = []
        pure_solids = []

        k = j
        while k < L and not next_is_page_sep(k):
            st = lines[k].strip()

            if st.startswith("PHASE:"):
                ph_name = st.split("PHASE:")[1].strip()
                k += 2  # skip header line too
                rows = []
                while k < L:
                    s2 = lines[k].strip()
                    if (not s2) or s2.startswith("PHASE:") or s2.startswith("PURE SOLIDS:") or s2.startswith("System Thermodynamics:"):
                        break
                    if re.match(r'^\s*TOTAL', s2):
                        k += 1
                        continue
                    mrow = re.match(
                        r'^\s*(\S+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s*$',
                        s2
                    )
                    if mrow:
                        rows.append({
                            "name": mrow.group(1),
                            "g": float(mrow.group(2)),
                            "W": float(mrow.group(3)),
                            "n": float(mrow.group(4)),
                            "X": float(mrow.group(5)),
                            "a": float(mrow.group(6)),
                        })
                    k += 1
                phases.append({"name": ph_name, "rows": rows})
                continue

            if st.startswith("PURE SOLIDS:"):
                k += 2  # skip header lines
                while k < L:
                    s2 = lines[k].strip()
                    if (not s2) or s2.startswith("PHASE:") or s2.startswith("System Thermodynamics:"):
                        break
                    mps = re.match(r'^\s*(\S+)\s+([0-9.\-]+)\s+([0-9.\-]+)\s*$', s2)
                    if mps:
                        pure_solids.append({"name": mps.group(1), "g": float(mps.group(2))})
                    k += 1
                continue

            k += 1

        pages.append({
            "page_num": page_num,
            "T_C": temp_C,
            "Fe_g": fe_mass_g,
            "phases": phases,
            "pure_solids": pure_solids
        })
        i = k

    return pages


# ============== Phase classification helpers ==============

def is_metal_liquid_phase(name: str) -> bool:
    n = (name or "").lower()
    return ("fe-liq" in n) or ("liquid" in n)

def is_slag_liquid_phase(name: str) -> bool:
    n = (name or "").lower()
    return ("slag" in n) and ("liq" in n or "liquid" in n)

def is_metal_intermetallic_phase(name: str) -> bool:
    n = (name or "").lower()
    if "liq" in n or "liquid" in n or "slag" in n:
        return False
    return bool(re.search(r'\b(bcc|fcc|hcp)\b', n)) or ("fe" in n and "si" in n and "o" not in n)

def is_metal_intermetallic_puresolid(name: str) -> bool:
    raw = (name or "").replace("(s)", "")
    return ("Fe" in raw and "Si" in raw and "O" not in raw)

def is_slag_solid_phase(name: str) -> bool:
    n = (name or "").lower()
    if "liq" in n or "liquid" in n:
        return False
    if "slag" in n:
        return True
    return "o" in n  # heuristic: oxide solids

def is_slag_puresolid(name: str) -> bool:
    return "O" in (name or "")


# ============== Composition formatting ==============

def format_comp(comp_dict: dict) -> str:
    """Return composition string WITH wt% + activity (no total)."""
    def get_w(x):
        return x[0] if isinstance(x, (tuple, list)) and len(x) >= 1 else x

    items = sorted(comp_dict.items(), key=lambda kv: get_w(kv[1]), reverse=True)
    parts = []
    for name, val in items:
        if isinstance(val, (tuple, list)):
            W = val[0] if len(val) >= 1 else 0.0
            a = val[1] if len(val) >= 2 else 0.0
            parts.append(f"{name}:  {fmt3(W)}  (a={fmt4(a)})  ")
        else:
            parts.append(f"{name}:  {fmt3(float(val))}  ")
    return "  |  ".join(parts)

def comp_total_mass(comp_dict: dict) -> float:
    total = 0.0
    for v in comp_dict.values():
        if isinstance(v, (tuple, list)) and len(v) >= 3:
            total += float(v[2])
    return total

def comp_text_with_total(comp_dict: dict) -> str:
    base = format_comp(comp_dict)
    tot = comp_total_mass(comp_dict)
    return f"{base}  |  TOTAL: {fmt3(tot)} g"


# ============== Group pages into simulations ==============

def group_simulations(pages: list, fe_tol: float = 1e-3):
    sims = []
    cur = []
    last_fe = None
    for p in pages:
        reset = (p.get("page_num") == 1 and bool(cur)) or (
            last_fe is not None and p.get("Fe_g") is not None and abs(p["Fe_g"] - last_fe) > fe_tol
        )
        if reset:
            sims.append(cur)
            cur = []
        cur.append(p)
        if p.get("Fe_g") is not None:
            last_fe = p["Fe_g"]
    if cur:
        sims.append(cur)
    return sims


# ============== Analyze a simulation ==============

def analyze_simulation(sim_pages: list, g_threshold: float = G_MIN):
    # Work from high T to low T
    sim_pages = sorted(sim_pages, key=lambda x: x["T_C"], reverse=True)
    fe_g = sim_pages[0].get("Fe_g") if sim_pages else None

    collected = []  # Liquid#1 snapshots with Si wt%
    stop_T = None
    stop_phases = []

    prev_slag_comp = None
    slag_first_T = None
    slag_first_phases = []
    slag_comp_before = None

    for pg in sim_pages:
        # Slag solids first appearance
        if slag_first_T is None:
            slag_phases_here = []
            for ph in pg.get("phases", []):
                if ph.get("rows") and is_slag_solid_phase(ph.get("name", "")):
                    slag_phases_here.append(ph["name"])
            for ps in pg.get("pure_solids", []):
                if ps.get("g", 0) >= g_threshold and is_slag_puresolid(ps.get("name", "")):
                    slag_phases_here.append(ps["name"])
            if slag_phases_here:
                slag_first_T = pg["T_C"]
                slag_first_phases = sorted(set(slag_phases_here))
                slag_comp_before = prev_slag_comp

        # Metal precipitates (intermetallics)
        metal_phases_here = []
        for ph in pg.get("phases", []):
            if ph.get("rows") and is_metal_intermetallic_phase(ph.get("name", "")):
                metal_phases_here.append(ph["name"])
        for ps in pg.get("pure_solids", []):
            if ps.get("g", 0) >= g_threshold and is_metal_intermetallic_puresolid(ps.get("name", "")):
                metal_phases_here.append(ps["name"])
        if metal_phases_here and stop_T is None:
            stop_T = pg["T_C"]
            stop_phases = sorted(set(metal_phases_here))
            break  # stop at first precipitation

        # Liquid#1 (metallic liquid)
        liq = None
        for ph in pg.get("phases", []):
            if is_metal_liquid_phase(ph.get("name", "")) and ph.get("rows"):
                liq = ph
                break
        if liq:
            comp = {r["name"]: (r["W"], r["a"], r["g"]) for r in liq.get("rows", [])}
            si_w = comp.get("Si", (0.0, 0.0, 0.0))[0]
            collected.append({"T_C": pg["T_C"], "comp": comp, "Si_w": si_w})

        # Slag liquid comp (for "before first solids")
        slag_liq = None
        for ph in pg.get("phases", []):
            if is_slag_liquid_phase(ph.get("name", "")) and ph.get("rows"):
                slag_liq = ph
                break
        if slag_liq:
            prev_slag_comp = {r["name"]: (r["W"], r["a"], r["g"]) for r in slag_liq.get("rows", [])}

    if not collected:
        return {
            "Fe_g": fe_g,
            "best_T": None,
            "best_Si": None,
            "best_comp": {},
            "stop_T": stop_T,
            "stop_phases": stop_phases,
            "slag_first_T": slag_first_T,
            "slag_first_phases": slag_first_phases,
            "slag_comp_before": slag_comp_before or {},
        }

    best = max(collected, key=lambda r: r["Si_w"])
    return {
        "Fe_g": fe_g,
        "best_T": best["T_C"],
        "best_Si": best["Si_w"],
        "best_comp": best["comp"],
        "stop_T": stop_T,
        "stop_phases": stop_phases,
        "slag_first_T": slag_first_T,
        "slag_first_phases": slag_first_phases,
        "slag_comp_before": slag_comp_before or {},
    }


def run_liquid_analysis(txt_path: str, top_k: int = 3):
    pages = parse_txt_pages(txt_path)
    sims = group_simulations(pages)
    results = [analyze_simulation(sim) for sim in sims]
    # keep only valid runs with Liquid#1
    results = [r for r in results if r["best_Si"] is not None and r["best_comp"]]
    results.sort(key=lambda r: -(r["best_Si"]))
    return results[:top_k]


# ============== Public API for app.js (file-based) ==============

def analyze_xml_b64_to_file(b64data: str, outfile: str = "analysis.csv") -> str:
    """
    Preferred entry from JS. Writes CSV to `outfile` in Pyodide FS and returns the path.
    This avoids returning huge strings through the Pyodide bridge (prevents JS call stack errors).
    """
    xml_bytes = base64.b64decode(b64data)
    xml_text = xml_bytes.decode("utf-8", "ignore")

    xml_filename = "input.xml"
    txt_out = "input_parsed.txt"
    with open(xml_filename, "w", encoding="utf-8") as f:
        f.write(xml_text)

    convert_factsage_xml_to_txt(xml_filename, txt_out)
    hits = run_liquid_analysis(txt_out, top_k=3)

    import csv
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank",
            "Fe mass (g)",
            "Best T (°C)",
            "Si wt% (max)",
            "Liquid#1 composition",
            "Stopped at T",
            "First precipitates",
            "Slag first T",
            "Slag phases",
            "Slag liquid (before first solids)"
        ])
        for idx, r in enumerate(hits, 1):
            writer.writerow([
                idx,
                "not found" if r["Fe_g"] is None else fmt3(r["Fe_g"]),
                "not found" if r["best_T"] is None else f"{r['best_T']:.2f}",
                "not found" if r["best_Si"] is None else fmt3(r["best_Si"]),
                comp_text_with_total(r["best_comp"]) if r["best_comp"] else "not found",
                "not found" if r["stop_T"] is None else f"{r['stop_T']:.2f}",
                ", ".join(r["stop_phases"]) if r["stop_phases"] else "not found",
                "not found" if r["slag_first_T"] is None else f"{r['slag_first_T']:.2f}",
                ", ".join(r["slag_first_phases"]) if r["slag_first_phases"] else "not found",
                comp_text_with_total(r["slag_comp_before"]) if r["slag_comp_before"] else "not found"
            ])

    return outfile


# ============== Back-compat (string return, not recommended for huge CSV) ==============

def analyze_xml_b64(b64data: str) -> str:
    """
    Kept for backward compatibility. Builds the same CSV as a big string.
    For large outputs, prefer analyze_xml_b64_to_file and read the file from JS.
    """
    out_path = analyze_xml_b64_to_file(b64data, outfile="analysis.csv")
    # Read the just-written file back (this may be large!)
    with open(out_path, "r", encoding="utf-8") as f:
        return f.read()

def analyze_xml_path_to_file(xml_path: str, outfile: str = "analysis.csv") -> str:
    """
    Preferred for large files when calling from JS.
    Reads XML from `xml_path` (already in Pyodide FS), writes CSV to `outfile`.
    Returns the outfile path.
    """
    txt_out = "input_parsed.txt"
    convert_factsage_xml_to_txt(xml_path, txt_out)
    hits = run_liquid_analysis(txt_out, top_k=3)

    import csv
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank",
            "Fe mass (g)",
            "Best T (°C)",
            "Si wt% (max)",
            "Liquid#1 composition",
            "Stopped at T",
            "First precipitates",
            "Slag first T",
            "Slag phases",
            "Slag liquid (before first solids)"
        ])
        for idx, r in enumerate(hits, 1):
            writer.writerow([
                idx,
                "not found" if r["Fe_g"] is None else fmt3(r["Fe_g"]),
                "not found" if r["best_T"] is None else f"{r['best_T']:.2f}",
                "not found" if r["best_Si"] is None else fmt3(r["best_Si"]),
                comp_text_with_total(r["best_comp"]) if r["best_comp"] else "not found",
                "not found" if r["stop_T"] is None else f"{r['stop_T']:.2f}",
                ", ".join(r["stop_phases"]) if r["stop_phases"] else "not found",
                "not found" if r["slag_first_T"] is None else f"{r['slag_first_T']:.2f}",
                ", ".join(r["slag_first_phases"]) if r["slag_first_phases"] else "not found",
                comp_text_with_total(r["slag_comp_before"]) if r["slag_comp_before"] else "not found"
            ])

    return outfile
