"""
parse_ripl_gdr_slo.py
Parses RIPL-3 'gdr-parameters&errors-exp-SLO.dat' using 2-line paired fixed-width format.
Exports extracted parameters and uncertainties to a CSV file.
"""

import csv
from typing import Optional, List, Dict, Any


def parse_float(slice_str: str) -> Optional[float]:
    """Converts fixed-width string slice to float, returning None if empty."""
    s = slice_str.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def parse_slice(line: str, start: int, end: int) -> str:
    """Safely extracts a slice from a line that may be shorter than standard length."""
    if len(line) <= start:
        return ""
    return line[start : min(end, len(line))]


def parse_ripl_gdr_file(
    input_filepath: str, output_csvpath: str
) -> List[Dict[str, Any]]:
    # Step 1: Read and strip comment lines (starting with '#') and blank lines
    cleaned_lines: List[str] = []
    with open(input_filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.rstrip("\r\n")
            if (
                not stripped
                or stripped.lstrip().startswith("#")
                or stripped.lstrip().startswith("-")
            ):
                continue
            cleaned_lines.append(stripped)

    # Step 2: Ensure even number of lines (2 lines per record)
    if len(cleaned_lines) % 2 != 0:
        raise ValueError(
            f"Expected an even number of data lines (2 lines per record), got {len(cleaned_lines)} lines."
        )

    records: List[Dict[str, Any]] = []

    # Step 3: Iterate through pairs of (Line 1: Values, Line 2: Errors)
    for i in range(0, len(cleaned_lines), 2):
        l1 = cleaned_lines[i]
        l2 = cleaned_lines[i + 1]

        # Basic properties from Line 1
        z_str = parse_slice(l1, 2, 4).strip()
        a_str = parse_slice(l1, 5, 8).strip()
        el = parse_slice(l1, 9, 11).strip()
        reac_str = parse_slice(l1, 12, 14).strip()

        if not (z_str.isdigit() and a_str.isdigit()):
            continue

        z = int(z_str)
        a = int(a_str)
        reac = int(reac_str) if reac_str.isdigit() else None

        # First peak: (E1, CS1, W1) from Line 1; (Err) from Line 2
        e1 = parse_float(parse_slice(l1, 15, 20))
        e1_err = parse_float(parse_slice(l2, 15, 20))

        cs1 = parse_float(parse_slice(l1, 22, 27))
        cs1_err = parse_float(parse_slice(l2, 22, 27))

        w1 = parse_float(parse_slice(l1, 29, 34))
        w1_err = parse_float(parse_slice(l2, 29, 34))

        # Second peak: (E2, CS2, W2) from Line 1; (Err) from Line 2 (optional)
        e2 = parse_float(parse_slice(l1, 36, 41))
        e2_err = parse_float(parse_slice(l2, 36, 41))

        cs2 = parse_float(parse_slice(l1, 43, 48))
        cs2_err = parse_float(parse_slice(l2, 43, 48))

        w2 = parse_float(parse_slice(l1, 51, 56))
        w2_err = parse_float(parse_slice(l2, 51, 56))

        # Fitting range and reference code from Line 1
        e_min = parse_float(parse_slice(l1, 57, 61))
        e_max = parse_float(parse_slice(l1, 64, 68))
        ref = parse_slice(l1, 70, 75).strip()

        # Effective centroid energy: weighted average for 2-peak nuclei, E1 for 1-peak
        if cs1 is not None and e1 is not None:
            if e2 is not None and cs2 is not None and (cs1 + cs2) > 0:
                e_gdr_effective = (e1 * cs1 + e2 * cs2) / (cs1 + cs2)
            else:
                e_gdr_effective = e1
        else:
            e_gdr_effective = None

        record = {
            "Z": z,
            "A": a,
            "El": el,
            "Reac": reac,
            "E_GDR_eff": (
                f"{e_gdr_effective:.3f}" if e_gdr_effective is not None else ""
            ),
            "E1": f"{e1:.2f}" if e1 is not None else "",
            "E1_err": f"{e1_err:.2f}" if e1_err is not None else "",
            "CS1": f"{cs1:.1f}" if cs1 is not None else "",
            "CS1_err": f"{cs1_err:.1f}" if cs1_err is not None else "",
            "W1": f"{w1:.2f}" if w1 is not None else "",
            "W1_err": f"{w1_err:.2f}" if w1_err is not None else "",
            "E2": f"{e2:.2f}" if e2 is not None else "",
            "E2_err": f"{e2_err:.2f}" if e2_err is not None else "",
            "CS2": f"{cs2:.1f}" if cs2 is not None else "",
            "CS2_err": f"{cs2_err:.1f}" if cs2_err is not None else "",
            "W2": f"{w2:.2f}" if w2 is not None else "",
            "W2_err": f"{w2_err:.2f}" if w2_err is not None else "",
            "E_min": f"{e_min:.1f}" if e_min is not None else "",
            "E_max": f"{e_max:.1f}" if e_max is not None else "",
            "Ref": ref,
        }
        records.append(record)

    # Step 4: Write structured records to CSV
    fieldnames = [
        "Z",
        "A",
        "El",
        "Reac",
        "E_GDR_eff",
        "E1",
        "E1_err",
        "CS1",
        "CS1_err",
        "W1",
        "W1_err",
        "E2",
        "E2_err",
        "CS2",
        "CS2_err",
        "W2",
        "W2_err",
        "E_min",
        "E_max",
        "Ref",
    ]
    with open(output_csvpath, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"[✓] Extracted {len(records)} GDR entries into: {output_csvpath}")
    return records


if __name__ == "__main__":
    # Test on the provided snippet
    sample_filename = "gdr_snippet_sample.dat"
    with open(sample_filename, "w", encoding="utf-8") as f:
        f.write("""#   gdr-parameters&errors-exp-SLO.dat
# Z   A El Reac  E1     CS1     W1    E2     CS2    W2  E_{min}-E_{max} Ref
#-----------------------------------------------------------------------
   6  12 C   0 22.79   15.7   3.62                       13.0 - 25.0  63Bur
                0.07    0.7   0.28
   6  12 C   2 22.64   19.9   3.39                       20.1 - 25.0  03Var
                0.10    0.8   0.27
   6  13 C   2 24.40   13.3   7.70                       14.0 - 29.0  02Ish
                0.17    0.8   0.63
   6  14 C   2 16.68    7.5   3.10  25.87    9.4   6.84   9.0 - 30.0  02Ish
                0.57    1.6   0.67   0.17    0.5   0.84
""")

    parse_ripl_gdr_file(sample_filename, "gdr_parameters_slo.csv")


def main():
    ripl_file_path = "data/gdr-parameters&errors-exp-SLO.dat"
    df = parse_ripl_gdr_file(ripl_file_path, "data_processed/ripl3_gdr_exp.csv")
