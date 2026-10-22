import csv
import math
import re
from typing import Optional, Tuple

# Multipliers to convert half-life units into SI seconds
TIME_UNIT_TO_SECONDS = {
    "ys": 1e-24,
    "zs": 1e-21,
    "as": 1e-18,
    "fs": 1e-15,
    "ps": 1e-12,
    "ns": 1e-9,
    "us": 1e-6,
    "ms": 1e-3,
    "s": 1.0,
    "m": 60.0,
    "h": 3600.0,
    "d": 86400.0,
    "y": 365.2425 * 86400.0,
    "ky": 1e3 * 365.2425 * 86400.0,
    "my": 1e6 * 365.2425 * 86400.0,
    "gy": 1e9 * 365.2425 * 86400.0,
    "py": 1e12 * 365.2425 * 86400.0,
}

LN2 = math.log(2.0)


def clean_number(raw_val: str) -> Optional[float]:
    """Strips systematic '#' flags, whitespace, and returns float."""
    cleaned = re.sub(r"[#*~><=\s]", "", raw_val)
    if not cleaned:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def parse_spin_parity(jp_str: str) -> Tuple[Optional[float], Optional[int]]:
    """
    Extracts ground-state spin J (as float) and parity pi (+1 or -1).
    Example inputs: '1/2+*', '0+', '3/2-#', '(1/2+)', '2-'
    """
    # Remove metadata flags, parentheses, brackets, and whitespace
    cleaned = re.sub(r"[\(\)#*\[\]\s]", "", jp_str)
    if not cleaned:
        return None, None

    # Parity check
    parity = None
    if "+" in cleaned:
        parity = 1
    elif "-" in cleaned:
        parity = -1

    # Spin check
    spin_candidate = cleaned.replace("+", "").replace("-", "")
    # If isospin / multiple states exist, take the first valid state
    spin_candidate = spin_candidate.split(";")[0].split(",")[0].split("T=")[0]

    spin = None
    try:
        if "/" in spin_candidate:
            num, den = spin_candidate.split("/")[:2]
            spin = float(num) / float(den)
        elif spin_candidate != "":
            spin = float(spin_candidate)
    except ValueError:
        spin = None

    return spin, parity


def parse_decay_constant(t_val_str: str, unit_str: str) -> Optional[float]:
    """
    Computes decay constant lambda = ln(2) / T_1/2 in s^-1.
    Stable nuclei ('stbl') have lambda = 0.0.
    Particle unstable nuclei ('p-unst') have unmeasured / resonance half-lives.
    """
    raw_val = t_val_str.strip().lower()
    raw_unit = unit_str.strip().lower()

    if "stbl" in raw_val or raw_unit == "stbl":
        return 0.0

    if "p-unst" in raw_val or not raw_val:
        return None

    val = clean_number(raw_val)
    if val is None or val <= 0:
        return None

    factor = TIME_UNIT_TO_SECONDS.get(raw_unit)
    if factor is None:
        return None

    t_seconds = val * factor
    return LN2 / t_seconds


def parse_nubase(input_path: str, output_csv_path: str) -> None:
    extracted_rows = []

    with open(input_path, "r", encoding="latin-1") as f:
        for line in f:
            # Skip comment lines and lines shorter than required fixed columns
            if line.startswith("#") or len(line) < 100:
                continue

            try:
                # Fortran 1-based index mapping:
                # 1:3   -> line[0:3]   AAA
                # 5:8   -> line[4:8]   ZZZi (ZZZ = line[4:7], i = line[7])
                # 19:31 -> line[18:31] Mass Excess in keV
                # 70:78 -> line[69:78] Half-life T
                # 79:80 -> line[78:80] Half-life unit
                # 89:102-> line[88:102] Jpi

                a_str = line[0:3].strip()
                zzz_str = line[4:7].strip()
                isomer_flag = line[7:8].strip()

                # Filter: must be integer A and Z, and i == '0' (ground state only)
                if not (a_str.isdigit() and zzz_str.isdigit()):
                    continue
                if isomer_flag != "0":
                    continue

                a = int(a_str)
                z = int(zzz_str)
                n = a - z
                asymmetry = (n - z) / a if a > 0 else 0.0

                # Mass Excess (keV)
                me_raw = line[18:31]
                me = clean_number(me_raw)

                # Decay constant (s^-1)
                t_val = line[69:78]
                unit_t = line[78:80]
                decay_const = parse_decay_constant(t_val, unit_t)

                # Spin and Parity
                jpi_raw = line[88:102]
                spin, parity = parse_spin_parity(jpi_raw)

                extracted_rows.append(
                    {
                        "Z": z,
                        "A": a,
                        "N": n,
                        "ME": f"{me:.4f}" if me is not None else "",
                        "Spin": f"{spin:.1f}" if spin is not None else "",
                        "Parity": parity if parity is not None else "",
                        "Decay_Const": (
                            f"{decay_const:.6e}" if decay_const is not None else ""
                        ),
                        "Asymmetry": f"{asymmetry:.6f}",
                    }
                )

            except Exception:
                continue

    # Write target attributes to CSV
    fieldnames = ["Z", "A", "N", "ME", "Spin", "Parity", "Decay_Const", "Asymmetry"]
    with open(output_csv_path, "w", newline="", encoding="utf-8") as out_file:
        writer = csv.DictWriter(out_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(extracted_rows)

    print(
        f"Extraction complete. {len(extracted_rows)} ground-state nuclides written to {output_csv_path}"
    )


if __name__ == "__main__":
    parse_nubase(input_path="nubase_4.mas20.txt", output_csv_path="nubase_features.csv")
