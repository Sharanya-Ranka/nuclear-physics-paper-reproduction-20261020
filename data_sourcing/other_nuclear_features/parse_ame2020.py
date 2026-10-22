"""
parse_ame2020_mass.py
Parses AME2020 mass.mas20 / mass_1.mas20.txt into a tabular CSV.
Computes Sn, Sp, S2n, S2p, asymmetry, and total binding energy.
"""

import csv
import re
from typing import Optional

# Standard reference mass excesses from AME2020 (keV)
DELTA_NEUTRON_KEV = 8071.31806
DELTA_1H_KEV = 7288.971064


def clean_num(raw_str: str) -> Optional[float]:
    """Strips non-numeric flags (# estimates, *, whitespace) and returns float."""
    cleaned = re.sub(r"[#*~><=\s]", "", raw_str)
    if not cleaned:
        return None
    try:
        return float(cleaned)
    except ValueError:
        print(f"[!] Error converting to float: {cleaned}")
        return None


def parse_mass_table(input_path: str, output_csv_path: str) -> None:
    raw_records = []

    # Pass 1: Parse all raw rows based on AME2020 fixed-width Fortran specifications
    with open(input_path, "r", encoding="latin-1") as f:
        for line in f:
            # Skip page headers, comments, and empty lines
            if (
                len(line) < 90
                or line.startswith("1")
                or "MASS LIST" in line
                or "BINDING ENERGY" in line
            ):
                continue

            try:
                # Fortran 1-based format specification:
                # 2-4:   i3   NZ (N - Z)
                # 5-9:   i5   N
                # 10-14: i5   Z
                # 15-19: i5   A
                # 21-23: a3   EL (Element symbol)
                # 24-27: a4   O (Origin/flag)
                # 29-42: f14  MASS EXCESS (keV)
                # 43-54: f12  MASS EXCESS uncertainty (keV)
                # 55-67: f13  BINDING ENERGY/A (keV)
                # 69-78: f10  BINDING ENERGY uncertainty (keV)
                # 80-81: a2   BETA DECAY TYPE (B-)
                # 82-94: f13  BETA DECAY ENERGY (keV)
                # 95-105: f11 BETA DECAY uncertainty (keV)
                # 107-109: i3 ATOMIC MASS integer (u)
                # 111-123: f13 ATOMIC MASS fraction (micro-u)
                # 124-135: f12 ATOMIC MASS uncertainty (micro-u)

                n_str = line[4:9].strip()
                z_str = line[9:14].strip()
                a_str = line[14:19].strip()

                if not (n_str.isdigit() and z_str.isdigit() and a_str.isdigit()):
                    continue

                # breakpoint()
                n = int(n_str)
                z = int(z_str)
                a = int(a_str)
                el = line[20:23].strip()
                origin = line[23:27].strip()

                me = clean_num(line[28:42])
                me_unc = clean_num(line[42:54])
                bea = clean_num(line[54:67])
                bea_unc = clean_num(line[68:78])

                # breakpoint()

                beta_mode = line[79:81].strip()
                q_beta = clean_num(line[81:94])
                q_beta_unc = clean_num(line[94:105])

                # Parse Atomic Mass: integer + fraction in micro-u (10^-6 u)
                am_int_str = line[106:109].strip() if len(line) > 109 else ""
                am_frac_str = line[110:123].strip() if len(line) > 123 else ""
                am_frac = clean_num(am_frac_str)

                atomic_mass = None
                if am_int_str.isdigit() and am_frac is not None:
                    atomic_mass = float(am_int_str) + (am_frac * 1e-6)

                raw_records.append(
                    {
                        "Z": z,
                        "A": a,
                        "N": n,
                        "Element": el,
                        "Origin_Flag": origin,
                        "ME": me,
                        "ME_unc": me_unc,
                        "BEA": bea,
                        "BEA_unc": bea_unc,
                        "Beta_Mode": beta_mode,
                        "Q_Beta": q_beta,
                        "Q_Beta_unc": q_beta_unc,
                        "Atomic_Mass_u": atomic_mass,
                    }
                )
            except Exception:
                print(f"[!] Error parsing line: {line.strip()}")
                raise

    # Build fast lookup hashmap: (Z, A) -> Mass Excess (keV)
    me_map = {(r["Z"], r["A"]): r["ME"] for r in raw_records if r["ME"] is not None}

    # Pass 2: Calculate separation energies, total binding energy, and asymmetry
    enriched_records = []
    for r in raw_records:
        z = r["Z"]
        a = r["A"]
        n = r["N"]
        me = r["ME"]

        # Asymmetry (N - Z) / A
        asymmetry = (n - z) / a if a > 0 else 0.0

        # Total Binding Energy = (B/A) * A (in keV)
        total_be = (r["BEA"] * a) if (r["BEA"] is not None and a > 0) else None

        # Separation energies in keV (divide by 1000 for MeV)
        sn = None
        sp = None
        s2n = None
        s2p = None

        if me is not None:
            # Sn = ME(Z, A-1) + ME(n) - ME(Z, A)
            if (z, a - 1) in me_map:
                sn = me_map[(z, a - 1)] + DELTA_NEUTRON_KEV - me

            # Sp = ME(Z-1, A-1) + ME(1H) - ME(Z, A)
            if (z - 1, a - 1) in me_map:
                sp = me_map[(z - 1, a - 1)] + DELTA_1H_KEV - me

            # S2n = ME(Z, A-2) + 2*ME(n) - ME(Z, A)
            if (z, a - 2) in me_map:
                s2n = me_map[(z, a - 2)] + (2 * DELTA_NEUTRON_KEV) - me

            # S2p = ME(Z-2, A-2) + 2*ME(1H) - ME(Z, A)
            if (z - 2, a - 2) in me_map:
                s2p = me_map[(z - 2, a - 2)] + (2 * DELTA_1H_KEV) - me

        # Having a f"{xxx:.5f}" entered here may cause a loss of precision, but I think its okay for paper reproduction purposes
        enriched_records.append(
            {
                "Z": z,
                "A": a,
                "N": n,
                "Element": r["Element"],
                "Origin_Flag": r["Origin_Flag"],
                "ME": f"{me:.5f}" if me is not None else "",
                "ME_unc": f"{r['ME_unc']:.5f}" if r["ME_unc"] is not None else "",
                "BEA": f"{r['BEA']:.5f}" if r["BEA"] is not None else "",
                "BEA_unc": f"{r['BEA_unc']:.5f}" if r["BEA_unc"] is not None else "",
                "Total_BE": f"{total_be:.3f}" if total_be is not None else "",
                "Sn": f"{sn:.5f}" if sn is not None else "",
                "Sp": f"{sp:.5f}" if sp is not None else "",
                "S2n": f"{s2n:.5f}" if s2n is not None else "",
                "S2p": f"{s2p:.5f}" if s2p is not None else "",
                "Q_Beta": f"{r['Q_Beta']:.4f}" if r["Q_Beta"] is not None else "",
                "Atomic_Mass_u": (
                    f"{r['Atomic_Mass_u']:.8f}"
                    if r["Atomic_Mass_u"] is not None
                    else ""
                ),
                "Asymmetry": f"{asymmetry:.6f}",
            }
        )

    # Write output to CSV
    fieldnames = [
        "Z",
        "A",
        "N",
        "Element",
        "Origin_Flag",
        "ME",
        "ME_unc",
        "BEA",
        "BEA_unc",
        "Total_BE",
        "Sn",
        "Sp",
        "S2n",
        "S2p",
        "Q_Beta",
        "Atomic_Mass_u",
        "Asymmetry",
    ]

    with open(output_csv_path, "w", newline="", encoding="utf-8") as out_file:
        writer = csv.DictWriter(out_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(enriched_records)

    print(
        f"[✓] Successfully wrote {len(enriched_records)} entries to {output_csv_path}"
    )


if __name__ == "__main__":
    parse_mass_table("mass_1.mas20.txt", "ame2020_mass_features.csv")
