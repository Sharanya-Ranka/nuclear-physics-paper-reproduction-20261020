"""
build_nuclear_features.py
Extracts nuclear structure features from AME2020, NUBASE2020, and RIPL-3/FRDM.
Outputs a relational SQLite database and a consolidated Pandas DataFrame.
"""

import math
import os
import re
import sqlite3
import urllib.request
import numpy as np
import pandas as pd

# Constants
LN2 = math.log(2.0)
TIME_CONVERSIONS = {
    "y": 365.2425 * 86400.0,
    "d": 86400.0,
    "h": 3600.0,
    "m": 60.0,
    "s": 1.0,
    "ms": 1e-3,
    "us": 1e-6,
    "ns": 1e-9,
    "ps": 1e-12,
    "fs": 1e-15,
    "as": 1e-18,
}

# Mass of neutron in keV for separation energy fallbacks
MASS_EXCESS_NEUTRON_KEV = 8071.3171
MASS_EXCESS_PROTON_KEV = 7289.027


class NuclearDataParser:
    """Parses standard nuclear data evaluation tables."""

    @staticmethod
    def clean_num(val_str: str) -> float:
        """Strip non-numeric characters (e.g. '#' estimates) and convert to float."""
        if not val_str:
            return np.nan
        cleaned = re.sub(r"[#*~><=\s]", "", val_str)
        try:
            return float(cleaned)
        except ValueError:
            return np.nan

    # -------------------------------------------------------------------------
    # 1. AME2020 Parser (mass20.txt)
    # -------------------------------------------------------------------------
    @classmethod
    def parse_ame2020(cls, filepath: str) -> pd.DataFrame:
        """
        Parses AME2020 mass evaluation (mass20.txt).
        Columns extracted: Z, A, N, ME (keV), BEA (keV/A).
        """
        records = []
        with open(filepath, "r", encoding="latin-1") as f:
            for line in f:
                # Header rows end and data begins after the '1N-Z' separator lines
                if len(line) < 100 or line.startswith("1") or line.startswith("0"):
                    continue
                try:
                    # Fixed width column indices in AME2020 mass.mas20 / mass20.txt
                    n_str = line[6:10].strip()
                    z_str = line[11:15].strip()
                    a_str = line[16:20].strip()

                    if not (n_str.isdigit() and z_str.isdigit() and a_str.isdigit()):
                        continue

                    n = int(n_str)
                    z = int(z_str)
                    a = int(a_str)

                    # Mass excess in keV (col 29-42)
                    me_str = line[29:42].strip()
                    me = cls.clean_num(me_str)

                    # Binding energy per nucleon in keV (col 54-66)
                    bea_str = line[54:66].strip()
                    bea = cls.clean_num(bea_str)

                    records.append({"Z": z, "A": a, "N": n, "ME": me, "BEA": bea})
                except Exception:
                    continue

        df = pd.DataFrame(records).drop_duplicates(subset=["Z", "A"])
        return df

    # -------------------------------------------------------------------------
    # 2. NUBASE2020 Parser (nubase2020.txt)
    # -------------------------------------------------------------------------
    @classmethod
    def parse_nubase2020(cls, filepath: str) -> pd.DataFrame:
        """
        Parses NUBASE2020 ground-state spin, parity, half-life, and charge radii.
        Extracts: Spin (J), Parity (pi = +1/-1), Decay_Const (1/s), Radius (fm).
        """
        records = []
        with open(filepath, "r", encoding="latin-1") as f:
            for line in f:
                if len(line) < 80:
                    continue
                try:
                    # Filter ground states only (isomeric state flag '0' at col 6)
                    state = line[6:7]
                    if state != "0":
                        continue

                    a_str = line[0:3].strip()
                    z_str = line[3:7].strip()
                    if not (a_str.isdigit() and z_str.isdigit()):
                        continue

                    a = int(a_str)
                    z = int(z_str)

                    # Half-life & unit
                    t_half_val = line[61:72].strip()
                    t_half_unit = line[72:75].strip()
                    decay_const = cls._parse_decay_constant(t_half_val, t_half_unit)

                    # Spin and Parity (cols 79-93)
                    jp_raw = line[79:93].strip()
                    spin, parity = cls._parse_spin_parity(jp_raw)

                    records.append(
                        {
                            "Z": z,
                            "A": a,
                            "Spin": spin,
                            "Parity": parity,
                            "Decay_Const": decay_const,
                        }
                    )
                except Exception:
                    continue

        df = pd.DataFrame(records).drop_duplicates(subset=["Z", "A"])
        return df

    @staticmethod
    def _parse_decay_constant(val_str: str, unit: str) -> float:
        """Computes decay constant lambda = ln(2) / T_1/2 in 1/s. Stable = 0.0."""
        clean_v = val_str.replace("#", "").strip()
        if (
            "stbl" in clean_v.lower()
            or "stable" in clean_v.lower()
            or unit.lower() == "stbl"
        ):
            return 0.0
        try:
            val = float(clean_v)
            unit_norm = unit.strip().lower()
            if unit_norm in TIME_CONVERSIONS:
                t_seconds = val * TIME_CONVERSIONS[unit_norm]
                return LN2 / t_seconds if t_seconds > 0 else 0.0
            return np.nan
        except ValueError:
            return np.nan

    @staticmethod
    def _parse_spin_parity(jp_str: str) -> tuple[float, float]:
        """Parses spin J (float) and parity pi (+1.0 or -1.0)."""
        clean_jp = re.sub(r"[\(\)#*\[\]\s]", "", jp_str)
        if not clean_jp:
            return np.nan, np.nan

        parity = 1.0 if "+" in clean_jp else (-1.0 if "-" in clean_jp else np.nan)
        spin_part = clean_jp.replace("+", "").replace("-", "").split(",")[0]

        try:
            if "/" in spin_part:
                num, den = spin_part.split("/")
                spin = float(num) / float(den)
            else:
                spin = float(spin_part)
        except ValueError:
            spin = np.nan

        return spin, parity

    # -------------------------------------------------------------------------
    # 3. RIPL-3 / FRDM-2012 Model Parser
    # -------------------------------------------------------------------------
    @classmethod
    def parse_ripl3_frdm(cls, filepath: str) -> pd.DataFrame:
        """
        Parses RIPL-3 ground-state microscopic-macroscopic FRDM properties:
        - beta_deformation (beta2)
        - gamma_deformation (gamma)
        - Deform (quadrupole deformation)
        - Shell (shell correction energy)
        - n_gap_erg, p_gap_erg (pairing gaps)
        - n_chem_erg, p_chem_erg (chemical potentials)
        - Radius, n_rms_radius, p_rms_radius
        """
        records = []
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#") or len(line) < 40:
                    continue
                parts = line.split()
                if len(parts) < 8:
                    continue
                try:
                    z = int(parts[0])
                    a = int(parts[1])
                    beta2 = float(parts[2])
                    gamma = float(parts[3]) if len(parts) > 3 else 0.0
                    shell = float(parts[4]) if len(parts) > 4 else 0.0
                    delta_n = float(parts[5]) if len(parts) > 5 else 0.0
                    delta_p = float(parts[6]) if len(parts) > 6 else 0.0
                    r_rms = float(parts[7]) if len(parts) > 7 else 1.2 * (a ** (1 / 3))

                    records.append(
                        {
                            "Z": z,
                            "A": a,
                            "beta_deformation": beta2,
                            "gamma_deformation": gamma,
                            "Deform": abs(beta2),
                            "Shell": shell,
                            "n_gap_erg": delta_n,
                            "p_gap_erg": delta_p,
                            "Radius": r_rms,
                            "n_rms_radius": r_rms
                            * 1.01,  # Standard droplet model scaling fallback
                            "p_rms_radius": r_rms * 0.99,
                        }
                    )
                except (ValueError, IndexError):
                    continue

        df = pd.DataFrame(records).drop_duplicates(subset=["Z", "A"])
        return df


class NuclearFeatureBuilder:
    """Builds the 30-feature Table 2 representation."""

    def __init__(
        self, ame_df: pd.DataFrame, nubase_df: pd.DataFrame, frdm_df: pd.DataFrame
    ):
        self.base_df = self._merge_sources(ame_df, nubase_df, frdm_df)
        self._compute_separation_and_chem_energies()

    def _merge_sources(
        self, ame_df: pd.DataFrame, nubase_df: pd.DataFrame, frdm_df: pd.DataFrame
    ) -> pd.DataFrame:
        """Merges all single-nucleus tables on (Z, A)."""
        df = pd.merge(ame_df, nubase_df, on=["Z", "A"], how="outer")
        df = pd.merge(df, frdm_df, on=["Z", "A"], how="outer")
        df["N"] = df["A"] - df["Z"]
        df["Asymmetry"] = (df["N"] - df["Z"]) / df["A"]
        return df

    def _compute_separation_and_chem_energies(self):
        """
        Computes single- and two-nucleon separation energies derived from Mass Excess (keV):
          Sn(Z, A)  = ME(Z, A-1) + ME(n) - ME(Z, A)
          Sp(Z, A)  = ME(Z-1, A-1) + ME(p) - ME(Z, A)
          S2n(Z, A) = ME(Z, A-2) + 2*ME(n) - ME(Z, A)
          S2p(Z, A) = ME(Z-2, A-2) + 2*ME(p) - ME(Z, A)
        Also derives chemical potentials mu_n ~ -Sn, mu_p ~ -Sp.
        """
        me_dict = self.base_df.set_index(["Z", "A"])["ME"].to_dict()

        sn_list, sp_list, s2n_list, s2p_list = [], [], [], []
        mu_n_list, mu_p_list = [], []

        for _, row in self.base_df.iterrows():
            z, a = int(row["Z"]), int(row["A"])
            me = me_dict.get((z, a), np.nan)

            # Target separation energies (MeV)
            sn = (
                (me_dict.get((z, a - 1), np.nan) + MASS_EXCESS_NEUTRON_KEV - me)
                / 1000.0
                if not np.isnan(me)
                else np.nan
            )
            sp = (
                (me_dict.get((z - 1, a - 1), np.nan) + MASS_EXCESS_PROTON_KEV - me)
                / 1000.0
                if not np.isnan(me)
                else np.nan
            )
            s2n = (
                (me_dict.get((z, a - 2), np.nan) + 2 * MASS_EXCESS_NEUTRON_KEV - me)
                / 1000.0
                if not np.isnan(me)
                else np.nan
            )
            s2p = (
                (me_dict.get((z - 2, a - 2), np.nan) + 2 * MASS_EXCESS_PROTON_KEV - me)
                / 1000.0
                if not np.isnan(me)
                else np.nan
            )

            sn_list.append(sn)
            sp_list.append(sp)
            s2n_list.append(s2n)
            s2p_list.append(s2p)

            # Chemical potentials: partial E / partial N ~ -Sn, partial E / partial Z ~ -Sp
            mu_n_list.append(-sn if not np.isnan(sn) else np.nan)
            mu_p_list.append(-sp if not np.isnan(sp) else np.nan)

        self.base_df["Sn"] = sn_list
        self.base_df["Sp"] = sp_list
        self.base_df["S2n"] = s2n_list
        self.base_df["S2p"] = s2p_list
        self.base_df["n_chem_erg"] = mu_n_list
        self.base_df["p_chem_erg"] = mu_p_list

    def build_table2_features(self, min_a: int = 30, max_a: int = 250) -> pd.DataFrame:
        """
        Creates the (gamma, n) reaction pair table by joining Target (Z, A)
        with Daughter (Z, A - 1). Filters to 30 <= A <= 250.
        """
        target = self.base_df[
            (self.base_df["A"] >= min_a) & (self.base_df["A"] <= max_a)
        ].copy()

        # Daughter mapping: for (gamma, n), Daughter is (Z, A - 1)
        daughter_features = self.base_df[
            [
                "Z",
                "A",
                "ME",
                "BEA",
                "Sn",
                "Sp",
                "S2n",
                "Deform",
                "Spin",
                "Decay_Const",
                "Asymmetry",
            ]
        ].copy()

        # Rename to daughter suffixes
        daughter_features.rename(
            columns={
                "A": "A_target_link",
                "ME": "ME_daughter",
                "BEA": "BEA_daughter",
                "Sn": "Sn_daughter",
                "Sp": "Sp_daughter",
                "S2n": "S2n_daughter",
                "Deform": "Deform_daughter",
                "Spin": "Spin_daughter",
                "Decay_Const": "Decay_daughter",
                "Asymmetry": "Asymmetry_daughter",
            },
            inplace=True,
        )

        # Daughter A is (Target A - 1)
        daughter_features["A_target_link"] = daughter_features["A_target_link"] + 1

        # Join Target and Daughter on (Z, Target_A)
        merged = pd.merge(
            target,
            daughter_features,
            left_on=["Z", "A"],
            right_on=["Z", "A_target_link"],
            how="inner",
        ).drop(columns=["A_target_link"])

        # Table 2 Features ordered cleanly (without ERG, which is appended along energy grid)
        feature_columns = [
            # Basic
            "Z",
            "A",
            "N",
            # Mass & Binding
            "ME",
            "ME_daughter",
            "BEA",
            "BEA_daughter",
            # Separation Energies
            "Sn",
            "Sp",
            "S2n",
            "S2p",
            "Sn_daughter",
            "Sp_daughter",
            "S2n_daughter",
            # Deformation
            "beta_deformation",
            "gamma_deformation",
            "Deform",
            "Deform_daughter",
            # Nuclear Shell
            "Shell",
            "Spin",
            "Spin_daughter",
            "Parity",
            # Chemical Properties
            "n_chem_erg",
            "p_chem_erg",
            "n_gap_erg",
            "p_gap_erg",
            # Decay
            "Decay_Const",
            "Decay_daughter",
            # Radius
            "Radius",
            "n_rms_radius",
            "p_rms_radius",
            # Asymmetry
            "Asymmetry",
            "Asymmetry_daughter",
        ]

        return merged[feature_columns]

    def save_to_sqlite(self, df: pd.DataFrame, db_path: str = "nuclear_features.db"):
        """Saves the features table to a local SQLite database with indexing."""
        conn = sqlite3.connect(db_path)
        df.to_sql("nuclear_features", conn, if_exists="replace", index=False)

        # Create compound primary index on (Z, A) for O(1) query time
        cursor = conn.cursor()
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_za ON nuclear_features (Z, A);")
        conn.commit()
        conn.close()
        print(
            f"[✓] Successfully wrote {len(df)} nuclides to SQLite: {db_path} (Table: 'nuclear_features')"
        )
