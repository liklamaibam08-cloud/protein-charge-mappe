def export_to_csv(positions, profile, filename="charge_profile.csv", sequence=None):
    """Save positions + charge profile to a .csv file.

    Columns: Position, [Residue], Net_Charge. Returns the filename written.
    """
    if len(positions) != len(profile):
        raise ValueError("positions and profile must have the same length")
    if sequence is not None and len(sequence) != len(positions):
        raise ValueError("sequence length must match positions length")

    with open(filename, "w") as f:
        if sequence is not None:
            f.write("Position,Residue,Net_Charge\n")
        else:
            f.write("Position,Net_Charge\n")
        for i in range(len(positions)):
            charge = round(profile[i], 6)
            if sequence is not None:
                f.write(f"{positions[i]},{sequence[i]},{charge}\n")
            else:
                f.write(f"{positions[i]},{charge}\n")
    return filename


def export_summary_csv(summary_rows, filename="summary.csv"):
    """summary_rows: list of dicts (same keys). Written with the keys as header."""
    columns = list(summary_rows[0].keys())
    with open(filename, "w") as f:
        f.write(",".join(columns) + "\n")
        for row in summary_rows:
            f.write(",".join(str(row[c]) for c in columns) + "\n")
    return filename
