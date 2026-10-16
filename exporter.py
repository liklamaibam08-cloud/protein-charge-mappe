def export_to_csv(positions, profile, sequence, filename="charge_profile.csv"):
    if len(positions) != len(profile) or len(sequence) != len(positions):
        raise ValueError("positions, profile and sequence must have the same length")
    with open(filename, "w") as f:
        f.write("Position,Residue,Net_Charge\n")
        for i in range(len(positions)):
            f.write(f"{positions[i]},{sequence[i]},{round(profile[i], 6)}\n")
    return filename
