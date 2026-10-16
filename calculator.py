PKA_VALUES = {
    'D': 3.9,
    'E': 4.3,
    'H': 6.0,
    'C': 8.3,
    'Y': 10.1,
    'K': 10.5,
    'R': 12.5
}

N_TERM_PKA = 9.0
C_TERM_PKA = 3.1

def compute_terminal_charges(ph):
    n_charge = 1.0 / (1.0 + 10.0**(ph - N_TERM_PKA))
    c_charge = -1.0 / (1.0 + 10.0**(C_TERM_PKA - ph))
    return n_charge, c_charge

def compute_amino_acid_charge(aa, ph):
    aa = aa.upper()
    charge = 0.0
    if aa in ['K', 'R', 'H']:
        charge = 1.0 / (1.0 + 10.0**(ph - PKA_VALUES[aa]))
    elif aa in ['D', 'E', 'C', 'Y']:
        charge = -1.0 / (1.0 + 10.0**(PKA_VALUES[aa] - ph))
    return charge

def compute_charge_profile(sequence, window_size=11, ph=7.4, include_termini=False):
    sequence = sequence.upper()
    n_charge, c_charge = compute_terminal_charges(ph) if include_termini else (0.0, 0.0)
    positions = []
    profile = []
    half_window = window_size // 2
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window_seq = sequence[start:end]
        total_charge = 0.0
        for aa in window_seq:
            total_charge += compute_amino_acid_charge(aa, ph)
        if start == 0:
            total_charge += n_charge
        if end == len(sequence):
            total_charge += c_charge
        positions.append(i + 1)
        profile.append(total_charge)
    return positions, profile
