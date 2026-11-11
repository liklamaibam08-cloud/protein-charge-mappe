PKA_VALUES = {
    'D': 3.9,
    'E': 4.3,
    'H': 6.0,
    'C': 8.3,
    'Y': 10.1,
    'K': 10.5,
    'R': 12.5
}

def compute_amino_acid_charge(aa, ph):
    aa = aa.upper()
    charge = 0.0
    if aa in ['K', 'R', 'H']:
        charge = 1.0 / (1.0 + 10.0**(ph - PKA_VALUES[aa]))
    elif aa in ['D', 'E', 'C', 'Y']:
        charge = -1.0 / (1.0 + 10.0**(PKA_VALUES[aa] - ph))
    return charge

def compute_charge_profile(sequence, window_size=11, ph=7.4):
    sequence = sequence.upper()
    positions = []
    profile = []
    half_window = window_size // 2
    
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        
        total_charge = sum(compute_amino_acid_charge(aa, ph) for aa in sequence[start:end])
        
        positions.append(i + 1)
        profile.append(total_charge)
        
    return positions, profile
