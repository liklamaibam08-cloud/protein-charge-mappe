# calculator.py - Protein Charge Calculation Logic

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
    if aa in ['K', 'R']:
        charge = 1.0 / (1.0 + 10.0**(ph - PKA_VALUES[aa]))
    elif aa == 'H':
        charge = 1.0 / (1.0 + 10.0**(ph - PKA_VALUES[aa]))
    elif aa in ['D', 'E']:
        charge = -1.0 / (1.0 + 10.0**(PKA_VALUES[aa] - ph))
    elif aa in ['C', 'Y']:
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
        window_seq = sequence[start:end]
        
        total_charge = 0.0
        for aa in window_seq:
            total_charge += compute_amino_acid_charge(aa, ph)
            
        positions.append(i + 1)
        profile.append(total_charge)
        
    return positions, profile
