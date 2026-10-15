PKA_VALUES = {
    'N_term': 8.0,
    'C_term': 3.1,
    'D': 3.9,
    'E': 4.1,
    'H': 6.0,
    'C': 8.5,
    'Y': 10.1,
    'K': 10.5,
    'R': 12.5
}
def calculate_residue_charge(residue, position, length, ph):
  res = residue.upper()
  charge = 0.0

    if position == 0:
        pka = PKA_VALUES['N_term']
        charge = charge + (1.0 / (1.0 + 10 ** (ph - pka)))

    if position == length - 1:
        pka = PKA_VALUES['C_term']
        charge = charge - (1.0 / (1.0 + 10 ** (pka - ph)))

    if res == 'K' or res == 'R':
        pka = PKA_VALUES[res]
        charge = charge + (1.0 / (1.0 + 10 ** (ph - pka)))
    elif res == 'H':
        pka = PKA_VALUES['H']
        charge = charge + (1.0 / (1.0 + 10 ** (ph - pka)))

    elif res == 'D' or res == 'E':
        pka = PKA_VALUES[res]
        charge = charge - (1.0 / (1.0 + 10 ** (pka - ph)))
    elif res == 'C' or res == 'Y':
        pka = PKA_VALUES[res]
        charge = charge - (1.0 / (1.0 + 10 ** (pka - ph)))
        
    return charge



def sliding_window_charge(sequence, window_size=11, ph=7.4):
    length = len(sequence)
    individual_charges = []
    
    for i in range(length):
        res = sequence[i]
        chg = calculate_residue_charge(res, i, length, ph)
        individual_charges.append(chg)
        
    half_window = window_size // 2
    profile = []
    
    for i in range(length):
        start = max(0, i - half_window)
        end = min(length, i + half_window + 1)
        window_slice = individual_charges[start:end]
        total_window_charge = sum(window_slice)
        profile.append(total_window_charge)
        
    positions = np.arange(1, length + 1)
    return positions, profile
