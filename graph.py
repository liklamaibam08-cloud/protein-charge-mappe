import matplotlib.pyplot as plt

def plot_charge_distribution(positions, profile, protein_name="Protein Sequence", window_size=11, ph=7.4):
    
    plt.figure(figsize=(10, 4.5))
    
    plt.plot(positions, profile, color='blue', linewidth=2, label=f'Window Size: {window_size}')
    
    plt.axhline(0, color='gray', linestyle='--', alpha=0.7)
    
    plt.title(f'Protein Charge Distribution: {protein_name} (pH = {ph})', fontsize=12, fontweight='bold')
    plt.xlabel('Residue Position', fontsize=10)
    plt.ylabel('Net Charge (Window Sum)', fontsize=10)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    plt.show()
