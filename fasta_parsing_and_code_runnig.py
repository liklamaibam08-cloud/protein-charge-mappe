from calculator import compute_charge_profile
from visualizer import plot_charge_distribution

def parsefasta(text):
    
    lines = text.split("\n")
    clean_sequence = ""

    for line in lines:
        line = line.strip()
        
        if line.startswith(">"):
            continue
            
        clean_sequence = clean_sequence + line

    return clean_sequence.upper()
