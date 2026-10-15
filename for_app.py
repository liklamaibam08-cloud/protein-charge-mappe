import streamlit as st
import matplotlib.pyplot as plt
from calculator import compute_charge_profile

st.title("Protein Charge Distribution Mapper")
st.write("Calculation and Visualization of a Protein Charge Sequnece:")
st.header("Step 1: Enter the protein sequence")
uploaded_file = st.file_uploader("Upload a FASTA file (.txt or .fasta)")
protein_sequence = ""

if uploaded_file is not None:
    file_bytes = uploaded_file.read()
    file_text = file_bytes.decode("utf-8")

lines = file_text.split("\n")
    clean_text = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            continue
        clean_text = clean_text + line
        
    protein_sequence = clean_text.upper()
    st.success("DONE")

else:
    protein_sequence = typed_text.upper()
    protein_sequence = protein_sequence.replace(" ", "")
    protein_sequence = protein_sequence.replace("\n", "")

preset_window = 11
st.header("Adjust the pH")
chosen_ph = st.slider("Select pH Level", 1.0, 14.0, 7.4)

st.header("Generating the Charge Map")

if st.button("Generate Charge Map"):
    if protein_sequence == "":
        st.warning("Please enter a protein sequence or upload a file first!")
    else:
        positions, profile = compute_charge_profile(protein_sequence, window_size=preset_window, ph=chosen_ph)
        
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(positions, profile, color="blue", linewidth=2)
        ax.axhline(0, color="gray", linestyle="--")
        
        ax.set_title(f"Protein Charge Distribution Profile (pH = {chosen_ph})")
        ax.set_xlabel("Amino Acid Position")
        ax.set_ylabel("Net Charge")
        
        st.pyplot(fig)
