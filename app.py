import streamlit as st
import matplotlib.pyplot as plt
from calculator import compute_charge_profile
from fasta_parsing_and_code_running import parse_fasta_text

st.title("Protein Charge Distribution Mapper")
st.write("Calculation and Visualization of a Protein Charge Sequence:")

st.header("Step 1: Enter the protein sequence")
uploaded_file = st.file_uploader("Upload a FASTA file (.txt or .fasta)")
typed_text = st.text_area("Or type/paste your sequence here:", value="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")

if uploaded_file is not None:
    fasta_text = uploaded_file.read().decode("utf-8")
else:
    fasta_text = typed_text

records = parse_fasta_text(fasta_text)

protein_sequence = ""
if len(records) > 0:
    chosen_protein = st.selectbox("Protein to analyse", list(records.keys()))
    protein_sequence = records[chosen_protein]

preset_window = 11
st.header("Adjust the pH")
chosen_ph = st.slider("Select pH Level", 1.0, 14.0, 7.4)

include_termini = st.checkbox("Include N- and C-terminal charges", value=False)

st.header("Step 2: Generate the Graph")

if st.button("Generate Charge Map"):
    if protein_sequence == "":
        st.warning("Please enter a protein sequence or upload a file first!")
    else:
        positions, profile = compute_charge_profile(protein_sequence, window_size=preset_window, ph=chosen_ph, include_termini=include_termini)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(positions, profile, color="blue", linewidth=2)
        ax.axhline(0, color="gray", linestyle="--")
        ax.set_title(f"Protein Charge Distribution Profile (pH = {chosen_ph})")
        ax.set_xlabel("Amino Acid Position")
        ax.set_ylabel("Net Charge")
        st.pyplot(fig)
