import os
from calculator import compute_charge_profile
from fasta_parsing_and_code_runnig import parse_fasta_file
from gff_parser import parse_gff3, expected_protein_length
from exporter import export_to_csv, export_summary_csv
from graph import plot_charge_distribution

FASTA_PATH = "data/ref_proteins.faa"
GFF_PATH = "data/MR7_annotation.gff3"
OUT_DIR = "output"
PH = 7.4
WINDOW = 11


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    proteins = parse_fasta_file(FASTA_PATH)
    genes = parse_gff3(GFF_PATH)

    summary = []
    for gene_id, seq in proteins.items():
        positions, profile = compute_charge_profile(seq, window_size=WINDOW, ph=PH, include_termini=True)
        export_to_csv(positions, profile, f"{OUT_DIR}/{gene_id}_charge.csv", sequence=seq)
        plot_charge_distribution(positions, profile, protein_name=gene_id,
                                 window_size=WINDOW, ph=PH, save_path=f"{OUT_DIR}/{gene_id}_charge.png")
        expected = expected_protein_length(genes[gene_id]) if gene_id in genes else None
        summary.append({
            "gene": gene_id,
            "product": genes.get(gene_id, {}).get("product", ""),
            "length_aa": len(seq),
            "gff_expected_aa": expected,
            "length_matches_gff": len(seq) == expected,
            "max_window_charge": round(max(profile), 3),
            "min_window_charge": round(min(profile), 3),
        })

    export_summary_csv(summary, f"{OUT_DIR}/summary.csv")
    print(f"pH {PH}, window {WINDOW}")
    print(f"{'gene':<6}{'aa':>5}{'gff':>5}{'match':>7}{'max':>8}{'min':>8}  product")
    for r in summary:
        print(f"{r['gene']:<6}{r['length_aa']:>5}{r['gff_expected_aa']:>5}{str(r['length_matches_gff']):>7}"
              f"{r['max_window_charge']:>8}{r['min_window_charge']:>8}  {r['product']}")


if __name__ == "__main__":
    main()
