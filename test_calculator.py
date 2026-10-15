import math

import pytest

from calculator import (
    compute_amino_acid_charge,
    compute_charge_profile,
    compute_terminal_charges,
)
from exporter import export_to_csv

@pytest.mark.parametrize("aa", ["K", "R", "H"])
def test_basic_residues_positive_at_low_ph(aa):
    assert compute_amino_acid_charge(aa, 1.0) == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("aa", ["D", "E", "C", "Y"])
def test_acidic_residues_negative_at_high_ph(aa):
    assert compute_amino_acid_charge(aa, 14.0) == pytest.approx(-1.0, abs=1e-2)


def test_charge_is_half_at_pka():
    assert compute_amino_acid_charge("D", 3.9) == pytest.approx(-0.5)
    assert compute_amino_acid_charge("K", 10.5) == pytest.approx(0.5)


@pytest.mark.parametrize("aa", ["A", "G", "L", "X", "*"])
def test_neutral_or_unknown_residue_is_zero(aa):
    assert compute_amino_acid_charge(aa, 7.4) == 0.0


def test_lowercase_residue_handled():
    assert compute_amino_acid_charge("k", 7.4) == compute_amino_acid_charge("K", 7.4)

def test_net_charge_decreases_as_ph_rises():
    seq = "KRHDEKRHDE"
    totals = [sum(compute_charge_profile(seq, window_size=len(seq), ph=ph)[1][:1])
              for ph in (2, 5, 7.4, 10, 13)]
    assert totals == sorted(totals, reverse=True)


def test_extreme_ph_values_do_not_crash():
    for ph in (0.0, 14.0):
        _, profile = compute_charge_profile("MKRDEH", ph=ph, include_termini=True)
        assert all(math.isfinite(v) for v in profile)

def test_terminal_charges_at_physiological_ph():
    n, c = compute_terminal_charges(7.4)
    assert n > 0.9 and c < -0.99


def test_terminal_charges_at_extreme_ph():
    n_low, c_low = compute_terminal_charges(1.0)
    assert n_low == pytest.approx(1.0, abs=1e-3)
    assert c_low == pytest.approx(0.0, abs=0.01)
    n_high, c_high = compute_terminal_charges(14.0)
    assert n_high == pytest.approx(0.0, abs=1e-3)
    assert c_high == pytest.approx(-1.0, abs=1e-3)


def test_termini_off_by_default():
    _, profile = compute_charge_profile("AAAA", window_size=3, ph=7.4)
    assert profile == [0.0, 0.0, 0.0, 0.0]


def test_n_terminus_only_affects_first_window():
    n, c = compute_terminal_charges(7.4)
    _, profile = compute_charge_profile("A" * 30, window_size=5, ph=7.4, include_termini=True)
    assert profile[0] == pytest.approx(n)
    assert profile[-1] == pytest.approx(c)
    assert all(v == 0.0 for v in profile[3:-3])


def test_single_residue_gets_both_termini():
    n, c = compute_terminal_charges(7.4)
    _, profile = compute_charge_profile("A", window_size=11, ph=7.4, include_termini=True)
    assert profile == [pytest.approx(n + c)]

def test_empty_sequence():
    assert compute_charge_profile("", ph=7.4) == ([], [])


def test_profile_length_and_positions():
    positions, profile = compute_charge_profile("MKRDE", window_size=3)
    assert positions == [1, 2, 3, 4, 5]
    assert len(profile) == 5


def test_window_larger_than_sequence():
    _, profile = compute_charge_profile("KD", window_size=51, ph=7.4)
    assert profile[0] == pytest.approx(profile[1])


def test_lowercase_sequence_matches_uppercase():
    assert compute_charge_profile("mkrde") == compute_charge_profile("MKRDE")

def test_export_csv_roundtrip(tmp_path):
    positions, profile = compute_charge_profile("MKRDE", window_size=3)
    out = tmp_path / "out.csv"
    export_to_csv(positions, profile, str(out), sequence="MKRDE")
    with open(out) as f:
        rows = [line.strip().split(",") for line in f]
    assert rows[0] == ["Position", "Residue", "Net_Charge"]
    assert len(rows) == 6
    assert rows[1][:2] == ["1", "M"]
    assert float(rows[2][2]) == pytest.approx(profile[1], abs=1e-5)


def test_export_csv_length_mismatch(tmp_path):
    with pytest.raises(ValueError):
        export_to_csv([1, 2], [0.0], str(tmp_path / "bad.csv"))

from fasta_parsing_and_code_runnig import parse_fasta_file
from gff_parser import parse_gff3, expected_protein_length


def test_parse_fasta_file_multi_record(tmp_path):
    p = tmp_path / "t.faa"
    p.write_text(">a first\nmk\nrd\n>b second\nEE\n")
    assert parse_fasta_file(str(p)) == {"a": "MKRD", "b": "EE"}


def test_fasta_and_gff_agree_on_real_data():
    proteins = parse_fasta_file("data/ref_proteins.faa")
    genes = parse_gff3("data/MR7_annotation.gff3")
    assert set(proteins) == set(genes)
    for gene_id, seq in proteins.items():
        assert len(seq) == expected_protein_length(genes[gene_id])
