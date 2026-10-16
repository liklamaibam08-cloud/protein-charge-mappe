def parse_fasta_text(text):
    """Multi-record FASTA text -> {gene_id: sequence}. gene_id = first word of the header."""
    sequences = {}
    current_id = None
    for line in text.split("\n"):
        line = line.strip()
        if line == "":
            continue
        if line.startswith(">"):
            current_id = line[1:].split()[0]
            sequences[current_id] = ""
        else:
            if current_id is None:
                current_id = "sequence"
                sequences[current_id] = ""
            sequences[current_id] = sequences[current_id] + line.upper()
    return sequences


def parse_fasta_file(path):
    with open(path) as f:
        return parse_fasta_text(f.read())
