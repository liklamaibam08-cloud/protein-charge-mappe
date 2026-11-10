def parse_fasta_text(text):
    sequences = {}
    current_id = None
    for line in text.split("\n"):
        line = line.strip()
        if line == "":
            continue
        if line.startswith(">"):
            words = line[1:].split()
            current_id = words[0] if len(words) > 0 else "sequence"
            sequences[current_id] = ""
        else:
            if current_id is None:
                current_id = "sequence"
                sequences[current_id] = ""
            sequences[current_id] = sequences[current_id] + line.upper().replace(" ", "")
    return sequences
