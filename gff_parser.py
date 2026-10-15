def parse_attributes(attr_text):
    """'ID=a;Name=b' -> {'ID': 'a', 'Name': 'b'}"""
    attributes = {}
    for pair in attr_text.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            attributes[key] = value
    return attributes


def parse_gff3(path):
    """Return {gene_id: {"cds_start": int, "cds_end": int, "product": str}}"""
    genes = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.rstrip("\n").split("\t")
            seqid, feature = fields[0], fields[2]
            start, end = int(fields[3]), int(fields[4])
            attributes = parse_attributes(fields[8])

            if seqid not in genes:
                genes[seqid] = {}
            if feature == "mRNA":
                genes[seqid]["product"] = attributes.get("product", "")
            elif feature == "CDS":
                genes[seqid]["cds_start"] = start
                genes[seqid]["cds_end"] = end
    return genes


def expected_protein_length(gene_info):
    """CDS length / 3 minus 1 stop codon."""
    cds_length = gene_info["cds_end"] - gene_info["cds_start"] + 1
    return cds_length // 3 - 1
