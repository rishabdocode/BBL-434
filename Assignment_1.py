import argparse
import random
from collections import defaultdict

# =====================================
# FASTA I/O
# =====================================

def read_fasta(filename):
    seq = ""
    with open(filename) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


def write_fasta(filename, seq):
    with open(filename, "w") as f:
        f.write(">Synthetic_Plasmid\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i+70] + "\n")


# =====================================
# ORI Detection (Frequent k-mers + AT bias)
# =====================================

def find_ori(genome, k=9, window=500):
    best_score = -1
    best_region = genome[:window]

    for i in range(len(genome) - window):
        region = genome[i:i+window]
        freq = defaultdict(int)

        for j in range(len(region) - k):
            freq[region[j:j+k]] += 1

        max_repeat = max(freq.values())
        at_content = (region.count("A") + region.count("T")) / len(region)

        score = max_repeat + at_content * 10  # AT bias

        if score > best_score:
            best_score = score
            best_region = region

    return best_region


# =====================================
# Spacer DNA (Neutral Linker)
# =====================================

def generate_spacer(length=40):
    return "".join(random.choice(["A", "T"]) for _ in range(length))


# =====================================
# Markers & Design Parsing
# =====================================

def read_markers(filename):
    markers = {}
    with open(filename) as f:
        for line in f:
            if line.strip():
                name, seq = line.strip().split("\t")
                markers[name] = seq.upper()
    return markers


def read_design(filename):
    design = []
    with open(filename) as f:
        for line in f:
            if line.strip():
                part, name = line.strip().split(",")
                design.append(name.strip())
    return design


# =====================================
# Restriction Sites
# =====================================

RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC"
}


# =====================================
# Plasmid Construction
# =====================================

def build_plasmid(genome, design, markers):
    plasmid = ""

    # 1. ORI
    ori = find_ori(genome)
    plasmid += ori

    # 2. Spacer
    plasmid += generate_spacer()

    # 3. Restriction Sites Block
    for item in design:
        if item in RESTRICTION_SITES:
            plasmid += RESTRICTION_SITES[item]
            plasmid += generate_spacer(10)

    plasmid += generate_spacer()

    # 4. Marker Genes
    for item in design:
        if item in markers:
            plasmid += markers[item]
            plasmid += generate_spacer(20)

    return plasmid


# =====================================
# Validation Checks
# =====================================

def validate_plasmid(plasmid, ori):
    assert ori in plasmid, "ORI missing"
    assert len(plasmid) > len(ori), "Plasmid too small"
    print("✔ Validation passed")


# =====================================
# Main
# =====================================

def main():
    parser = argparse.ArgumentParser(description="Generalized Plasmid Designer")
    parser.add_argument("--input", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    genome = read_fasta(args.input)
    design = read_design(args.design)
    markers = read_markers(args.markers)

    ori = find_ori(genome)
    plasmid = build_plasmid(genome, design, markers)

    validate_plasmid(plasmid, ori)
    write_fasta(args.output, plasmid)

    print("✅ Plasmid successfully designed")
    print(f"Final size: {len(plasmid)} bp")


if __name__ == "__main__":
    main()
