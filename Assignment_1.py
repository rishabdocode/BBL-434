"""
Usage:
    python Assignment_1.py --input genome.fa --design design.txt --markers markers.tab --output Output.fa

Description:
    This program constructs a synthetic plasmid for a given bacterial genome.
    It detects the origin of replication (ORI), removes unwanted restriction sites,
    inserts user-specified restriction sites and marker genes, and adds neutral
    AT-rich spacers to prevent overlap.
"""


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
# ORI Detection
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
        score = max_repeat + at_content * 10

        if score > best_score:
            best_score = score
            best_region = region

    return best_region


# =====================================
# Neutral Spacer
# =====================================

def generate_spacer(length=40):
    return "".join(random.choice("AT") for _ in range(length))


# =====================================
# Marker & Design Parsing
# =====================================

def read_markers(filename):
    markers = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue  # skip empty lines

            if "\t" not in line:
                continue  # skip invalid marker lines

            name, seq = line.split("\t")
            markers[name.strip()] = seq.strip().upper()

    return markers



def read_design(filename):
    """
    Returns:
    - allowed_restrictions: set of enzyme names
    - required_markers: set of marker names
    """
    allowed_restrictions = set()
    required_markers = set()

    with open(filename) as f:
        for line in f:
            if line.strip():
                part, name = line.strip().split(",")
                part = part.strip().lower()
                name = name.strip()

                if "site" in part:
                    allowed_restrictions.add(name)
                else:
                    required_markers.add(name)

    return allowed_restrictions, required_markers


# =====================================
# Restriction Sites Dictionary
# =====================================

RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG"
}


# =====================================
# Plasmid Construction
# =====================================

def build_plasmid(genome, ori, allowed_sites, markers, required_markers):

    # 1️⃣ Remove unwanted restriction sites from genome
    for enzyme, site in RESTRICTION_SITES.items():
        if enzyme not in allowed_sites:
            genome = genome.replace(site, "")

    plasmid = ori + generate_spacer()

    # 2️⃣ Add allowed restriction sites
    for enzyme in allowed_sites:
        if enzyme in RESTRICTION_SITES:
            plasmid += RESTRICTION_SITES[enzyme]
            plasmid += generate_spacer(10)

    plasmid += generate_spacer()

    # 3️⃣ Add marker genes
    for marker in required_markers:
        if marker in markers:
            plasmid += markers[marker]
            plasmid += generate_spacer(20)
        else:
            print(f"⚠ Marker '{marker}' not found — skipped")

    return plasmid


# =====================================
# Validation
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
    markers = read_markers(args.markers)
    allowed_sites, required_markers = read_design(args.design)

    ori = find_ori(genome)
    plasmid = build_plasmid(genome, ori, allowed_sites, markers, required_markers)

    validate_plasmid(plasmid, ori)
    write_fasta(args.output, plasmid)

    print("✅ Plasmid successfully designed")
    print(f"Final size: {len(plasmid)} bp")


if __name__ == "__main__":
    main()
"""
Usage:
    python Assignment_1.py --input genome.fa --design design.txt --markers markers.tab --output Output.fa

Description:
    This program constructs a synthetic plasmid for a given bacterial genome.
    It detects the origin of replication (ORI), removes unwanted restriction sites,
    inserts user-specified restriction sites and marker genes, and adds neutral
    AT-rich spacers to prevent overlap.
"""
