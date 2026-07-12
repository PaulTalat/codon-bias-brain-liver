import ssl
import time
import csv
import json
import math
import os
from Bio import Entrez, SeqIO
from collections import Counter

# ============================================================
# 1. CONFIGURATION
# ============================================================
Entrez.email = "talatilyass07@gmail.com"
ssl._create_default_https_context = ssl._create_unverified_context

# ============================================================
# 2. REFERENCE GENE SET (for building CAI weights)
# ============================================================
# These are NOT brain/liver genes — this is a separate reference set used
# only to calculate what "preferred codon usage" looks like in highly-expressed genes.
REFERENCE_GENES = [
    ("RPL3", "NM_000967.4"),
    ("RPL4", "NM_000968.4"),
    ("RPL7", "NM_000971.4"),
    ("RPL8", "NM_001317782.2"),
    ("RPL11", "NM_000975.5"),
    ("RPL13", "NM_000977.4"),
    ("RPS2", "NM_002952.4"),
    ("RPS3", "NM_001005.5"),
    ("RPS6", "NM_001010.3"),
    ("RPS9", "NM_001013.4"),
]

# ============================================================
# 3. BIOLOGY REFERENCE — codon table
# ============================================================
# Note: 'M' (ATG) and 'W' (TGG) and stop codons are deliberately excluded.
# Each has no synonymous competition (single codon per amino acid)
# fixed single stop position), so RSCU/CAI are not meaningful for them.
aa_to_codons = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG']
}

# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================
def fetch_cds(accession):
    """Fetch the coding sequence (CDS) only — not full mRNA, not genomic DNA."""
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta_cds_na", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq).upper()

def sanity_check(seq):
    """The mandatory checkpoint: every valid CDS starts with ATG, has a length
    divisible by 3, and ends in a stop codon. If any of these fail, something
    was fetched wrong (e.g. full mRNA instead of CDS)."""
    if not seq.startswith('ATG'):
        return False, "Missing Start Codon"
    if len(seq) % 3 != 0:
        return False, "Length not divisible by 3"
    if seq[-3:] not in ['TAA', 'TAG', 'TGA']:
        return False, "Missing Stop Codon"
    return True, "Passed"

def build_reference_weights(reference_genes):
    """Pool codon counts across all reference (highly-expressed) genes, then
    compute relative adaptiveness w per codon = count / max count among its
    synonymous codons. This is the standard Sharp & Li (1987) CAI method,
    applied to a reference set you built and can fully explain."""
    pooled_counts = Counter()

    print("Building CAI reference weights from your ribosomal protein set...\n")
    for gene, accession in reference_genes:
        print(f"  Fetching reference gene {gene} ({accession})...")
        try:
            seq = fetch_cds(accession)
            ok, msg = sanity_check(seq)
            if not ok:
                print(f"    -> FAILED sanity check ({msg}), skipping this reference gene.")
                continue
            codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
            pooled_counts.update(codons)
        except Exception as e:
            print(f"    -> ERROR: {e}, skipping this reference gene.")
        time.sleep(0.4)  

    weights = {}
    for aa, syn_codons in aa_to_codons.items():
        counts = [pooled_counts[c] for c in syn_codons]
        max_count = max(counts) if counts else 0
        if max_count > 0:
            for codon in syn_codons:
                weights[codon] = round(pooled_counts[codon] / max_count, 3)

    with open('data/reference_weights.json', 'w') as f:
        json.dump(weights, f, indent=2)

    print(f"\nReference weights built from {len(reference_genes)} genes, saved to data/reference_weights.json\n")
    return weights

def calculate_metrics(sequence, weights):
    """Returns CAI, RSCU (per codon), GC content, and sequence length for one gene."""
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3) if len(sequence[i:i+3]) == 3]
    codon_counts = Counter(codons)

    # --- GC content ---
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = gc_count / len(sequence)

    # --- RSCU ---
    # Only calculated for amino acids with real synonymous competition
    # (aa_to_codons already excludes Met, Trp, and stop codons).
    rscu_scores = {}
    for aa, syn_codons in aa_to_codons.items():
        aa_total = sum(codon_counts[c] for c in syn_codons)
        if aa_total > 0:
            expected = aa_total / len(syn_codons)
            for codon in syn_codons:
                rscu_scores[codon] = round(codon_counts[codon] / expected, 3)

    # --- CAI (geometric mean of w-values) ---
    log_w_sum = 0
    valid_codon_count = 0
    for codon in codons[:-1]:  # skip the stop codon (last codon in sequence)
        if codon in ('ATG', 'TGG'):  # skip Met/Trp — no synonymous competition
            continue
        if codon in weights:
            log_w_sum += math.log(weights[codon])
            valid_codon_count += 1

    cai = math.exp(log_w_sum / valid_codon_count) if valid_codon_count > 0 else 0

    return round(cai, 3), rscu_scores, round(gc_content, 3), len(sequence)

# ============================================================
# 5. REFERENCE WEIGHTS
# ============================================================
if os.path.exists('data/reference_weights.json'):
    with open('data/reference_weights.json', 'r') as f:
        human_w = json.load(f)
    print("Loaded cached reference weights from data/reference_weights.json\n")
else:
    human_w = build_reference_weights(REFERENCE_GENES)

# ============================================================
# 6. MAIN PIPELINE — brain + liver gene list
# ============================================================
results = []
print("Starting Pipeline on brain/liver gene list...\n")

with open('data/gene_list.csv', 'r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        gene = row['Gene Name']
        group = row['Group']
        transcript = row['RefSeq Transcript ID']

        if not transcript or not transcript.startswith('NM_'):
            continue

        print(f"Fetching {gene} ({transcript}, {group})...")
        try:
            seq = fetch_cds(transcript)

            is_valid, msg = sanity_check(seq)
            if not is_valid:
                print(f"  -> FAILED: {msg}. Skipping.")
                continue

            cai, rscu, gc, length = calculate_metrics(seq, human_w)

            results.append({
                'gene': gene,
                'group': group,
                'cai': cai,
                'rscu_json': json.dumps(rscu),
                'seq_length': length,
                'gc_content': gc
            })
            print(f"  -> OK. CAI: {cai} | Length: {length}")

        except Exception as e:
            print(f"  -> ERROR: {e}")

        time.sleep(0.5)

# ============================================================
# 7. EXPORT
# ============================================================
with open('data/results.csv', 'w', newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=['gene', 'group', 'cai', 'rscu_json', 'seq_length', 'gc_content'])
    writer.writeheader()
    writer.writerows(results)

print(f"\nPipeline complete! Saved {len(results)} genes to data/results.csv")