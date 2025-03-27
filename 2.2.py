from collections import defaultdict
import pandas as pd
'''
for g in texts:
    g = g.split(".")[0]
    T = sum(kmer_index[g].values())
    print(f"Total k-mers for {g}: {T:,d}")
print(f"Theoretically possible 31-mers: {4**31:,d}")
'''

def build_kmer_index(genome_files, k=31):
    kmer_index = defaultdict(lambda: defaultdict(int))

    for genome_file in genome_files:
        with open(genome_file, "r") as f:
            genome_seq = "".join(line.strip() for line in f if not line.startswith(">"))  
        
        for i in range(len(genome_seq) - k + 1):
            kmer = genome_seq[i:i+k]
            kmer_index[kmer][genome_file] += 1

    return kmer_index

# Example usage

texts = [
    "GCF_000005845.2_ASM584v2_genomic.fna",
    "GCF_000009045.1_ASM904v1_genomic.fna",
    "GCF_000006765.1_ASM676v1_genomic.fna",
    "GCF_000013425.1_ASM1342v1_genomic.fna",
    "GCF_000195955.2_ASM19595v2_genomic.fna"
]
names = {
"texts/GCF_000005845.2_ASM584v2_genomic.fna": "E. col",
"texts/GCF_000009045.1_ASM904v1_genomic.fna": "B. subtilis",
"texts/GCF_000006765.1_ASM676v1_genomic.fna": "P. aeruginosa",
"texts/GCF_000013425.1_ASM1342v1_genomic.fna":"S. aureus",
"texts/GCF_000195955.2_ASM19595v2_genomic.fna":"M. tuberculosis"
}

texts = ["texts/"+i for i in texts]
kmer_index = build_kmer_index(texts, k=31)
COUNTS = [
4641622,
4215576,
6264374,
2821331,
4411502
]

# The index has been created, maybe we should just save it somewhere

'''
print(f"Total unique k-mers: {len(kmer_index)}")
print(f"Theoretically possible 31-mers: {4**31}")

genome_length = sum(len(open(g).read().replace("\n", "")) for g in texts)
print(f"Total genome length across all genomes: {genome_length}")


'''
Total_unique_kmers = 22104695
Theoretically_possible_31mers = 4611686018427387904

def classify_reads(reads_file, kmer_index, k=31):
    read_matches = defaultdict(lambda: defaultdict(int))  # {read_id: {genome_id: count}}

    with open(reads_file, "r") as f:
        read_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                read_id = line[1:]  # Remove ">" from FASTA header
            else:
                # Extract k-mers from the read
                for i in range(len(line) - k + 1):
                    kmer = line[i:i+k]
                    if kmer in kmer_index:
                        for genome_id in kmer_index[kmer]:
                            read_matches[read_id][genome_id] += kmer_index[kmer][genome_id]

    return read_matches

reads_file = "reads/reads_header_m.fasta"
classified_reads = classify_reads(reads_file, kmer_index, k=31)

genome_match_counts = defaultdict(int)

for read_id, genome_counts in classified_reads.items():
    for genome, count in genome_counts.items():
        genome_match_counts[genome] += count

for i, (genome, count) in enumerate(genome_match_counts.items()):
    print(f"{names[genome]}: {count:,d} ({round(100*count/COUNTS[i], 2)}%)")
