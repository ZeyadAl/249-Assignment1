from collections import defaultdict

def build_kmer_index(genome_files, k=31):
    kmer_index = defaultdict(lambda: defaultdict(int))

    for genome_file in genome_files:
        with open(genome_file, "r") as f:
            genome_seq = "".join(line.strip() for line in f if not line.startswith(">")) 
        
        # Extract k-mers
        for i in range(len(genome_seq) - k + 1):
            kmer = genome_seq[i:i+k]
            kmer_index[kmer][genome_file] += 1

    return kmer_index


texts = [
    "GCF_000005845.2_ASM584v2_genomic.fna",
    "GCF_000009045.1_ASM904v1_genomic.fna",
    "GCF_000006765.1_ASM676v1_genomic.fna",
    "GCF_000013425.1_ASM1342v1_genomic.fna",
    "GCF_000195955.2_ASM19595v2_genomic.fna"
]

texts = ["texts/"+i for i in texts]
kmer_index = build_kmer_index(texts, k=31)

print(f"Total unique k-mers: {len(kmer_index):,d}")
print(f"Theoretically possible 31-mers: {4**31:,d}")
