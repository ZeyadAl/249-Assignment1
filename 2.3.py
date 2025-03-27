from collections import defaultdict
names = {
"texts/GCF_000005845.2_ASM584v2_genomic.fna": "E. col",
"texts/GCF_000009045.1_ASM904v1_genomic.fna": "B. subtilis",
"texts/GCF_000006765.1_ASM676v1_genomic.fna": "P. aeruginosa",
"texts/GCF_000013425.1_ASM1342v1_genomic.fna":"S. aureus",
"texts/GCF_000195955.2_ASM19595v2_genomic.fna":"M. tuberculosis"
}

COUNTS = [
4641613,
4215567,
6264365,
2821322,
4411493
]

def build_minimizer_index(genome_files, k=31, w=10):
    minimizer_index = defaultdict(lambda: defaultdict(int))
    
    for genome_file in genome_files:
        with open(genome_file, "r") as f:
            genome_seq = "".join(line.strip() for line in f if not line.startswith(">"))
        
        window_length = w + k - 1
        for i in range(0, len(genome_seq) - window_length + 1):
            window_seq = genome_seq[i:i+window_length]
            kmer_candidates = [window_seq[j:j+k] for j in range(w)]
            minimizer = min(kmer_candidates)
            minimizer_index[minimizer][genome_file] += 1
            #minimizer_index[genome_file][minimizer] += 1
    
    return minimizer_index


def classify_reads_minimizer(reads_file, minimizer_index, k=31, w=10):
    read_matches = defaultdict(lambda: defaultdict(int))
    
    with open(reads_file, "r") as f:
        read_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                read_id = line[1:]  # Remove >
            else:
                window_length = w + k - 1
                for i in range(0, len(line) - window_length + 1):
                    window_seq = line[i:i+window_length]
                    kmer_candidates = [window_seq[j:j+k] for j in range(w)]
                    minimizer = min(kmer_candidates)
                    if minimizer in minimizer_index:
                        for genome_id in minimizer_index[minimizer]:
                            read_matches[read_id][genome_id] += 1
                            #read_matches[read_id][genome_id] += minimizer_index[minimizer][genome_id]
                            
    return read_matches

texts = [
    "texts/GCF_000005845.2_ASM584v2_genomic.fna",
    "texts/GCF_000009045.1_ASM904v1_genomic.fna",
    "texts/GCF_000006765.1_ASM676v1_genomic.fna",
    "texts/GCF_000013425.1_ASM1342v1_genomic.fna",
    "texts/GCF_000195955.2_ASM19595v2_genomic.fna"
]
minimizer_index = build_minimizer_index(texts, k=31, w=10)
'''
for g in texts:
    T = sum(minimizer_index[g].values())
    #print(f"Total minimizer for {g}: {T:,d}")
    print(T)
exit()
'''
classified_reads_min = classify_reads_minimizer("reads/reads_header.fasta", minimizer_index, k=31, w=10)

genome_match_counts = defaultdict(int)

for read_id, genome_counts in classified_reads_min.items():
    for genome, count in genome_counts.items():
        genome_match_counts[genome] += count


print("Matches per genome:")

for i, (genome, count) in enumerate(genome_match_counts.items()):
    #print(f"{names[genome]}: {count:,d}")
    print(f"{names[genome]}: {count:,d} ({round(100*count/COUNTS[i], 2)}%)")
