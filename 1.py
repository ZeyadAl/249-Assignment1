import numpy as np
import gzip
from AhoCorasick import AhoCorasick
import pprint 
pp = pprint.PrettyPrinter()

# we want to implement the Aho_Corasick algo

texts = [
    "GCF_000005845.2_ASM584v2_genomic.fna",
    "GCF_000009045.1_ASM904v1_genomic.fna",
    "GCF_000006765.1_ASM676v1_genomic.fna",
    "GCF_000013425.1_ASM1342v1_genomic.fna",
    "GCF_000195955.2_ASM19595v2_genomic.fna"
]

reads = [
    "simulated_reads_miseq_10k_R1.fastq.gz",
    "simulated_reads_miseq_10k_R2.fastq.gz",
    "simulated_reads_no_errors_10k_R1.fastq.gz",
    "simulated_reads_no_errors_10k_R2.fastq.gz"
]


string_size = 125
per_read = 5000
n = 2 * per_read

header_array = np.full(n, '', dtype=f'U{string_size}')
patterns = np.full(n, '', dtype=f'U{string_size}')

## extract reads
def extract_sequences(input_file, read_num, per_read=per_read):
    count = 0
    with gzip.open(input_file, 'rt', encoding='utf-8') as infile:
        while count < per_read:
            header_array[read_num*per_read + count] = infile.readline().strip()
            patterns[read_num*per_read + count] = infile.readline().strip()
            infile.readline()  # Skip '+'
            infile.readline()  # Skip quality scores
            count += 1

extract_sequences("reads_gz/"+reads[0], 0)
extract_sequences("reads_gz/"+reads[1], 1)



aho = AhoCorasick()

for i, pattern in enumerate(patterns):
    aho.add_pattern(pattern, i)

aho.build_automaton()

def load_genome(genome_file):
    with open(genome_file, 'r') as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))  # Ignore headers
    
# Read_Org_dic = {}
org_match = {}

for gcf in texts:

    genome_file = "texts/" + gcf
    genome_sequence = load_genome(genome_file)

    matches = aho.search(genome_sequence, patterns)
    print(f"{gcf}: number of matches:{len(matches)}")
    org_match[gcf] = set()
    for pattern_idx, pattern, start_pos in matches:  # Show first 10 matches
        org_match[gcf].add(header_array[pattern_idx])
        # if pattern_idx == 0:
        #print(f"{header_array[pattern_idx][1:]},{start_pos + 1}")
        
        # if header_array[pattern_idx] in Read_Org_dic:
        #     if gcf not in Read_Org_dic[header_array[pattern_idx]]:
        #         Read_Org_dic[header_array[pattern_idx]].append(gcf)
        # else:
        #     Read_Org_dic[header_array[pattern_idx]] = [gcf]

        # print(f"header: {header_array[pattern_idx]}")
        # print(f"from arr: {patterns[pattern_idx]}")
        # print(f"Pattern {pattern_idx} ({pattern}) found at position {start_pos}")
    # print(org_match)
    # pp.pprint(org_match)
    #break

# Print some results
# for pattern_idx, pattern, start_pos in matches[:10]:  # Show first 10 matches
#     print(f"Pattern {pattern_idx} ({pattern}) found at position {start_pos}")

