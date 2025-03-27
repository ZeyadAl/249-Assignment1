import os
import time
db = [
"makeblastdb -in ../texts/GCF_000005845.2_ASM584v2_genomic.fna -dbtype nucl -out Genome/genome1_db",
"makeblastdb -in ../texts/GCF_000009045.1_ASM904v1_genomic.fna -dbtype nucl -out Genome/genome2_db",
"makeblastdb -in ../texts/GCF_000006765.1_ASM676v1_genomic.fna -dbtype nucl -out Genome/genome3_db",
"makeblastdb -in ../texts/GCF_000013425.1_ASM1342v1_genomic.fna -dbtype nucl -out Genome/genome4_db",
"makeblastdb -in ../texts/GCF_000195955.2_ASM19595v2_genomic.fna -dbtype nucl -out Genome/genome5_db"
]


t0 = time.time()
for j in range(len(db)):
    print(f"Making the database for {db[j]}:")
    os.system(db[j])
t1 = time.time()
print(f"Time for creating the DBs is {t1-t0}")

for i in range(len(db)):
    print(f"Querying the database for {i}...")
    #r = f"blastn -query ../reads/reads_header.fasta -db Genome/genome{i+1}_db -out Results/results{i+1}.txt -outfmt 6 -task megablast -perc_identity 100 -word_size 100"
    # r = f"blastn -query ../reads/reads_header.fasta -db Genome/genome{i+1}_db -out Results/results{i+1}.csv -outfmt 6 -task blastn -perc_identity 99 -qcov_hsp_perc 100"
    r = f"blastn -query ../reads/reads_header_m.fasta -db Genome/genome{i+1}_db -out Results/Mresults{i+1}.csv -outfmt 6 -task blastn -perc_identity 99 -qcov_hsp_perc 100 -strand plus"
    #r = f"blastn -query ../reads/reads_header.fasta -db Genome/genome{i+1}_db -out Results/exact_results{i+1}.csv -outfmt 10 -task megablast -perc_identity 100 -qcov_hsp_perc 100"

    os.system(r)
    r = f"blastn -query ../reads/reads_header_m.fasta -db Genome/genome{i+1}_db -out Results/Eresults{i+1}.csv -outfmt 6 -task blastn -perc_identity 100 -qcov_hsp_perc 100 -strand plus"
    os.system(r)
t2 = time.time()
#print(f"Time for querying the DBs is {t2-t1}")

print(f"Total time for creating and querying is {t2-t0}")



