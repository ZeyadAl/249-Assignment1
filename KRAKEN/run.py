import os

texts = [
"../texts/GCF_000005845.2_ASM584v2_genomic.fna",
"../texts/GCF_000006765.1_ASM676v1_genomic.fna",
"../texts/GCF_000009045.1_ASM904v1_genomic.fna",
"../texts/GCF_000013425.1_ASM1342v1_genomic.fna",
"../texts/GCF_000195955.2_ASM19595v2_genomic.fna"
]


for i in range(len(texts)):
    #print("copying db ...")
    #os.system(f"cp -r OG mygendb{i}")
    os.system(f"../kraken2/kraken2-build --no-masking --add-to-library {texts[i]} --db mistask")
    #os.system(f"rm -rf mygendb{i}")


## ../kraken2/kraken2-build --build --db lab

os.system(f"../kraken2/kraken2-build --build --db mistask")
os.system(f"../kraken2/kraken2 --db mistask ../reads/reads_header_m.fastq --report final/Msample_report.txt --output final/Moutput.txt")
