# Report for Assignment 1
### By Zeyad Aljaali


## Task 1

### Task 1.1: Multiple Matches
We used the Aho-Corasick algorithm to search through all the genomes at once. We combined the reads into one file, for easier manipulation (of course the one with errors in a seperate file than the one without errors). Then we checked each read separately to see which organisms it matches. Aho-Corasick is an efficient algorithm that finds all matching sequences in one pass. This keeps the overall process efficient even when working with large collections of reads. 


### Task 1.2: Exact Matching

We implemented the Aho-Corasick algorithm in file called 1.py where we just need the genomes are in texts file, and reads in the reads file.

Results:
For no errors:
```
                matches  (ratio)
E. coli:        3289     (32.89%)

B. subtilis:     549     (5.49%)

P. aeruginosa:   512     (5.12%)

S. aureus:       520     (5.2%)

M. tuberculosis: 571     (5.71%)
```
For miseq:
```
E. coli:        2750    (27.5%)

B. subtilis:     454    (4.54%)

P. aeruginosa:   432    (4.32%)

S. aureus:       464    (4.64%)

M. tuberculosis: 445    (4.45%)
```
### Task 1.3: Approximate Matching
Not implemented 


### Task 1.4: Comparison with Bioinformatics tools
These times are for no errors:

Aho-Corasick time: 3.95s, peak memory: 466MB

BLAST: time: 27.56s , peak memory: 203MB

These times are for miseq:

Aho-Corasick time: 4.07s, peak memory: 462MB

BLAST: time: 55.7s, peak memory: 142MB

BLAST results:
For no errors:
```
                  matches  (ratio)
E. col:             3321   (33.21%)

B. subtilis:         552   (5.52%)

P. aeruginosa:       515   (5.15%)

S. aureus:           523   (5.23%)

M. tuberculosis:     574   (5.74%)
```
For miseq:
```
E. col:             1878   (18.78%)

B. subtilis:         304   (3.04%)

P. aeruginosa:       283   (2.83%)

S. aureus:           334   (3.34%)

M. tuberculosis:     278   (2.78%)
```
For the no errors reads, the differences in results is due to having up-to one mismatch from BLAST. I looked into exact matching and we get the same exact results.
As for the miseq reads, there's a stark difference. This is because BLAST depends on MSP scores and local alignment, so it's not really checking for exact matches, so this will lead to discrepancies. Also another reason is being that we are checking for up-to-one mismatches, we'd think this will give us a higher number of matches.

## Task 2

### Task 2.1: Multiple Matches

• Briefly describe your data structure and justify your choices
We used a nested defaultdict structure. The Outer Dictionary
Maps each unique 31-mer (a string of nucleotides) to an inner dictionary. This provides fast lookup for any k-mer. The Inner Dictionary: For each k-mer, it maps genome names (extracted from the file names) to the number that k-mer occurs in the genome. (I actually reverse the order of the outer and inner to calculate how many kmers in each genome).

Using defaultdict allows for automatic creation of nested dictionaries without needing to check if a key already exists. This structure efficiently groups data by k-mer and genome, enabling quick retrieval of number occurrences.

• How many k-mers are in your index?
Total unique k-mers: 22,104,695
• How many k-mers of length k = 31 are theoretically possible?
Theoretically possible 31-mers: 4,611,686,018,427,387,904
• Explain any discrepancy between these numbers


The actual genomes have a finite length, which means that only a very small subset of the  theoretically possible 31-mers appears in real biological data. Many regions of the genomes are similar or identical, causing repeated k-mers. This redundancy further reduces the number of unique k-mers compared to the astronomical number of all possible combinations. Evolutionary pressures and genomic structure limit the diversity of sequences. Not every theoretical combination is viable or present in nature. While the mathematical possibility of 31-mers is extremely high, the actual indexed k-mers are constrained by the real genomic sequences and their lengths.


### Task 2.2: Implement Classification


• Output the number of matching reads (k-mers) for each organism
For no error reads:
```
E. col:         322,046 (6.94%)

B. subtilis:     53,679 (1.22%)

P. aeruginosa:   51,157 (1.21%)

S. aureus:       51,810 (1.84%)

M. tuberculosis: 55,683 (0.89%)
```
Miseq:
```
E. col:          725,334 (15.63%)

B. subtilis:     123,839 (2.94%)

P. aeruginosa:   114,528 (2.6%)

S. aureus:       123,260 (1.97%)

M. tuberculosis: 119,539 (4.24%)
```
• Is the ratio identical to the approach based on string matching? Explain any discrepancy (if any):

No, it's not. The reason is that the kmers are a lot, so this increases the denominator, but the matches are relitavly few so the numerator doesn't increase much. Especially that we have the number of unique kmers is 99% of the total kmers we have.

• Justify how you handle k-mers that have multiple matches in one genome:

I just count them the number of the times they appear, this shouldn't cause a problem since the unique kmers are 99% of the total number of kmers.
However, for the minimizer we count them only once because we are going to have a lot more matches.


### Task 2.3: Mimimizers

No erro then miseq:

k-mer: time: 39.13s, memory: 9.39GB
       time: 26.61s, memory: 6.66GB
minimizer k-mer: time: 21.68s, memory: 2.69GB
                 time: 20.02s, memory: 1.48GB
reduction in memory by 71.35% and 77.77%


No error:
```
E. col:          263,648 (5.68%) (Accuracy: 5.68%  / 6.94%  = 81.84%)

B. subtilis:     44,107  (0.99%) (Accuracy: 0.99%  / 1.22%  = 81.15%)

P. aeruginosa:   43,648  (1.05%) (Accuracy: 1.05%  / 1.21%  = 86.78%)

S. aureus:       44,513  (1.5%)  (Accuracy: 1.5%   / 1.84%  = 81.52%)

M. tuberculosis: 43,741  (0.70%) (Accuracy: 0.70%  / 0.89%  = 78.65%)
```

Miseq:
```
E. col:          646,036 (13.92%) (Accuracy = 13.92% / 15.63% = 89.06%)

B. subtilis:     106,548 (2.53%)  (Accuracy = 2.53%  / 2.94%  = 86.05%)

P. aeruginosa:   107,686 (2.44%)  (Accuracy = 2.44%  / 2.6%   = 93.85%)

S. aureus:       112,879 (1.8%)   (Accuracy = 1.8%   / 1.97%  = 91.37%)

M. tuberculosis: 105,908 (3.75%)  (Accuracy = 3.75%  / 4.24%  = 88.44%)
```
# 3

## 3.1

We first download the taxanomy
kraken2-build --download-taxonomy --db $DBNAME

Then we add our genomes to the database
kraken2-build --no-masking --add-to-library $GENOME.fna --db $DBNAME

Then we build the database
kraken2-build --build --db $DBNAME

Then we classify the reads using the commnad:
kraken2 --db $DBNAME --report report.txt --output output.txt reads.fasta


Analyze difference in terms of:

• Species abundance estimates
Kraken2:

No error:
```
E. coli:         5983     (59.83%)

B. subtilis:      999     (9.9%)

P. aeruginosa:   1000     (10%)

S. aureus:       1000     (10%)

M. tuberculosis: 1000     (10%)
```
Miseq:
```
E. coli:         5984     (59.84%)

B. subtilis:     1000     (10%)

P. aeruginosa:   1000     (10%)

S. aureus:       1000     (10%)

M. tuberculosis: 1000     (10%)
```
• Classification speed and resource usage
Kraken2: 

no errors: time: 21.86s, memory: 7.62GB
miseq: time: 20.17s, memory: 7.62GB

So Kraken2 uses more memory and takes more time. 
    
• Impact of the minimizer scheme on accuracy
The estimates are far from the results we got in the previous parts.

Advantages and Disadvantages
Kraken2 is easier to run, and we can add a lot of genomes easily. Also the report it generates is very insightful and gives us a clear discreption of the species.

Aho-Corasick uses the least time, while BLAST uses the least memory.

## 3.2

first we download the standard-8 database from:
https://benlangmead.github.io/aws-indexes/k2

then we untar the database using the command:
tar -xzvf k2_standard_08gb_20241228.tar.gz -C standard-8

then we download the metagenomic samples from https://www.ncbi.nlm.nih.gov/sra
To do this, we first install sra-tools using the command:
conda install -c bioconda sra-tools
then we use the command prefetch with the accessions like so:
prefetch SRR11412973

Then we change the format of the file two fasta using the command:
fasterq-dump --split-files SRR11412973

this gives us a pair of files for each sample


last we use kraken2 to classify the reads and make the reports using the command:
kraken2 --db /path/to/kraken2_db/standard-8 \
        --threads 8 \
        --paired \
        --report reports/SRR11412973_report.txt \
        --output results/SRR11412973_output.txt \
        SRR11412973_1.fastq SRR11412973_2.fastq


Then we read the reports to classify the samples
```
SRR11412973: Phascolarctobacterium faecium (6.13%)

SRR11412976: Phocaeicola vulgatus (16.14%)

SRR11412979: Segatella copri (37.91%)

SRR11412980: Bacteroides uniformis (3.94%)

SRR11412984: Segatella copri (16.94%)

SRR21907296: Severe acute respiratory syndrome-related coronavirus (72.16%)

SRR21907303: Severe acute respiratory syndrome-related coronavirus (44.67%)

SRR21907307: Severe acute respiratory syndrome-related coronavirus (41.33%)

SRR21907330: Severe acute respiratory syndrome-related coronavirus (84.25%)

SRR21907332: Severe acute respiratory syndrome-related coronavirus (91.02%)
```

As for the data driven-driven way to seperate the samples, we find the lineage for each species
we see that they all belong to Phascolarctobacterium Genus
but five of them are Severe acute respiratory syndrome-related coronavirus
so that would the first class, and the second class will be the rest

```

SRR11412973: Phascolarctobacterium faecium (6.13%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412976: Phocaeicola vulgatus (16.14%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412979: Segatella copri (37.91%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412980: Bacteroides uniformis (3.94%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412984: Segatella copri (16.94%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907296: Severe acute respiratory syndrome-related coronavirus (72.16%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907303: Severe acute respiratory syndrome-related coronavirus (44.67%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907307: Severe acute respiratory syndrome-related coronavirus (41.33%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907330: Severe acute respiratory syndrome-related coronavirus (84.25%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907332: Severe acute respiratory syndrome-related coronavirus (91.02%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']


We could make this more automated than it already is, by embedding the lineage into a vector space, then doing 2NN clustring algorithm
Also, Phascolarctobacterium is found abundantly in the human gastrointestinal tract so this is consistent with the environments we found them in
```

```

SRR11412973: Phascolarctobacterium faecium (6.13%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412976: Phocaeicola vulgatus (16.14%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412979: Segatella copri (37.91%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412980: Bacteroides uniformis (3.94%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR11412984: Segatella copri (16.94%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907296: Severe acute respiratory syndrome-related coronavirus (72.16%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907303: Severe acute respiratory syndrome-related coronavirus (44.67%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907307: Severe acute respiratory syndrome-related coronavirus (41.33%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']

SRR21907330: Unclassified

SRR21907332: Severe acute respiratory syndrome-related coronavirus (91.02%)
Lineage: ['Phascolarctobacterium', 'Acidaminococcaceae', 'Acidaminococcales', 'Negativicutes', 'Bacillota', 'Bacteria']
```

