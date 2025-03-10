# Human genome task

## HPC
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip hg38.fa.gz

for i in {1..4}; do
    nohup python -u task$i.py > task$i.log 2>&1 &
done
```

## Python scripts

### Task 1

```python
import math

# Reads a FASTA file and returns a dictionary
def read_fasta(filename):
    genome = {}
    chrom = None
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chrom = line[1:].split()[0]
                genome[chrom] = []  
            elif chrom is not None:
                genome[chrom].append(line)
    return {chrom: ''.join(seq) for chrom, seq in genome.items()}  

# Splits the genome into fixed-size bins and calculates GC content
def split_genome_into_bins(genome, num_bins=1000000):
    gc_results = []
    for chrom, seq in genome.items():
        bin_size = math.ceil(len(seq) / num_bins)
        for i in range(0, len(seq), bin_size):
            bin_seq = seq[i:i + bin_size]
            gc_content = calculate_gc_content(bin_seq)
            gc_results.append((chrom, i, min(i + bin_size, len(seq)), gc_content))  
    return gc_results

# Calculates the GC content (%) of a given DNA sequence
def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in "GC")  
    return (gc_count / len(sequence) * 100) if sequence else 0  

# Counts the number of bins with GC content <30% and >60%
def count_gc_thresholds(gc_results, low=30, high=60):
    low_count = sum(gc < low for _, _, _, gc in gc_results)
    high_count = sum(gc > high for _, _, _, gc in gc_results)
    return low_count, high_count

# Writes GC content results to a BED format file
def write_bed_file(gc_results, output_file):
    with open(output_file, 'w') as f:
        f.write("chrom\tstart\tend\tGC_content\n")
        f.writelines(f"{chrom}\t{start}\t{end}\t{gc:.2f}\n" for chrom, start, end, gc in gc_results)  

# File paths
base_path = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/"
fasta_file = f"{base_path}hg38.fa"
output_bed = f"{base_path}gc_content.bed"
output_txt = f"{base_path}gc_content_stats.txt"

# Process genome
genome = read_fasta(fasta_file)
gc_results = split_genome_into_bins(genome)
write_bed_file(gc_results, output_bed)
low_gc, high_gc = count_gc_thresholds(gc_results)

# Writing GC threshold results to file
with open(output_txt, 'w') as f:
    f.write(f"GC content <30% bins: {low_gc}\nGC content >60% bins: {high_gc}\n")

print("Processing complete. Results saved in:")
print(f"  - BED file: {output_bed}")
print(f"  - GC statistics: {output_txt}")
```

### Task 2

```python
# Returns the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1] 

# Reads the genome and generates the reverse complement sequence
def generate_reverse_complement(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        chrom = None
        seq_lines = []
        for line in infile:
            if line.startswith(">"):
                if chrom:
                    outfile.write(f">{chrom}_reverse_complement\n")
                    outfile.write(reverse_complement("".join(seq_lines)) + "\n")
                chrom = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if chrom:
            outfile.write(f">{chrom}_reverse_complement\n")
            outfile.write(reverse_complement("".join(seq_lines)) + "\n")


# File paths
fasta_file = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/hg38.fa"
output_fasta = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/hg38_reverse_complement.fa"

# Generate reverse complement
generate_reverse_complement(fasta_file, output_fasta)

print(f"Reverse complement genome saved to: {output_fasta}")
```

### Task 3

```python
from collections import Counter

# Counts the frequency of 4-mer motifs
def count_4mer_frequencies(input_fasta, output_txt):
    kmer_counts = Counter()
    
    with open(input_fasta, 'r') as infile:
        chrom = None
        seq_lines = []
        for line in infile:
            if line.startswith(">"):  
                if chrom:  
                    genome_seq = "".join(seq_lines)
                    for i in range(len(genome_seq) - 3):
                        kmer = genome_seq[i:i+4] 
                        kmer_counts[kmer] += 1
                chrom = line[1:].strip()
                seq_lines = []  
            else:
                seq_lines.append(line.strip())
        if chrom:
            genome_seq = "".join(seq_lines)
            for i in range(len(genome_seq) - 3):
                kmer = genome_seq[i:i+4]
                kmer_counts[kmer] += 1
    
    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
    
    with open(output_txt, 'w') as outfile:
        outfile.write("4-mer\tFrequency\n")
        for kmer, count in sorted_kmers:
            outfile.write(f"{kmer}\t{count}\n")

# File paths
fasta_file = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/hg38.fa"
output_txt = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/4mer_frequencies.txt"

# Count 4-mer frequencies
count_4mer_frequencies(fasta_file, output_txt)

print(f"4-mer frequency analysis saved to: {output_txt}")
```

### Task 4

```python
# Counts CG dinucleotides and calculates their percentage per chromosome
def count_CG_dinucleotides(input_fasta, output_txt):
    total_dinucleotides = 0
    CG_count = 0
    
    with open(input_fasta, 'r') as infile:
        chrom = None
        seq_lines = []
        for line in infile:
            if line.startswith(">"):  
                if chrom:  
                    genome_seq = "".join(seq_lines)
                    for i in range(len(genome_seq) - 1):
                        dinucleotide = genome_seq[i:i+2] 
                        if dinucleotide == "CG":
                            CG_count += 1
                        total_dinucleotides += 1
                chrom = line[1:].strip()
                seq_lines = []  
            else:
                seq_lines.append(line.strip())
        
        if chrom:
            genome_seq = "".join(seq_lines)
            for i in range(len(genome_seq) - 1):
                dinucleotide = genome_seq[i:i+2]
                if dinucleotide == "CG":
                    CG_count += 1
                total_dinucleotides += 1
    
    CG_percentage = (CG_count / total_dinucleotides) * 100 if total_dinucleotides > 0 else 0
    
    with open(output_txt, 'w') as outfile:
        outfile.write(f"Total dinucleotides: {total_dinucleotides}\n")
        outfile.write(f"CG dinucleotide count: {CG_count}\n")
        outfile.write(f"CG percentage: {CG_percentage:.2f}%\n")

# File paths
fasta_file = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/hg38.fa"
output_txt = "/home/xh368/rds/rds-huang_xr-CGClDViiOBk/cuhk/CG_dinucleotide_stats.txt"

# Count CG dinucleotides
count_CG_dinucleotides(fasta_file, output_txt)

print(f"CG dinucleotide statistics saved to: {output_txt}")
```
