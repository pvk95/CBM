# GroupB aligner

Aligner created by GroupB Computational Biomedicine students.

**Authors:** Brynja Sigurpalsdottir, Hana Parizkova, Karthik Pattisapu, Modestas Filipavicius, Zhi Ye

### Strategies used
- **Genome Indexing:** Burrows-Wheeler transform
- **Seeding:** FM-index; variable length and number of seeds; position of seed in the read chosen randomly
- **Alignment:** Semi-global dynamic-programming alingment with affine gap penalties; scoring as in bowtie2

### Usage
**Requires Python version 3.6.6 or newer.**

Run `run_alignment.py`.

#### Required arguments
- `--fastq` ... names of one or more FastQ files with reads to be aligned
- `--genomeIndexDir` ... directory where genome index is stored, or where it should be stored if we run indexing as well
- `--outFile` ... output file (in SAM format)
- `--genomeFasta` ... name of Fasta file with the genome sequence

#### Other basic arguments
- `--runIndexing` ... run also indexing of the genome; by default indexing is **not** run

#### Advanced arguments
- `--overhang` ... float value specifiying how large overhang over the length of the read will be taken into account when running DP alignment; default value `0.33`, i.e. the alignment will be computed on a piece of genome of lenght `1.33*readLength`
- `--seedLength` ... length of seeds in base pairs; default `10`
- `--numberOfSeeds` ... number of seeds to be tried for each read; default value `3`

#### Examples
- Basic usage, genome index already present: `python run_alignment.py --fastq FILE.FASTQ --genomeIndexDir GENOMEDIR --genomeFasta GENOME.FASTA --outFile OUTFILE`
- If we want to run indexing of genome as well: `python run_alignment.py --fastq FILE.FASTQ --genomeIndexDir GENOMEDIR --genomeFasta GENOME.FASTA --outFile OUTFILE --runIndexing`
