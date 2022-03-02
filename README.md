# Bisulfite XM tagger
Python code for adding XM tag to bisulfite bam file.

Both Single-end and paired-end are available.

### Package requirement
+ pysam

### File requirement
The input bam file must be indexed.

## Usage
```
python XMtagger.py [reference genome file] [input bam file] [output bam file]
```

### Example
```
python XMtagger.py hg38.fa input.bam output.bam
```
