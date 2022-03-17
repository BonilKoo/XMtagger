# Bisulfite XM tagger
Python code for adding XM tag (methylation call string) to bisulfite bam file.

Both Single-end and paired-end are available.

### Motivation
Some bisulfite alignment tools do not provide XM tag. XM tag is required to compute levels of DNA methylation heterogeneity.

### Package requirement
+ pysam (https://pysam.readthedocs.io/en/latest/index.html)

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

## Additional Information
For faster results, we recommend using the tag command in Metheor implemented in the Rust programming language.
https://github.com/dohlee/metheor#miscellaneous
