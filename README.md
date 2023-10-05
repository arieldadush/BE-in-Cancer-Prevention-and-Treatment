# BE-in-Cancer-Prevention-and-Treatment
While existing research primarily focuses on applying BE to treat inherited diseases caused by single-nucleotide variants (SNVs), it may seem that BE has limited relevance in the field of oncology, given the multitude of mutations present in cancer cells. Nevertheless, we demonstrate that BE holds promise in the treatment and prevention of cancer across various tumor types.
# Genomic annotations
We assume you have the following files:
1. ```hg19.fa```
2. ```hg19_wgEncodeGencodeBasicV44lift37.fa``` 

```hg38_ncbiRefSeqCurated_cds.fa```
3. 
5. ```GRCh38_latest_protein_shortHeders.faa```
6. ```RefSeq Curated.gtf```

# Installation and Requirements
## Dependencies

[Python 3.10.4](https://www.python.org/downloads/release/python-3104/)

[BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat)

## Human cancer mutations data
We use these database:
1. [Pediatric Cancer Panel](https://www.ncbi.nlm.nih.gov/gtr/tests/562503/)
2. [Table S1 of "The majority of genetic point mutations can be reinstated by DNA and RNA base editors" paper](https://github.com/arieldadush/BE-on-genetic-point-mutations.git)
3. [MedGen](https://www.ncbi.nlm.nih.gov/medgen/)
4. [MSK-IMPACT](https://www.nature.com/articles/s41588-021-00949-1)
5. [PCAWG project](https://www.nature.com/articles/s41586-020-1969-6)


**Pipline on MSK-IMPACT**

```MSK_database_off_target.ipynb```

After downloading the database, we downloaded from [Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1723543828_rv8jAZ6jYoeMPtpo2d5W7THrMVa9&clade=mammal&org=&db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeqCurated&hgta_regionType=range&position=&hgta_outputType=primaryTable&hgta_outFileName=) a list of RefSeq genes from NCBI (in hg19 annotation).
We used this list to add a strand to the mutation and arranged the database to create a BED6 format.

```add_20_dna_seq_cancer_mutatuon.sh```
Adding a DNA sequence to the mutation - 20 bases before and 20 bases after
