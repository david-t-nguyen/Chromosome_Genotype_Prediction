# READ ME FILE

Call script using the line
```
   python3 snps_caller.py align_sort.bam putatative_snps.tsv
```

Note:
- .bam file is the file with the alignments and data
- .tsv file is the metadata file 

Format to use for other files would be 
```
   python3 snps_caller.py <.bam file>  <.tsv file>
```

3 files will be created
## .bai file 
- indexed bam file 

## output.tsv 
- contains the posterior probability of the most probable genotype for each chromosome

## all_poss_genotypes_output.tsv 
- contains the posteroir probability of all possible genotypes for each chromosome



At the end, we can the following line to get remove output files and clean the directory 
for future uses of snps_caller.py
```
python3 clean.py
```



