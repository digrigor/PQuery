# <img src="app/static/pquery_logo-03.png" width="400">

PQuery is a Flask-based web application enabling fast interactive visualisation, querying, flexible filtering and alignment inspection of multi-sample, annotated, genetic variant data. 
  
Without programming skills required, users can automatically convert their input VCF file (compressed or uncompressed, annotated or not annotated) to an indexed and structured SQLite database, fully adjusted to the imported VCF file. Once the database is created and stored, PQuery can repeatedly open in a browser window and instantly run on it with no need for re-importing every time it runs. The application allows users to query the genetic variants for specific samples and genomic regions and it then instantly returns a table with the queried variants (one variant per row) and all fields present in the original VCF file. It also dynamically calculates the allele counts and frequencies for the selected samples.
  
PQuery can run in two modes:
a. The 'cohort mode' in which no sample-specific information is included in the returned table
b. the 'sample-specific' mode in which the generated table will fully include the sample names and all the related information of the original VCF file.
  
Users can efficiently filter all columns in both modes by using a customisable filter pane and inspect the variant calls through the integrated IGV viewer.

# Getting Started

## Execute PQuery
- From a linux machine:
```
sh run_app.sh
```

- From a Windows machine:  
Double-click the PQuery.bat file or the PQuery.lnk shortcup (you can copy the shortcut to your Desktop).

## How to Guide
Read more in the [PQuery Manual](./PQuery_Manual.pdf)

# Preprocess an annotated VCF file so you can use in PQuery
Check the preprocessing directory.

# Further Information
- PQuery is currently operating on POlab hg19 data.

To do:
- Update the DB with the latest hg38 data.
- Update the back-end DB so can be faster. This scheme used right now (VCFdbR) is currently slow for querying whole exome for many samples. The major limitation is the time needed to calculate on-the-fly AC,AN and AF (Allele Count, Allele Number, Allele Frequency) for the queried cohort of samples and then filter these numbers to only include variants which have at least ONE alternative allele in at least ONE queried sample-set. Using NOSQL db management or GATK GenomicDBI might be good solutions.

# PQuery Architecture
