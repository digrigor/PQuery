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
