This directory has all the necessary files you need to take an annotated multi-sample VCF 
convert it into an sqlite3 db and make the db compatible with PQuery

Follow these steps:

1. Run the annotated VCF file through VCFdbR to convert it to an sqlite3 db. 
   You can find VCFdbR software in this directory.
   VCFdbR information can be found here: https://github.com/tkoomar/VCFdbR

2. Get the output sqlite3 database and use the vcfdbR_to_PQuery_db.py to make it PQuery compatible.
