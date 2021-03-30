# Developed by: Dionysios Grigoriadis
# Last Updated: September 2020


#############################################################################################################
# This scripts takes an sqlite3 db created by the VCFdbR software and makes it compatible to the PQuery app.
# Just change the db_path to your db and it should then work.
# Use the venv_py3 in PQuery directory to run it.
#############################################################################################################

# Import Libraries
import time 
start = time.time()
import sqlite3
import pandas as pd
import ntpath

# Dependencies
db_path = 'test.db'


# DB connection
db_name = ntpath.basename(db_path)
conn = sqlite3.connect(db_path, detect_types = sqlite3.PARSE_DECLTYPES)
cur = conn.execute('PRAGMA TABLE_INFO(variant_info)')

## If a column has more than 50 distinct values the column will become 'raw' type instead of 'select'. 
mnouv = 50

# variant_info table columns
cols = [x[1] for x in conn.execute('PRAGMA TABLE_INFO(variant_info)')]
df = pd.read_sql("SELECT * FROM variant_info ORDER BY RANDOM() LIMIT 1000",conn)

# Here we try to define a) the data-type of each column b) the filter-type of each column and 
## c) if it's already in the db or if it's going to be calculated on the fly 
## by PQuery when user performs a query: 
## a) The types can be:
## int64: for numerical columns (like chromosome column)
## float64: for numerical columns with float numbers (like cadd_phred_score column)
## text: for columns with text (like func.refGene column)
##
## b) The filter types can be:
## raw: if the column doesn't have unique values (N<50) which can be used in a dropdown multi-select button for example.
## select: the opposite.
##
## c) in_db: The column is already in the table.
##    out_dc: This column is empty and will be calculated on-the-fly
##            everytime the user performs a query on PQuery.

newtypes = {}
for dcol in cols:
    dfcol = df[dcol]
    cdtype = str(dfcol.dtype)
    col_list = set(df[dcol].to_list())
    if cdtype == "int64":
        newtypes[dcol]=["int64"]
    elif cdtype == "float64":
        newtypes[dcol]=["float64"]
    else:
        try:
            dfcol.replace(".",0).astype("float64")
            newtypes[dcol]=["float64"]
        except ValueError:
            newtypes[dcol]=["text"]
    if len(col_list)<mnouv and col_list!="None":
        newtypes[dcol].append("select")
    else:
        newtypes[dcol].append("raw")
    newtypes[dcol].append("in_db")


### These are the only out_db columns. 
### dynamic_ac: Allele Count for each variant on the USER-SELECTED COHORT OF SAMPLES calculated on-the-fly
###             every time the user performs a query.
### dynamic_an: Allele Number: Same as above. 
### dynamic_af: Allele Frequency. Same as above.
newtypes['dynamic_ac'] = ['float64',"raw","out_db"]
newtypes['dynamic_an'] = ['float64',"raw","out_db"]
newtypes['dynamic_af'] = ['float64',"raw","out_db"]

### List created above to dataframe
nt_df = pd.DataFrame.from_dict(newtypes, orient='index')
nt_df.columns = ['data_type','filter_type','db_status']

### Double-check that the select-columns have less unique values than the minimum limited number specified.
select_ids = list(nt_df[nt_df["filter_type"]=='select'].index)
select_to_raw=[]
for x in select_ids:
    counter = [x[0] for x in conn.execute('SELECT COUNT(DISTINCT "'+x+'") FROM variant_info')][0]
    if counter > mnouv:
        select_to_raw.append(x)


### Write nt_df to sql
nt_df.to_sql('processed_field_types',conn)

# Get the select ids and delimit them with this delimiter: "xdio3" (just chose this random string)
select_ids = list(nt_df[nt_df["filter_type"]=='select'].index)
unique_vals_dict = {}
for x in select_ids:
    unique_vals_dict[x] = 'xdio3'.join([str(x[0]) for x in conn.execute('SELECT DISTINCT "'+x+'" FROM variant_info')])

# Create a dataframe with all the 'select' columns and their unique_values 
uvd_df = pd.DataFrame.from_dict(unique_vals_dict, orient='index')
uvd_df.columns = ["unique_values"]
uvd_df.to_sql('select_unique_values',conn)

# Write the column order in a file
current_order = list(nt_df.index)
with open(db_path+'.columns-order.csv', 'a') as the_file:
    for co in current_order:
        the_file.write(co+'\n')

allcolumns = list(nt_df.index)

prexisting_idxs = [x[1] for x in conn.execute('PRAGMA INDEX_LIST(variant_info)')]
prexisting_idxs_cols = []
for pi in prexisting_idxs:
    prexisting_idxs_cols+=[x[2] for x in conn.execute('PRAGMA INDEX_INFO('+pi+')')]

to_index_cols = [x for x in allcolumns if x not in prexisting_idxs_cols+['dynamic_ac','dynamic_an','dynamic_af']]

for tic in to_index_cols:
    try:
        #print('CREATE INDEX idx_info_'+tic.replace('.','_')+' '+"ON variant_info ('"+tic+"')")
        conn.execute('CREATE INDEX idx_info_'+tic.replace('.','_')+' '+"ON variant_info ('"+tic+"')")
    except OperationalError:
        print("Couldn't create index for "+tic)


geno_cols = [x[1] for x in conn.execute('PRAGMA TABLE_INFO(variant_geno)')]

#DB Indexing
conn.execute('CREATE INDEX idx_geno_2'+' '+"ON variant_geno ("+geno_cols[0]+", "+geno_cols[2]+")")
conn.execute('CREATE INDEX idx_geno_3'+' '+"ON variant_geno ("+geno_cols[1]+", "+geno_cols[0]+", "+geno_cols[2]+")")

end = time.time()
print(end - start)
