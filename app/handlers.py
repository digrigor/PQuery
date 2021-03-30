from flask import render_template, request, redirect
from sqlalchemy import distinct, func, and_, or_, cast, Numeric, case, String, Float
from datatables import ColumnDT, DataTables
import pandas as pd
import os
import time
from io import TextIOWrapper
import sys
from optparse import OptionParser
import csv
#import paramiko
import re
import base64
import datetime
import io
import sqlite3
import functools
from app.dependencies import *
#import time
#import sys
#from optparse import OptionParser
#import os
#import sys
#import csv
##import paramiko
#import subprocess
#import re
#import getpass
#import labkey

def parse_contents(fileobject):
    """Takes a file object as its input and it outputs the lines of the input file in a list."""
    csv_file = TextIOWrapper(fileobject, encoding='utf-8')
    csv_reader = csv.reader(csv_file, delimiter=',')
    clines=[]
    for row in csv_reader:
        try:
            clines.append(row[0])
        except IndexError:
            pass
    return(clines)

flatten = lambda l: [item for sublist in l for item in sublist]

def csvlines2list(csv_path):
    """Function that takes a path for a csv file, it opens it and returns each line as a
    python list element."""
    csv_in = csv_path.strip()
    with open(csv_in, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        outlist=[]
        for row in csv_reader:
            row = [x.strip() for x in row]
            row = [x.split("\xef\xbb\xbf")[1] if x.startswith("\xef\xbb\xbf") else x for x in row]
            outlist+=row
    return(outlist)


def string_none(in_option):
    if in_option==None:
        out_option=["None"]
        return (out_option)
    else:
        return(in_option)

def recin(instring):
    """Function that takes a string as its input and returns a list of boolean values indicating
    what type this input is"""
    recdict={}
    instring = instring.strip()
    if bool(re.compile(repattern).match(instring)): otype="loc"
    elif bool((re.compile("ENSG\d{11}").match(instring))): otype="ens"
    elif (instring.endswith(".csv") | instring.endswith(".txt")): otype="csv"
    elif instring=="all": otype="all"
    else:
        otype="gene"
    return(otype)

def genes2locs(genes, g2f):
    """Function that takes a list of genes or ensembl identifiers as its input and returns
    its corresponding coordinates"""
    flist=[]
    for gen in genes:
        gen = gen.strip()
        if recin(gen)=="ens":
            flist+=g2f[g2f["ensembl_gene_id"]==gen][["chromosome_name","start_position","end_position"]].values.tolist()
        elif recin(gen)=="gene":
            flist+=g2f[g2f["hgnc_symbol"]==gen][["chromosome_name","start_position","end_position"]].values.tolist()
    flist = [list(x) for x in set(tuple(x) for x in flist)]
    return(flist)

def check_coords(coord):
    """Function that takes coordinates and a regex pattern string as its input,
       checks if the coordinate matches the coordinate input criteria and
       returns a list of the genomic coordinates"""
    #Compile the checking pattern
    pattern=re.compile(repattern)
    #Check if given coordinates follow the pattern
    correct_pattern = bool(pattern.match(coord))
    #Accept or complain
    if correct_pattern == True:
        pass
    else:
        print("ERROR: Sorry. Please check your genomic coordinates.")
        print("       This program only accepts genomic coordinates")
        print("       of this form:")
        print("       chrN:XXXXXX-XXXXXX")
        print("       N:XXXXXX-XXXXXX")
        print("       chrN:XXXXXX")
        print("       N:XXXXXX")
        print("       chrN")
        print("       N")
        raise
    #Split the coordinates (: and -)
    splcoord = re.split('; |, |\:|\-',coord)
    splcoord = [x for x in splcoord if x!=""]
    #Get the three variables
    c1=splcoord[0]
    c1=c1.split("chr")[-1]
    if len(splcoord)>1:
        c2=splcoord[1]
    else:
        c2=""
    if len(splcoord)>2:
        c3=splcoord[2]
        if(int(c2)>int(c3)):
            print("Please check your range. The end position is before the start position")
            raise
    else:
        c3=""
    coordlist=[c1,c2,c3]
    return(coordlist)

def return_status(in_list, many):
    """Function that takes a list as its input and a boolean variable and returns NoneError if the len(list)==0,
    ConflictError if len(in_list)>=2 and many==True or OK if nothing of the above applies.
    This will help handling issues when the user uses more than of the proposed ways to input samples or genomic coordinates."""
    if len(in_list)==0:
        return("NoneError")
    elif len(in_list)>=2 and many==True:
        return("ConflictError")
    else:
        return("OK")

def handle_input_error_messages(all_samples, coords):
    """Function that takes two lists (each of them is direct output of the return_status function) and returns a list with error messages"""
    error_list = []
    if all_samples == "ConflictError":
        error_list.append(
            "Please select samples by using only ONE of the given options (Select Samples OR Select Cohort OR Upload CSV)")
    if all_samples == "NoneError":
        error_list.append(
            "No sample selected. Please select samples to proceed. To select all samples just select All in the Select Samples dropdown menu.")
    if coords == "ConflictError":
        error_list.append(
            "Please select Genomic Location(s) by using only ONE of the given options (Select Genomic Location OR Enter Gene or Location OR Upload CSV)")
    if coords == "NoneError":
        error_list.append(
            "No Genomic Location(s) selected. To select the Whole Exome sequencing regions just select All in the Select Genomic Location dropdown menu.")
    return(error_list)



def handle_filter_req(req, reqf):
    """Function that takes a request and a request.files POST OBJECT and returns 2 lists with the selection of the user
    for a)samples b)genes and a dictionary with c)filtering_preferences(not currently used)."""
    in_samples = [y for y in [req.getlist(key) for key in req if key.startswith("samples")] if y != ['']]
    in_genes = [y for y in [req.getlist(key) for key in req if key.startswith("genes")] if y != ['']]
    in_filt = [y[0] for y in [req.getlist(key) for key in req if key.startswith("filt-")] if y != ['']]
    in_filt_names = [key for key in req if key.startswith("filt-")]
    in_filt_dict = dict(zip(in_filt_names,in_filt))
    if reqf["samples-upload"].filename!="":
        upsamples = parse_contents(reqf["samples-upload"])
        if len(upsamples)!=0:
            in_samples.append(upsamples)
    if reqf["genes-upload"].filename!="":
        upgenes = parse_contents(reqf["genes-upload"])
        if len(upgenes)!=0:
            in_genes.append(upgenes)
    in_samples = [x for x in in_samples if (x != "" or x != [""])]
    in_genes = [x for x in in_genes if (x != "" or x != [""])]
    return(in_samples,in_genes,in_filt_dict)

def handle_genes_input(in_genes):
    """Function that takes a list with the user genes input and returns a list with genomic coordinates corresponding to the output"""
    in_genes = flatten(in_genes)
    in_genes = flatten([x.split(",") for x in in_genes])
    if "all" not in in_genes:
        coords=[]
        coords += genes2locs([x for x in in_genes if recin(x) in ["ens", "gene"]], g2f)
        coords += [check_coords(x) for x in in_genes if recin(x) == "loc"]
        coords = [list(y) for y in set([tuple(x) for x in coords])]
    else:
        coords = "all"
    return(coords)

def handle_samples_input(in_samples):
    """Function that takes a list with the user samples input. Iterates through the input list and creates a list with the to-be-queried samples.
    Through this iteration the function also checks if any of the input words correspond to the specific keywords (i.e. lymphoedema) and if yes
    then adds the corresponding samples of the cohort to the output."""
    in_samples = flatten(in_samples)
    new_samples = []
    for x in in_samples:
        if x in coh2file_dict.keys():
            new_samples+=csvlines2list(coh2file_dict[x])
        else:
            new_samples.append(x)
    return(new_samples)


def create_connection(db_file=vardb):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn

def db_colnames(conn):
    dbcolnames = []
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(var);")
    rows = cur.fetchall()
    for row in rows:
        dbcolnames.append(row[1])
    return(dbcolnames)

def db_serverschema(in_cols):
    counter = 1
    server_list = []
    for dbc in in_cols:
        tempdict = {"data_name": dbc,
                    "column_name": dbc,
                    "default": "",
                    "order": counter,
                    "searchable": True}
        server_list.append(tempdict)
        counter+=1
    return(server_list)

def df_dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def get_yadcf_type(db_col, db_sams, field_type_dict):
    """Function that takes as its input:
    a) db_col: String: Column ID/Name
    b) db_sams: List: List of Samples.
    c) Field_type_dict: output of field_type_dict_builder(ftquery)
    And returns a string showing the filtering type of the input column name:(multi-select,float,int,text or sample)"""
    if db_col in field_type_dict.keys():
        if field_type_dict[db_col]['filter_type']=='select':
            return('multi_select')
        else:
            if field_type_dict[db_col]['data_type'] in ['float64','int64']:
                return('range_number')
            if field_type_dict[db_col]['data_type']=='text':
                return('text')
    elif db_col in db_sams:
        return('sample')
    else:
        return(None)

def column_dt_builder(main_cols, db_sams, Alleles, subquery, mode="nosamples"):
    """Function that takes as its input:
    a) main_cols: List of samples which will be returned to the var_Table page
    b) db_sams: List, List of Sample ids.
    c) Alleles: Model of genotypes sqlite3 table
    e) subquery: Sqlite3 query
    f) mode: String: "samples" or "nosamples"

    And returns the list of ColumnDT columns which will be returned by the main query.
    """
    dt_columns = []
    if mode=="samples":
        all_cols = main_cols + db_sams
    else:
        all_cols = main_cols
    for ac in all_cols:
        if ac in main_cols:
            if ac == "dynamic_ac": dt_columns.append(ColumnDT(func.sum(Alleles.gt).label(ac)))
            elif ac == "dynamic_an": dt_columns.append(ColumnDT((func.count(Alleles.gt)*2).label(ac)))
            elif ac == "dynamic_af": dt_columns.append(ColumnDT((func.sum(Alleles.gt)/(func.count(Alleles.gt)*2)).label(ac)))
            else:
                dt_columns.append(ColumnDT(subquery.columns[ac].label(ac)))
        else:
            samcol = func.max(case([(Alleles.sample == ac, Alleles.gt_raw+":"+Alleles.sample),],))
            dt_columns.append(ColumnDT(samcol.label(ac)))
    return(dt_columns)



def select_dict_builder(svquery):
    """Function that takes the output of the select_unique_values table of the sqlite database
    and returns a dictionary with:
    INFO column names a keys. Each key points to a mini dictionary containing a 'value' key which
    corresponds to a unique value for this field and a 'label' key which corresponds to the
    label for each unique value (the one that will appear in the drop-down buttons)."""
    select_values = {}
    for f in svquery:
        dis_vals = f.unique_values.split(unique_delimiter)
        select_values[f.index] = [{"value":x, "label":x.replace('\\x3b', '-')} for x in dis_vals]
    return(select_values)

def field_type_dict_builder(ftquery):
    """Function that takes the output of the processed_field_types table of the sqlite database
    and returns a dictionary with:
    INFO column names a keys. Each key points to a mini dictionary with sqlite datatype (float,text,etc),
    fliter type (multiselect, raw, etc.) and if it's already in the db or it has to be
    calculated on the fly when executing the query """
    field_types = {}
    for f in ftquery:
        field_types[f.index] = {'data_type': f.data_type,
                                'filter_type': f.filter_type,
                                'db_status': f.db_status}
    return(field_types)


def filter_ldict_builder(main_cols, db_sams, field_type_dict, select_dict):
    """Function that takes these inputs:
    main_cols: List with the INFO column names.
    db_sams: List with selected sample names.
    field_type_dict: output of field_type_dict_builder(ftquery)

    And returns a list of diotionaries as its output. Each dictionary
    corresponds to a value from the main_cols + db_sams columns and it has 4 keys:
    id -> name of the input column
    value -> position in the final table
    type -> column type (one of: multi_select,range_number,data_type,text,sample)
    div_id -> unique_id for that column"""
    counter = 0
    filter_ldict = []
    for col in (main_cols+db_sams):
        temp_dict={}
        temp_dict['id'] = col
        temp_dict['value'] = counter
        temp_dict['type'] = get_yadcf_type(col, db_sams, field_type_dict)
        temp_dict['div_id'] = 'data_'+str(counter)
        if temp_dict['type'] == 'multi_select':
            dis_vals = select_dict[col]
            temp_dict['data'] = dis_vals
        filter_ldict.append(temp_dict)
        #print(temp_dict)
        counter+=1
    return(filter_ldict)

def initial_query_builder(inquery, coords, Var):
    """Function that takes: a) inquery: an initial query of the info table
    b) coords: user selected list of coords as its input (this is the output of handle_genes_input() function
    c) Var: the model of the sqlite3 table with the info columns
    It then builds a query by applying a filter on the inquery (filters the table based on the input genomic coordinates)
    And returns this query as its output."""
    tempquery = inquery
    conditions = []
    if coords!='all':
        for coord in coords:
            if coord[1]=='':
                conditions.append(and_(Var.chr == str(coord[0])))
            elif coord[2]!='':
                conditions.append(and_(Var.chr == str(coord[0]),Var.start >= str(coord[1]),Var.end <= str(coord[2])))
            else:
                conditions.append(and_(Var.chr == str(coord[0]), Var.start == str(coord[1])))
        tempquery = inquery.filter(or_(*conditions))
        outquery = tempquery
    else:
        outquery = inquery
    return(outquery)


def fpane_query_builder_var(varquery, fpane_dict, db_sams, Var, field_type_dict):
    """Function that takes as its input:
    a) varquery: initial query (usually output of the initial_query_builder)
    b) fpane_dict: dictionary of req object with user filter preferences.
    c) db_sams: list of sample ids
    d) Var: the model of the sqlite3 table with the info columns
    e) field_type_dict: output of field_type_dict_builder(ftquery).

    It then processes the fpane_dict dictionary,
    checks if the user filters are not on samples (if sample it will be processed by the fpane_query_builder_having below
    and calls filters_builder function
    that performs each specified filter sequentially.

    It then outputs the final filtered query for the info columns.
    """
    tempquery = varquery
    for key, value in fpane_dict.items():
        keyid = key.split('_from')[0].split('_to')[0]
        if keyid not in db_sams and field_type_dict[keyid]['db_status'] == 'in_db':
            filterable = Var.__table__.c[keyid]
            tempquery = filters_builder(tempquery, type='where', dkey=key, dvalue=value, db_sams=db_sams, field_type_dict=field_type_dict, filterable=filterable)
        else:
            pass
    return(tempquery)

def fpane_query_builder_having(mainquery, fpane_dict, db_sams, Alleles, field_type_dict):
    """Function that takes as its input:
    a) mainquery: initial query of the ALleles table.
    b) fpane_dict: dictionary of req object with user filter preferences.
    c) db_sams: list of sample ids
    d) Alleles: the model of the sqlite3 table with the per variant per sample genotypes
    e) field_type_dict: output of field_type_dict_builder(ftquery).

    It then processes the fpane_dict dictionary,
    and checks if the user filters are on samples
    or on the dynamically calculated info columns (dynamic ac,af and an)
    and calls filters_builder function
    that performs each specified filter sequentially (With sqlite3 HAVING statement).

    It then outputs the final filtered query for the info columns.
    """
    tempquery = mainquery
    for key, value in fpane_dict.items():
        print(key)
        keyid = key.split('_from')[0].split('_to')[0]
        print(keyid)
        if keyid in db_sams:
            filterable = func.max(case([(Alleles.sample == keyid, Alleles.gt_raw+":"+Alleles.sample),],))
            tempquery = filters_builder(tempquery, type='having', dkey=key, dvalue=value, db_sams=db_sams,
                                        field_type_dict=field_type_dict, filterable=filterable)
        elif keyid not in db_sams and field_type_dict[keyid]['db_status'] == 'in_db':
            pass
        else:
            if keyid == 'dynamic_ac':
                filterable = func.sum(Alleles.gt)
                tempquery = filters_builder(tempquery, type='having', dkey=key, dvalue=value, db_sams=db_sams,
                                            field_type_dict=field_type_dict, filterable=filterable)
            elif keyid == 'dynamic_an':
                filterable = func.count(Alleles.gt)*2
                tempquery = filters_builder(tempquery, type='having', dkey=key, dvalue=value, db_sams=db_sams,
                                            field_type_dict=field_type_dict, filterable=filterable)
                filterable = func.sum(Alleles.gt)/(func.count(Alleles.gt)*2)
                tempquery = filters_builder(tempquery, type='having', dkey=key, dvalue=value, db_sams=db_sams,
                                            field_type_dict=field_type_dict, filterable=filterable)
    return(tempquery)



def filters_builder(tempquery,type,dkey,dvalue,db_sams,field_type_dict,filterable):
    """Function that takes as its input:
    a) tempquery: sqlite3 query
    and adds filters to it based on the other inputs:
    b) type: String: Can be 'where' for adding a where statement or 'having' for adding a having statement.
    c) dkey: String: Column name (either info or sample) If range column then will ends with "_from" or "_to"
    d) dvalue: List: Filter value(s).
    e) dbsams: List: List of samples.
    f) field_type_dict: output of field_type_dict_builder(ftquery).
    e) filterable: Column which will be queried in sqlite format: i.e. Var.__table__.c[keyid]
    and returns the filtered query as its output
    """
    #If the to-be-applied filter corresponds to a column which is of type: range (from x to y):
    if dkey.endswith("_from") or dkey.endswith("_to"):
        if dkey.endswith("_from") and type=='where':
            dkeyid = dkey.split("_from")[0]
            tempquery = tempquery.filter(cast(filterable, Float) >= str(dvalue[0]))
        elif dkey.endswith("_from") and type=='having':
            dkeyid = dkey.split("_from")[0]
            tempquery = tempquery.having(filterable >= str(dvalue[0]))
        if dkey.endswith("_to") and type == 'where':
            dkeyid = dkey.split("_from")[0]
            tempquery = tempquery.filter(cast(filterable, Float) <= str(dvalue[0]))
        if dkey.endswith("_to") and type == 'having':
            dkeyid = dkey.split("_from")[0]
            tempquery = tempquery.having(filterable <= str(dvalue[0]))
    # If the to-be-applied filter corresponds to a column which is of type: multiselect:
    elif get_yadcf_type(dkey, db_sams, field_type_dict) == "multi_select":
        tempquery = tempquery.filter(filterable.in_([str(x) for x in dvalue]))
    # If the to-be-applied filter corresponds to a column which is a sample:
    elif get_yadcf_type(dkey, db_sams, field_type_dict) == "sample":
        if dvalue[0] == "refhom":
            tempquery = tempquery.having(filterable.startswith("0/0"))
        elif dvalue[0] == "althom":
            tempquery = tempquery.having(filterable.startswith("1/1"))
        elif dvalue[0] == "het":
            tempquery = tempquery.having(or_(filterable.startswith("0/1"), filterable.startswith("1/0")))
        elif dvalue[0] == "car":
            tempquery = tempquery.having(or_(filterable.startswith("0/1"),
                                             filterable.startswith("1/0"),
                                             filterable.startswith("1/1")))
    # If the to-be-applied filter corresponds to a column which is of type: free-form text:
    elif get_yadcf_type(dkey, db_sams, field_type_dict)=="text":
        tempquery = tempquery.filter(filterable == dvalue[0])
    outquery = tempquery
    return(outquery)






































###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
#######################################                                                      ##############################################
#######################################   UNUSED FUNCTIONS - COULD BE USEFUL IN THE FUTURE   ##############################################
#######################################                         :)                           ##############################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################



#UNUSED FUNCTIONS - COULD BE USEFUL IN THE FUTURE
def pd_intervals_subset(df, coords):
    ranges = coords
    index_flist = []
    for rg in ranges:
        nrg = ["", "", ""]
        nrg[0] = rg[0]
        if rg[1] == "":
            nrg[1] = 0
        else:
            nrg[1] = rg[1]
        if rg[2] == "":
            nrg[2] = 100000000000000
        else:
            nrg[2] = rg[2]
        index_flist = index_flist + list(df[(df['CHROM'] == nrg[0]) & (df['POS'] <= max([nrg[1], nrg[2]])) & (
                    df['POS'] >= min([nrg[1], nrg[2]]))].index.values)
    cdf = df.loc[index_flist, :]
    return(cdf)


def slider_handler(rcolname, slider_value):
    sv = int(slider_value)
    conv = gnomad_filt_conv_dict[slider_value]
    return [rcolname,0,conv]

def filter_subset_df(df, grange):
    fdf = df[(df[grange[0]]>=grange[1]) & (df[grange[0]]<=grange[2])]
    return(fdf)

def subset_df(df, coords, all_samples, gnex, gnge, sgaf):
    ##filters_subset
    ###gnomad_exome
    if gnex[1:3] == [0,1]:
        f1df = df
    else:
        f1df = filter_subset_df(df, gnex)
    ##gnomad_genome
    if gnge[1:3] == [0,1]:
        f2df = f1df
    else:
        f2df = filter_subset_df(df,gnge)
    ##sgul_aff
    if sgaf[1:3] == [0,1]:
        f3df = f2df
    else:
        f3df = filter_subset_df(df,sgaf)
    ##coords subset
    if coords=="all":
        cdf = f3df
    else:
        cdf = pd_intervals_subset(f3df,coords)
    ##samples_subset
    if all_samples=="all":
        asdf = cdf
    else:
        cols = list(df.columns)
        main_cols = cols[0:cols.index("FORMAT")+1]
        sams = cols[cols.index("FORMAT")+1:len(cols)]
        new_sams = [x for x in sams if x in all_samples]
        new_cols = main_cols+new_sams
        asdf = cdf[new_cols]
    return(asdf)


    ##samples subset
def split_filter_part(filter_part):
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find('{') + 1: name_part.rfind('}')]
                value_part = value_part.strip()
                v0 = value_part[0]
                if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
                    value = value_part[1: -1].replace('\\' + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part
                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value
    return [None] * 3


def main_raw_query_builder(v_info, v_geno, coords, samples, info_cols, format_cols, table_name):
    ac = "sum("+v_geno+".gt) AS dynamic_ac"
    an = "count("+v_geno+".gt)*2 AS dynamic_an"
    af = "sum("+v_geno+".gt)/count("+v_geno+".gt)*2 AS dynamic_af"
    if samples!="all":
        leftjoins, sampleindi = samples2leftjoins(v_geno, samples, format_cols, info_cols)
        sampleindi2=","+sampleindi
    else:
        leftjoins=''
        sampleindi2=''
    query = "CREATE TABLE "+table_name+" "+"AS SELECT "+v_geno+"."+info_cols[0]+","+\
            ",".join(["t.'"+x+"'" for x in info_cols[1:]])+","+ac+","+an+","+af+' '+sampleindi2+\
            "FROM "+v_geno+" "+\
            "JOIN "+"(SELECT * FROM "+v_info+") t ON "+v_geno+"."+info_cols[0]+"=t."+info_cols[0]+" "+\
            (leftjoins if len(samples)<=cohortsampleslimit else '')+" "+\
            "WHERE "+v_geno+"."+info_cols[0]+" IN "+"(SELECT "+info_cols[0]+" FROM "+v_info+" "+coords2where(v_info,coords,info_cols)+") "+\
            ("AND "+v_geno+".sample"+" IN ("+",".join(["'"+x+"'" for x in samples])+") " if samples!="all" else '')+\
            "GROUP BY "+v_geno+"."+info_cols[0]
    return(query+";")

def coords2where(v_info, coords, tcolumns):
    if coords!='all':
        where = []
        for coord in coords:
            if coord[2]!='':
                #conditions.append(and_(Var.chr == str(coord[0]),Var.start >= str(coord[1]),Var.end <= str(coord[2])))
                where.append("("+v_info+"."+tcolumns[1]+"="+str(coord[0])+\
                       " AND "+v_info+"."+tcolumns[2]+">="+str(coord[1])+ \
                       " AND "+ v_info+"."+tcolumns[3]+"<="+str(coord[2])+")")
            else:
                where.append("(" + v_info + "." + tcolumns[1] + "=" + str(coord[0]) + \
                             " AND " + v_info + "." + tcolumns[2] + "=" + str(coord[1])+")")
        return("WHERE "+" OR ".join(where))
    else:
        return(" ")

def samples2leftjoins(v_geno, lsams, format_cols, tcolumns):
    leftjoins = []
    stables = []
    counter = 1
    for s in lsams:
        ljname = 's'+str(counter)
        scd1 = "LEFT OUTER JOIN "+v_geno+" "+ljname+" ON "+v_geno+"."+tcolumns[0]+"="+ljname+"."+tcolumns[0]+" AND "+ljname+".sample='"+s+"'"
        leftjoins.append(scd1)
        stables.append("||':'||".join([ljname+"."+x for x in ['sample']+format_cols])+" AS "+'"'+s+'"')
        counter+=1
    return(' '.join(leftjoins), ','.join(stables))

def column_dt_builder_all(main_cols, db_sams, Alleles, subquery, mode="nosamples"):
    dt_columns = []
    if mode=="samples":
        all_cols = main_cols + db_sams
    else:
        all_cols = main_cols
    for ac in all_cols:
        if ac in main_cols:
            if ac == "dynamic_ac": dt_columns.append(ColumnDT(func.sum(Alleles.columns.gt).label(ac)))
            elif ac == "dynamic_an": dt_columns.append(ColumnDT((func.count(Alleles.columns.gt)*2).label(ac)))
            elif ac == "dynamic_af": dt_columns.append(ColumnDT((func.sum(Alleles.columns.gt)/(func.count(Alleles.columns.gt)*2)).label(ac)))
            else:
                dt_columns.append(ColumnDT(subquery.__table__.c[ac].label(ac)))
        else:
            samcol = func.max(case([(Alleles.columns.sample == ac, Alleles.columns.gt_raw+":"+Alleles.columns.sample),],))
            dt_columns.append(ColumnDT(samcol.label(ac)))
    return(dt_columns)

def subquery_builder(lsams, db_class, db):
    l = lsams
    lis = [getattr(db_class, x) for x in l]
    lis_af = [func.ifnull(x,0) for x in lis]
    lis_an = [func.ifnull((x+1)/(x+1),0) for x in lis]
    all_the_sum_together = functools.reduce(lambda a,b: a+b, lis_af)
    all_the_an_together = functools.reduce(lambda a,b: a+b, lis_an)*2
    all_the_af_together = all_the_sum_together/(all_the_an_together)
    dquery = db.session.query(db_class.p_id, all_the_sum_together.label('CURRENT_AC'), all_the_af_together.label('CURRENT_AF'),
                              all_the_an_together.label("CURRENT_AN"))
    return(dquery.subquery())


def unique_column_values_dict_builder(sqltablecols, db_sams, Var):
    res_ldict = {}
    for col in sqltablecols:
        if filterable(col.key, db_sams):
            if get_yadcf_type(col.key, db_sams, field_type_dict) == 'multi_select':
                dis_vals = Var.query.with_entities(Var.__table__.c[col.key]).distinct()
                dis_vals = [{"value":x[0], "label":x[0].replace('\\x3b', '-')} for x in dis_vals]
                res_ldict[col.key] = dis_vals
    return(res_ldict)