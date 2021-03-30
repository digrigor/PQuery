#Load Libraries
from flask import render_template, request, redirect, flash, jsonify, session
import sqlite3
import webbrowser
from app import app
from app.handlers import *
from app.dependencies import *
from datatables import ColumnDT, DataTables
from flask_sqlalchemy import SQLAlchemy, get_debug_queries
from sqlalchemy import distinct, func, and_, or_, cast, Numeric, case, Float
from app.models import db, Var, Alleles, FieldTypes, SelectValues, Samples
from werkzeug.datastructures import ImmutableMultiDict

webbrowser.open('http://127.0.0.1:5000/', new=2)

#Get the FieldTypes and Select values tables from the sqlite3 database
ftquery = FieldTypes.query.all()
svquery = SelectValues.query.all()

#setting the order of the output table columns
info_cols_ordered = csvlines2list(columns_order)
real_db_info_cols = [x.index for x in FieldTypes.query.all()]
if len([x for x in real_db_info_cols if x not in info_cols_ordered])>0:
    info_cols_ordered = real_db_info_cols

#Get a list with the output columns which are not dynamically created (the ones directly read by the db
info_cols_no_dynamic = [x for x in info_cols_ordered if x not in ['dynamic_af','dynamic_an','dynamic_ac']]

#LOADING DICTIONARIES CONTAINING INFO ABOUT
#INFO COLUMNS THEIR DATA TYPES AND THE
#UNIQUE VALUES OF THE SELECT COLUMNS
ftquery = FieldTypes.query.all() # This table contains all the INFO column names,
                                 # their SQLite datatype (text,float,etc),
                                 # their filtering method(raw,multiselect etc) and
                                 # whether they are part of the original INFO table
                                 # or they will be processed on the fly later.
svquery = SelectValues.query.all() # This table contains all the INFO column names which
                                   # can be filtered by multi-select and their unique values
                                   # which can be placed in dropdown menus

#Converting the above tables to dictionaries:
select_values = select_dict_builder(svquery)
field_types = field_type_dict_builder(ftquery)

#Format the alleles table and columns
alleles_cols = Alleles.__table__.columns.keys()
format_cols = alleles_cols[2:]

#Get all the sample names
db_samples = [x.value for x in Samples.query.all()]

#Home-Page
@app.route('/', methods=["GET", "POST"])
@app.route('/index', methods=["GET", "POST"])
def index():
    #get all samples
    colnames = db_samples

    #handle form input from the user
    if request.method == 'POST':
        req = request.form
        reqf = request.files

        #Get the user input
        insams,ingenes,infilt = handle_filter_req(req, reqf)

        #Check for errors (if user hasn't provided input, or if user has selected more than one proposed
        #ways to select samples and genes.
        errlist = handle_input_error_messages(return_status(insams, True),return_status(ingenes, True))

        #If any errors: prints out warning message
        if len(errlist)!=0:
            flash(errlist,"warning")
            return redirect(request.url)
        else:
            #No errors detected:
            coords = handle_genes_input(ingenes) #input genes to coordinates conversion
            all_samples = handle_samples_input(insams) #input samples to sample list
            errlist = handle_input_error_messages(return_status(all_samples, False), return_status(coords, False)) #check input for errors
            #print(coords)
            #If any errors (maybe the user gene was mispelled so no coordinates could be identified):
            if len(errlist) != 0:
                flash(errlist, "warning")
                return redirect(request.url)
            else:
                #Successfull submission - Store input in 'session' objects
                #so they can be passed to other routes.
                flash(["Successfull Submission!"],"success")
                session["dt_samples"] = all_samples
                session["dt_coords"] = coords
                session["dt_filt"] = infilt
                session['fpane_clean_dict'] = ''
                return redirect("/var_table")

    return render_template('index.html', title='Home', samplenames_list=colnames, slider_vals=gnomad_filt_conv_dict,
                           cohorts_dict=coh2show_dict, genregs_list=genes_drop_options)



#Output page - Queried Vatiants output:
@app.route("/var_table", methods=["GET", "POST"])
def var_table():
    #Get user input from the filter pane (top of the page) IF THE USER decides to filter:
    if request.method == 'POST':
        #print("POST")
        req = request.form

        #Store user filter selection in a dictionary
        fpane_dict = req.to_dict(flat=False) #request object to dict
        fpane_dict_clean = {x:fpane_dict[x] for x in fpane_dict if fpane_dict[x]!=['']}

        #print(fpane_dict_clean)
        #Save user's filtering preferences in a session object
        session['fpane_clean_dict'] = fpane_dict_clean
        return redirect(request.url)
    """List Variants."""
    try:
        # if queried samples < cohort_cutoff, PQuery will return all info columns and the individual genotypes per sample.
        if len(session["dt_samples"])<=cohort_cutoff:
            tab_cols = info_cols_ordered + session["dt_samples"]
        else: #If queried samples > cohort_cutoff, then PQuery will run in cohort mode without showing individual genotypes.
            tab_cols = info_cols_ordered
        #Get the final dictionary which holds the info for the columns that will be returned:
        divs = filter_ldict_builder(info_cols_ordered, session["dt_samples"], field_types, select_values)
    except KeyError:
        #If any error happens with the session object, try this alternative:
        tab_cols = info_cols_ordered + db_samples
        divs = filter_ldict_builder(info_cols_ordered, db_samples, field_types, select_values)
    #print(list(range(int(db_cols.index("FORMAT") + 1), len(tab_cols))))
    return render_template('var_table.html', title='SGUL Genetics Exome Data', db_cols=tab_cols, divs=divs, sfos=sample_filters_options)

#IGV page:
@app.route("/igv", methods=["GET","POST"])
def igv_page():
    try:
        #User can select which samples to show IGV from the list of
        #samples initially selected.
        sams = session["dt_samples"]
    except KeyError:
        #If any errors with the session object, show IGV for ALL samples.
        sams = db_samples
    if request.method == 'POST':
        #Handle user selection of samples:
        #print("POST")
        req = request.form
        req_dict = req.to_dict(flat=False)
        try:
            #Save the diotionary with the sample names and bam locations for IGV input.
            session["igv_dict"] = [{'name':x, 'url':igv_samples_bams_bais_dict[x][0], 'indexURL':igv_samples_bams_bais_dict[x][1], 'format':'bam'} for x in req_dict["igv_samples"]]
        except KeyError:
            session["igv_dict"] = []
        #print(session["igv_dict"])
        return redirect(request.url)
    try:
        igv_dict =  session["igv_dict"]
    except KeyError:
        igv_dict =  []
    return render_template('igv.html', igv_sams=sams, igv_dict=igv_dict)

# Function that creates the final query, executes it and then
# returns the final object which is being returned by var_Table page.
@app.route('/data')
def data():
    """Return server side data."""
    # Getting html datatables parameters
    params = request.args.to_dict()

    # Getting samples selected
    lsams = session["dt_samples"]

    #Get all the whole table queried as a subquery (So can be used in SELECT * FROM this_subquery) for example
    all_info_query = db.session.query(Var).subquery()

    # Build the initial query:
    ## 1. Query all variants first -> db.session.query(Var.variant_id)
    ## 2. Then pass to the builder the user-selected coords -> session["dt_coords"]
    ## 3. Build the query
    filter_info_query_in = initial_query_builder(db.session.query(Var.variant_id),session["dt_coords"],Var)

    #Do the same but as a subquery (o can be used in SELECT * FROM this_subquery)
    filter_info_subquery_in = initial_query_builder(db.session.query(Var.variant_id), session["dt_coords"], Var).subquery()

    # Query the table with the genotypes (Alleles model)
    ## select_from(Alleles): which is the basic table
    ##
    ## filter(Alleles.variant_id.in_(filter_info_subquery_in): Filter the super
    ## long Alleles table to only include variants which are also in the above-queried info table.
    ##
    ## Alleles.sample.in_(lsams): Filter the alleles table to only include rows with samples that the user has selected.
    ##
    ## func.sum(Alleles.gt)>0: Only include samples which have at least one variant in the given coordinates
    geno_query = db.session.query().select_from(Alleles).filter(Alleles.variant_id.in_(filter_info_subquery_in), Alleles.sample.in_(lsams)).group_by(Alleles.variant_id).having(func.sum(Alleles.gt)>0)

    nondynamics_state=False

    #If user has performed filter on the var_Table:
    if session['fpane_clean_dict']!='':
        fpd = session['fpane_clean_dict']
        nondynamics_state = (len([y for y in ['dynamic_ac', 'dynamic_af', 'dynamic_an'] for x in list(fpd.keys()) if y not in x]) > 0)
        # True if AC,AN.AF should be calculated on the fly for the selected samples, FALSE is not. It's always TRUE

        # Filter the initial query according to the user-defined filters:
        filter_info_query_in_filtered = fpane_query_builder_var(filter_info_query_in,session['fpane_clean_dict'], lsams, Var, field_types)
        filter_info_subquery_in = filter_info_query_in_filtered.subquery()

        #Re-do the geno_query as above but now filter to have same variants with the filter_info_query_in filtered
        geno_query = db.session.query().select_from(Alleles).filter(Alleles.variant_id.in_(filter_info_subquery_in), Alleles.sample.in_(lsams)).group_by(Alleles.variant_id).having(func.sum(Alleles.gt)>0)
        #Filter sample columns (hom,het etc..) and dynamic_af,dynamic_an,dynamic_ac columns:
        geno_query_filtered = fpane_query_builder_having(geno_query, session['fpane_clean_dict'], lsams, Alleles, field_types)
        geno_query = geno_query_filtered


    #Join varquery and geno query
    query = geno_query.join(all_info_query, Alleles.variant_id == all_info_query.columns["variant_id"])

    # Here we control which columns of the above-built queries will be returned in the final table.
    #'# if queried samples < cohort_cutoff, PQuery will return all info columns and the individual genotypes per sample.
    if len(lsams) <= cohort_cutoff:
        cols = column_dt_builder(info_cols_ordered, lsams, Alleles, all_info_query, "samples")
    ## If queried samples > cohort_cutoff, then PQuery will run in cohort mode without showing individual genotypes.
    elif len(lsams) > cohort_cutoff:
        cols = column_dt_builder(info_cols_ordered, lsams, Alleles, all_info_query)

    # Here we specify the number of columns returned and building the final object which will be returned
    # on the var_table page:
    if session["dt_coords"]!="all" or nondynamics_state==False:
        # HERE IS A BUG:
        # We just count the total rows before the filtering as otherwise it
        # would be super slow:
        rowcount = filter_info_query_in.count()
        rowTable = DataTables(params, query, cols, rowcount)
        return jsonify(rowTable.output_result())

    # If the user selects to return whole exome we rebuild the queries to deal with this huge task:
    else:
        print("2exome")
        geno_query = db.session.query(Alleles).filter(Alleles.sample.in_(lsams)).subquery()
        infoquery = db.session.query().select_from(Var).group_by(geno_query.columns.variant_id).having(func.sum(geno_query.columns.gt)>0)
        query = infoquery.join(geno_query, geno_query.columns.variant_id == Var.variant_id)
        if len(lsams) <= cohort_cutoff:
            cols = column_dt_builder_all(info_cols_ordered, lsams, geno_query, Var, mode="samples")
        elif len(lsams) > cohort_cutoff:
            cols = column_dt_builder_all(info_cols_ordered, lsams, geno_query, Var)
        rowcount_query =  db.session.query(Alleles.variant_id).filter(Alleles.sample.in_(lsams)).group_by(Alleles.variant_id).having(func.sum(Alleles.gt)>0)
        rowcount = rowcount_query.count()
        rowTable = DataTables(params, query, cols, rowcount)
        return jsonify(rowTable.output_result())


#PREVIOUS CHUNKS OF CODE, MIGHT BE USEFUL
#     try:
#         if session['fpane_clean_dict']!='':
#             print(session['fpane_clean_dict'])
#             fpd = session['fpane_clean_dict']
#             dynamics_state = len([y for y in ['dynamic_ac', 'dynamic_af', 'dynamic_an'] for x in list(fpd.keys()) if y in x])>0
#             nondynamics_state = len([y for y in ['dynamic_ac', 'dynamic_af', 'dynamic_an'] for x in list(fpd.keys()) if y not in x])>0
#
# def final_query_builder(db, Var, Alleles, coords, filt, insams, fpd, params, dynamics_state, nondynamics_state):
#
# general_genoquery = db.session.query(A)
# genoquery_pre = db.session.query(Alleles.variant_id, func.sum(Alleles.gt).label('dynamic_AC'),
#                                  (func.count(Alleles.gt) * 2).label('dynamic_AN'),
#                                  (func.sum(Alleles.gt) / (func.count(Alleles.gt) * 2)).label('dynamic_AF')).filter(
#     Alleles.sample.in_(lsams)).group_by(
#     Alleles.variant_id)
#     if dynamics_state == False & nondynamics_state == False:
#         dt_columns = dt_columns_builder()
#         varquery = db.session.query().select_from(Var)
#         genoquery_final = genoquery_pre.limit(params.get('length')).offset(params.get('start')).subquery()
#         varquery_final = varquery
#         varquery_final2 = initial_query_builder(varquery_final, session["dt_coords"], session["dt_filt"], Var)
#         varquery_final3 = fpane_query_builder_var(varquery_final2, fpd, lsams, Var, field_types)
#         final_query = varquery_final3.join(genoquery_final, genoquery_final.columns.variant_id == Var.variant_id)
#         row_count = varquery_final3.add_columns(*[c.sqla_expr for c in [ColumnDT(Var.variant_id)]]).count()
#     elif dynamics_state == False & nondynamics_state == True:
#         dt_columns = dt_columns_builder()
#         genoquery_final = genoquery_pre.subquery()
#         varquery_sq1 = db.session.query(Var)
#         varquery_sq1 = initial_query_builder(varquery_sq1, session["dt_coords"], session["dt_filt"], Var)
#         varquery_sq1 = fpane_query_builder_var(varquery_sq1, fpd, lsams, Var, field_types)
#         varquery_sq1 = varquery_sq1.limit(params.get('length')).offset(params.get('start')).subquery()
#         varquery_sq2 = db.session.query().select_from(varquery_sq1).add_columns(*[c.sqla_expr for c in dt_columns])
#
# main_geno = db.session.query(Alleles)
#
# main_geno = db.session.query(mc.variant_id, func.sum(mc.gt).label('AC'),
#                  (func.count(mc.gt)*2).label('AN'),
#                  (func.sum(mc.gt)/(func.count(mc.gt)*2)).label('AF'),
#                  func.group_concat(case([(mc.sample == "S2201",  mc.gt+":"+mc.dp+":"+mc.gq+":"+mc.pgt+":"+mc.pid+":"+mc.sample)])).label('S2201'),
#                  func.group_concat(case([(mc.sample == "S2202",  mc.gt+":"+mc.dp+":"+mc.gq+":"+mc.pgt+":"+mc.pid+":"+mc.sample)])).label('S2202'),
#                  func.group_concat(case([(mc.sample == "S2203",  mc.gt+":"+mc.dp+":"+mc.gq+":"+mc.pgt+":"+mc.pid+":"+mc.sample)])).label('S2203')).group_by(main_sub.columns.variant_id)
#
#         final_query = varquery_sq2.join
#
#
#
#             if dynamics_state == True:
#                 if nondynamics_state == True:
#                     print('case1')
#                     genoquery_fpane_filt_pre = fpane_query_builder_geno(genoquery_pre, fpd, lsams, Alleles, field_types)
#                     row_count = genoquery_fpane_filt_pre.count()
#                     genoquery_final = genoquery_fpane_filt_pre.subquery()
#                     varquery_final = varquery
#                     final_query = varquery_final.join(genoquery_final, genoquery_final.columns.variant_id == Var.variant_id)
#                     final_query_infilt = initial_query_builder(final_query, session["dt_coords"], session["dt_filt"], Var)
#                     final_query_infilt_fpanefilt = fpane_query_builder_var(final_query_infilt, fpd, lsams, Var, field_types).limit(params.get('length')).offset(params.get('start'))
#                 elif nondynamics_state == False:
#                     print('case2')
#                     genoquery_fpane_filt_pre = fpane_query_builder_geno(genoquery_pre, fpd, lsams, Alleles, field_types)
#                     row_count = genoquery_fpane_filt_pre.count()
#                     genoquery_final = genoquery_fpane_filt_pre.limit(params.get('length')).offset(params.get('start')).subquery()
#                     varquery_final = varquery
#                     final_query = varquery_final.join(genoquery_final,genoquery_final.columns.variant_id == Var.variant_id)
#                     final_query_infilt = initial_query_builder(final_query, session["dt_coords"], session["dt_filt"],Var)
#                     final_query_infilt_fpanefilt = final_query_infilt
#                 session['fpane_clean_dict'] = ''
#             else:
#                 print('case3')
#                 varquery_final = varquery
#                 varquery_final2 = initial_query_builder(varquery_final, session["dt_coords"], session["dt_filt"], Var)
#                 varquery_final3 = fpane_query_builder_var(varquery_final2, fpd, lsams, Var, field_types)
#                 varquery_sq = varquery_final3.add_columns(*[c.sqla_expr for c in [ColumnDT(Var.variant_id)]]).subquery()
#                 genoquery_fpane_filt_pre = genoquery_pre
#                 genoquery_final = genoquery_fpane_filt_pre.join(varquery_sq,  = ).limit(params.get('length')).offset(params.get('start')).subquery()
#                 final_query = varquery_final3.join(genoquery_final, genoquery_final.columns.variant_id == Var.variant_id)
#                 final_query_infilt_fpanefilt = fpane_query_builder_var(final_query, fpd, lsams, Var, field_types)
#                 row_count = varquery_final3.add_columns(*[c.sqla_expr for c in [ColumnDT(Var.variant_id)]]).count()
#                 session['fpane_clean_dict'] = ''
#         else:
#             print('case4')
#             #print("fpane_query_else")
#             #print(len(session['fpane_clean_dict']))
#             #print(session['fpane_clean_dict'])
#             dynamics_state = False
#             genoquery_final = genoquery_pre.limit(params.get('length')).offset(params.get('start')).subquery()
#             varquery_final = varquery
#             final_query = varquery_final.join(genoquery_final, genoquery_final.columns.variant_id == Var.variant_id)
#             final_query_infilt = initial_query_builder(final_query, session["dt_coords"], session["dt_filt"], Var)
#             final_query_infilt_fpanefilt = final_query_infilt
#             row_count = varquery_final.add_columns(*[c.sqla_expr for c in [ColumnDT(Var.variant_id)]]).count()
#             session['fpane_clean_dict'] = ''
#     except KeyError:
#         print('case5')
#         genoquery_final = genoquery_pre.limit(params.get('length')).offset(params.get('start')).subquery()
#         varquery_final = varquery
#         final_query = varquery_final.join(genoquery_final, genoquery_final.columns.variant_id == Var.variant_id)
#         final_query_infilt = initial_query_builder(final_query, session["dt_coords"], session["dt_filt"], Var)
#         final_query_infilt_fpanefilt = final_query_infilt
#         row_count = final_query_infilt_fpanefilt.add_columns(*[c.sqla_expr for c in [ColumnDT(Var.variant_id)]]).count()
# dt_columns = []
# for x in info_cols_ordered+lsams:
#     if x in lsams: dt_columns.append(ColumnDT(Var.__table__.c[x]))
#     else:
#         if field_types[x]['db_status'] == 'in_db': dt_columns.append(ColumnDT(Var.__table__.c[x]))
#         elif field_types[x]['db_status'] == 'out_db':
#             if '_ac' in x: dt_columns.append(ColumnDT(func.sum(Alleles.gt)))
#             elif '_an' in x: dt_columns.append(ColumnDT(func.count(Alleles.gt)*2))
#             elif '_af' in x: dt_columns.append(ColumnDT(func.sum(Alleles.gt)/(func.count(Alleles.gt)*2)))
#     #print(fquery)
#     rowTable = DataTables(params, final_query_infilt_fpanefilt, dt_columns, row_count=row_count)
#     #print(rowTable)
#     return jsonify(rowTable.output_result())