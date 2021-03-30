import pandas as pd
import os



#CONFIGURATION
maindir="/home/dgrigoriadis/postergaard/NGS/Genomes_from_GEL/hg38/"
#annotfile=bcfdir+"60k_GRCH38_germline_mergedgVCF_all_chr.vep.annotated.vcf.gz"
upstream=10000
downstream=10000
cohort_cutoff=40

#DEPENDENCIES
repattern="(^chr(\d{1,2}|[XY])|(^\d{1,2}|[XY]))($|:\d+)($|-\d+$)"

#SETTINGS
unique_delimiter = 'xdio3'

#SQLITE3 DB AND DB RELATED FILES

columns_order = "app/resources/poster_072019.db.columns-order.csv"



#Cohorts_dict
coh2show_dict = {'lymphoedema':'coh_lymphoedema_coh',
                'lipoedema':'coh_lipoedema_coh'}

coh2file_dict = {'coh_lymphoedema_coh':'app/resources/poster_072019_Lymphoedema_samples.csv',
                'coh_lipoedema_coh':'app/resources/poster_072019_Lipoedema_samples.csv'}





vardb = "app/resources/varDB"

gnomad_filt_conv_dict = {"No Filter" : 1, 0.01 : 0.01, 0.02 : 0.02, 0.03 : 0.03, 0.04 : 0.04, 0.05 : 0.05, 0.075 : 0.075, 0.1 : 0.1, 0.3 : 0.3, 0.5 : 0.5}


gene2loc = "app/resources/ensembl_hg19_genes_to_locations.txt"
g2f = pd.read_csv(gene2loc, sep="\t")

genes_drop_options = {y:y for y in ['chr'+str(x) for x in list(range(1,23))+['X','Y']]}
gnomad_filt_numbers_list = [0,1,2,3,4,5,7,10]
gnomad_genome_colname =  "gnomad_genome_all"
gnomad_exome_colname = "gnomad_genome_all"
sgul_af_colname = "AF"
operators = [['ge ', '>='],
             ['le ', '<='],
             ['lt ', '<'],
             ['gt ', '>'],
             ['ne ', '!='],
             ['eq ', '='],
             ['contains '],
             ['datestartswith ']]

initial_af_filters_dict = {'filt-gnex-drop':gnomad_exome_colname,
                           'filt-gnge-drop':gnomad_genome_colname,
                           'filt-sgex-drop':sgul_af_colname}

sample_filters_options = {"refhom" : "0/0",
                          "althom" : "1/1",
                          "het" : "0/1",
                          "car" : "1/1 or 0/1"}

multi_select_columns = ["CHROM","FILTER","Func_refGene","Func_ensGene",
                        "ExonicFunc_refGene","ExonicFunc_ensGene","LRT_pred",
                        "MetaLR_pred","MetaSVM_pred","MutationAssessor_pred","POSITIVE_TRAIN_SITE","PROVEAN_pred",
                        "Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","SIFT_pred","culprit",
                        "fathmm_MKL_coding_pred","gnomAD_exome_Filter",
                        "gnomAD_genome_Filter","Eigen_coding_or_noncoding"]
range_columns = ["POS","QUAL","AF","1000g","AN","BaseQRankSum","CADD","DANN","Eigen-PC-raw","Eigen-raw"
                 "EXAC","ExcessHet","FATHMM_converted_rankscore","FS","GERP_RS","GERP_RS_r_mammalian"
                 "GME","GenoCanyon_score","GenoCanyon_score_rankscore","HRC","InbreedingCoeff","Kaviar",
                 "LRT_score","LRT_converted_rankscore"]
text_exact = ['p_id', 'ID', 'avsnp150']
text_contains = ['Gene_refGene','Gene_ensGene','GeneDetail_refGene','GeneDetail_ensGene']

bams_walk = os.walk("app/static/bams/")
bams_walk2 = list([d[2] for d in bams_walk])[0]
igv_samples_bams_bais_dict = {x.split("_sorted_")[0]:['static/bams/'+x,'static/bams/'+x.split(".bam")[0]+".bai"] for x in bams_walk2 if x.endswith(".bam")}

#bams = glob.glob("app/static/bams/*.bam")
#igv_samples_bams_bais_dict = {x.split("/")[-1].split("bamout")[0].split("_sorted_unique")[0]:[x.split("app/")[-1],x.split(".bam")[0].split("app/")[-1]+".bai"] for x in bams}


PAGE_SIZE=1000

