#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
from string import *

# --------------------------
# custom package
# --------------------------


### tool function
from Drseqpipe.Utility      import (sp,
                                   sperr,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   CMD
                                   )
# --------------------------
# main 
# --------------------------

def step0_integrate_data(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    wlog("Start Drseq",logfile)
    wlog("Step0: Data integrate",logfile)
    
    ### check data path , format ,
    if not "/" in conf_dict['General']['barcode_file'].strip("/"):
        conf_dict['General']['barcode_file'] = conf_dict['General']['startdir'] + conf_dict['General']['barcode_file']
    if not "/" in conf_dict['General']['reads_file'].strip("/"):
        conf_dict['General']['reads_file'] = conf_dict['General']['startdir'] + conf_dict['General']['reads_file']
    
    if not os.path.isfile(conf_dict['General']['barcode_file']):
        ewlog("barcode file %s not found"%(conf_dict['General']['barcode_file']),logfile)
    if not os.path.isfile(conf_dict['General']['reads_file']):
        ewlog("reads file %s not found"%(conf_dict['General']['reads_file']),logfile)
        
    
    if not conf_dict['General']['barcode_file'].endswith('.fastq') :
        if conf_dict['General']['barcode_file'].endswith('.txt'):
            wlog('barcode file is reformed txt file',logfile)
            conf_dict['General']['format1'] = 'txt'   
        else:
            ewlog("barcode file is not a fastq file: %s"%(conf_dict['General']['barcode_file']),logfile)
    else:
        conf_dict['General']['format1'] = 'fastq'
    if conf_dict['General']['reads_file'].endswith('.fastq') or conf_dict['General']['reads_file'].endswith('.fq'):
        conf_dict['General']['format'] = 'fastq'
        wlog('Detected input file format is fastq',logfile)
    elif conf_dict['General']['reads_file'].endswith('.sam'):
        conf_dict['General']['format'] = 'sam'
        wlog('Detected input file format is sam',logfile)
    else:
        ewlog("reads file is not a fastq or sam file: %s"%(conf_dict['General']['reads_file']),logfile)
    ### check barcode length    
    try : 
        conf_dict['General']['cell_barcode_length'] = int(conf_dict['General']['cell_barcode_length'])
        conf_dict['General']['umi_length'] = int(conf_dict['General']['umi_length'])
    except:
        ewlog('barcode length should be int',logfile)
    ### check gene annotation file
    if conf_dict['General']['gene_annotation'] == "":
        ewlog("gene annotation file cannot be empty",logfile)
    if not "/" in conf_dict['General']['gene_annotation'] : 
        ewlog("absolute path for gene annotation file required",logfile)        
    if not os.path.isfile(conf_dict['General']['gene_annotation'] ):
        ewlog("cannot find gene annotation file : %s"%(conf_dict['General']['gene_annotation'] ),logfile)
        
    ### mapping index
    if conf_dict['General']['format'] == 'fastq':
#        if not conf_dict['Step1_Mapping']['mapindex'].endswith("/"):
#            conf_dict['Step1_Mapping']['mapindex'] += "/"
        if conf_dict['Step1_Mapping']['mapping_software_main'] == "STAR":
            wlog('use STAR as alignment tools',logfile)
#            conf_dict['Step1_Mapping']['mapindex'] +='%s.star'%(conf_dict['General']['genome_version'])
            if not os.path.isdir(conf_dict['Step1_Mapping']['mapindex']):
                ewlog("cannot find STAR index folder : %s"%(conf_dict['Step1_Mapping']['mapindex']),logfile)
        elif conf_dict['Step1_Mapping']['mapping_software_main'] == "bowtie2":
            wlog('use bowtie2 as alignment tools',logfile)
#            indexdir = conf_dict['General']['mapindex'] + '%s.bowtie2/'%(conf_dict['General']['genome_version'])
#            conf_dict['Step1_Mapping']['mapindex'] = indexdir + conf_dict['General']['genome_version']
            indexfile1 = conf_dict['Step1_Mapping']['mapindex']+'.1.bt2'
#           if not os.path.isdir(indexdir):
#               ewlog("cannot find bowtie2 index folder : %s "%(indexdir),logfile)
            if not os.path.isfile(indexfile1):
                ewlog("cannot find bowtie2 index file : %s "%(indexfile1),logfile)
        else:
            ewlog("alignment tools can only be STAR and bowtie2",logfile)


    ### check options
    wlog('option setting: ',logfile)
    try:
        wlog('mapping thread is %s'%(str(int(conf_dict['Step1_Mapping']['mapping_p']))),logfile)
    except:
        ewlog('mapping_p should be int, current value is %s'%(conf_dict['Step1_Mapping']['mapping_p']),logfile)
        
    if not int(conf_dict['Step1_Mapping']['q30filter']) in [0,1]:
        ewlog('q30filter measurement can only be 0/1, current value is %s'%(conf_dict['Step1_Mapping']['q30filter']),logfile)

    if not int(conf_dict['Step3_QC']['bulk_qc']) in [0,1]:
        ewlog('bulk_qc measurement can only be 0/1, current value is %s'%(conf_dict['Step3_QC']['q30filter']),logfile)

    if not int(conf_dict['Step2_ExpMat']['filterttsdistance']) in [0,1]:
        ewlog('filterttsdistance measurement can only be 0/1, current value is %s'%(conf_dict['Step2_ExpMat']['filterttsdistance']),logfile)
    
    if not int(conf_dict['Step2_ExpMat']['ttsdistance']) > 0:
        ewlog('ttsdistance value should greater than 0, current value is %s'%(conf_dict['Step2_ExpMat']['ttsdistance']),logfile)

    if	 int(conf_dict['Step2_ExpMat']['covergncutoff']) > 10000:
        ewlog('covergncutoff value cannot be greater than 10000, current value is %s'%(conf_dict['Step2_ExpMat']['covergncutoff']),logfile)
    
    if not int(conf_dict['Step2_ExpMat']['duplicate_measure']) in [0,1,2,3]:
        ewlog('duplicate_measure value can only be 0~3, current value is %s'%(conf_dict['Step2_ExpMat']['duplicate_measure']),logfile)

    if not int(conf_dict['Step3_QC']['select_cell_measure']) in [1,2]:
        ewlog('select_cell_measure value can only be 1 or 2, current value is %s'%(conf_dict['Step3_QC']['select_cell_measure']),logfile)
    
    if int(conf_dict['Step3_QC']['select_cell_measure']) == 1:
        try:
            int(conf_dict['Step3_QC']['covergncluster'])
        except:
            ewlog('covergncluster value should be integer, current value is %s'%(conf_dict['Step3_QC']['covergncluster']),logfile)
    elif int(conf_dict['Step3_QC']['select_cell_measure']) == 2:
        try:
            int(conf_dict['Step3_QC']['topumicellnumber'])
        except:
            ewlog('topumicellnumber value should be integer, current value is %s'%(conf_dict['Step3_QC']['covergncluster']),logfile)   
    else:
        ewlog('select_cell_measure value can only be 1 or 2, current value is %s'%(conf_dict['Step3_QC']['select_cell_measure']),logfile)

    if not int(conf_dict['Step3_QC']['remove_non_dup_cell']) in [0,1]:
        ewlog('remove_non_dup_cell measurement can only be 0/1, current value is %s'%(conf_dict['Step3_QC']['remove_non_dup_cell']),logfile)
    if float(conf_dict['Step3_QC']['non_dup_cutoff']) <= 0  or float(conf_dict['Step3_QC']['non_dup_cutoff']) >=1 :
        ewlog('non_dup_cutoff measurement should be in 0~1, current value is %s'%(conf_dict['Step3_QC']['non_dup_cutoff']),logfile)
    if float(conf_dict['Step4_Analysis']['highvarz']) <= 0  :
        ewlog('non_dup_cutoff measurement cannot be <= 0, current value is %s'%(conf_dict['Step4_Analysis']['highvarz']),logfile)
    if float(conf_dict['Step4_Analysis']['selectpccumvar']) <= 0 or float(conf_dict['Step4_Analysis']['selectpccumvar']) >=1 :
        ewlog('selectpccumvar measurement should be in 0~1, current value is %s'%(conf_dict['Step4_Analysis']['selectpccumvar']),logfile)
    if not int(conf_dict['Step4_Analysis']['clustering_method']) in [1,2,3,4] :
        ewlog('clustering_method measurement should be chosen from 1,2,3 and 4, current value is %s'%(conf_dict['Step4_Analysis']['clustering_method']),logfile)


    ### check RseQC 
    for qc in ['gb_cover','read_gc','read_nvc','read_qul']:
        if int(conf_dict['Step3_QC']['bulk_qc']) == 1 and sp(conf_dict['Step3_QC'][qc])[0] == "" :
            ewlog('require %s, check whether you install RseQC correctly'%(conf_dict['Step3_QC'][qc]),logfile)
    
    ### check samtools
    if not 'view' in sperr(conf_dict['Step1_Mapping']['samtools_main'])[1]:
        ewlog('require samtools',logfile)
        
    ### check samtools
    if sp(conf_dict['Step1_Mapping']['bedtools_main'])[0] == "":
        ewlog('require bedtools',logfile)

    ### check Rscript
    if not 'Usage' in sperr('Rscript')[1] and not 'version' in sperr('Rscript')[1]:
        ewlog('require Rscript',logfile)
    
    ### check pdflatex
    if sp('pdflatex --help')[0] == "":
        wlog('pdflatex was not installed, Dr.seq is still processing but no summary QC report generated',logfile)

    wlog('Step0 Data integrate DONE',logfile)



    return conf_dict
    
    
    
