#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string
import argparse
import subprocess

# -----------------------------------
# custom package
# -----------------------------------
import Drseqpipe

### tool function
from Drseqpipe.Utility      import          (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   readAnnotation,
                                   textformat,
                                   CMD
                                   )
### read and generate config file
from Drseqpipe.parse_config import (gen_conf,
                                   read_conf,
                                   make_conf
                                   )     


                                   
# -------------------
# main step
# -------------------
from Drseqpipe.step0_integrate_data   import step0_integrate_data
from Drseqpipe.step1_generate_matrix         import step1_generate_matrix
from Drseqpipe.step3_QC    import step3_QC
from Drseqpipe.step4_analysis import step4_analysis
from Drseqpipe.step5_summary import step5_summary

# ------------------------------------
# Misc functions
# ------------------------------------

    
#### read options

class ChiLinParser(argparse.ArgumentParser):
    '''
    Ex version for argparse(parameter parser) , add raise error function .
    '''
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()

def parse_args():
    '''
    Read parameter 
    '''
    description = ">.< Drseq pipeline "
    parser = ChiLinParser(description = description)
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    ### generate config file
    template_parser = sub_parsers.add_parser("gen",  help = "generate a template of config file",
                                             description = "Drseq config file generation")
    template_parser.add_argument("-n","--name", dest="config_name",required = True,help="name of your config file : config_name.conf")

    ### run config file
    pipe_parser = sub_parsers.add_parser("run", help = "run pipeline with a config file input",
                                         description = "Run Drseq pipeline with a config file input")
    pipe_parser.add_argument("-c","--config", required = True,
                             help = "specify the config file, -c config_name.conf" )
    pipe_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "The pipeline will over write output result if the output folder is already exist " )
    pipe_parser.add_argument("--clean",dest='Clean' , default=True, action='store_true',
                             help = "remove intermediate result generated during Dr.seq,default is YES" )
    ### simple mode
    simple_parser = sub_parsers.add_parser("simple", help = "run Drseq using simple mode",
                                         description = "(Run Drseq pipeline using simple mode/command line mode) Usage: Drseq.py -a barcode.fastq -b reads.fastq -n outname -s mm10 -g mm10_refgenes.txt --mapindex /yourmapindexfolder/mm10.star")
    simple_parser.add_argument("-b","--barcode", dest = 'barcode',required = True,
                             help = "barcode fastq file before any filtering step, only accept .fastq format" )
    simple_parser.add_argument("-r","--reads",dest='reads',required = True,
                             help = "reads fastq file, accept raw fastq input or aligned sam format, format(fastq,sam) fixed by extension(.fastq or .sam)" )
    simple_parser.add_argument("-n","--name", dest="name",required = True,help="name of you config file and output dir, name only , no extension")
    simple_parser.add_argument("-s","--species",  choices = ("hg38", "mm10"), required = True,
                             help = "species ,choose from hg38 and mm10" )
    simple_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "specify the config file to create output folder , this cmd will rm existing result if set True ~!! " )
    simple_parser.add_argument("--cellbarcodelength",dest='CBL' ,default='12', 
                             help = "specify the length of yoru cell barcode , default is 12(cellbarcode) + 8(umi) = 20 (barcodefastq)" )
    simple_parser.add_argument("--umilength",dest='UMIL',  default='8', 
                             help = "specify the length of yoru UMI , default is 12(cellbarcode) + 8(umi) = 20 (barcodefastq)" )
    simple_parser.add_argument("-g","--gene_annotation",dest='GA', required = False,
                             help = "gene annotation file, the annotation file can be download from UCSC, full annotation text format(see documents for detail), or users can download gene annotation file in hg38 and mm10 version from our homepage" )
    simple_parser.add_argument("--mapindex",dest='mapindex',required = False,
                             help = "mapping index folder, there should be a mm10.star folder under mapindex folder if you use STAR to map to mm10 genome, mm10.bowtie2 folder if use bowtie2. for bowtie2, the index file should named like mm10.1.bt2(see documents for detail)" )
    simple_parser.add_argument("--thread",dest='P' ,default='8', 
                             help = "number of alignment threads to launch, ignored for sam input" )
    simple_parser.add_argument("--clean",dest='Clean' , default=True, action='store_true',
                             help = "remove intermediate result generated during Dr.seq,default is YES" )
    
    
    args = parser.parse_args()
    ## generate config file template 
    if args.sub_command == "gen":
        gen_conf(args.config_name)
        sys.exit(0)

    ## run Dr.seq pipeline with config file input
    if args.sub_command == "run":
        if os.path.isfile(args.config):
            return args
        else:
            print ('ERROR : -c input is not a config file\n')
            print pipe_parser.print_help()
            sys.exit()
    
    ## run Dr.seq pipeline with a simple mode, input parameter in command line
    if args.sub_command == "simple":
        if args.name.endswith('.conf'):
            args.name = args.name[:-5]
        make_conf(args.barcode,args.reads,args.species,args.name,args.fover,args.CBL,args.UMIL,args.GA,args.P,args.mapindex)
        args.config = args.name + '.conf'
        return args

# ------------------------------------
# Main function
# ------------------------------------

def main():

    args = parse_args()
    conf_dict = read_conf(args.config)
    ### read raw path of output dir, the startdir will be used when the input file is not in absolute path
    conf_dict['General']['startdir'] = os.getcwd()+'/'
    
    ### check output name and dir from input parameter
    if conf_dict['General']['outname'] == "":
        print 'your outname cannot be left blank,exit'
        sys.exit(1)
    if conf_dict['General']['outputdirectory'] == "":
        conf_dict['General']['outputdirectory'] = conf_dict['General']['outname']
    if "~" in conf_dict['General']['outputdirectory']:
        print 'require absolute path for outputdirectory'
        sys.exit(1)
    if not conf_dict['General']['outputdirectory'].endswith('/'):
        conf_dict['General']['outputdirectory'] += '/'
    if not "/" in conf_dict['General']['outputdirectory'].rstrip("/"):
        conf_dict['General']['outputdirectory'] = conf_dict['General']['startdir'] + conf_dict['General']['outputdirectory']


    ### creat output dir
    if os.path.isfile(conf_dict['General']['outputdirectory'].rstrip("/")):
        print 'name of your output dir is exist as a file, cannot create a dir,exit'
        sys.exit(1)
    elif os.path.isdir(conf_dict['General']['outputdirectory']):
        if not args.fover:
            print 'name of your output dir is exist as a dir, and overwrite is turned off,exit'
            sys.exit(1)
        else: 
            print 'name of your output dir is exist as a dir, overwrite if turned on, write output result in existing dir'
    else:
		os.system("mkdir %s"%(conf_dict['General']['outputdirectory']))
     
    ### move to output dir
    os.chdir(conf_dict['General']['outputdirectory'])
    ## cp config file to output folder
    cmd = 'cp %s .'%(conf_dict['General']['startdir']+args.config)
    CMD(cmd)
    ### specify the main progress log file
    logfile = conf_dict['General']['outputdirectory']+'progress_log.txt'
    ### remove existing log file. 
    if os.path.isfile(logfile):
        CMD('rm %s'%logfile)
        
    ### Rscript location 
    #CONFIG_TEMPLATE = os.path.join(Drseq_pipe.__path__[0], "Config/Drseq_template.conf")
    conf_dict['rscript'] = os.path.join(Drseqpipe.__path__[0], "Rscript/")#'/mnt/Storage3/CR/Dropseq/drseq/Rscript/'
    conf_dict['clean'] = args.Clean
    ### main step for Dr.seq , see individual script for detail note.
    # preparing step, integrate parameter, prepare for following step
    t = time.time()
    step0_integrate_data(conf_dict,logfile)
    # main data processing step, including mapping, generate expression matrix and QC matrix which is used in next step
    step1_generate_matrix(conf_dict,logfile)
    step1time = time.time() -t
    wlog("running time for matrix generation: %s"%(step1time),logfile)
    # QC step, including bulk RNAseq QC(option), individual cell QC
    t = time.time()
    step3_QC(conf_dict,logfile)
    step3time = time.time()-t
    wlog("running time for QC: %s"%(step3time),logfile)
    # analysis step, including  select cell, filter high variance gene, pca + t-SNE dimentional reduction, k-means + Gap stat clustering
    t = time.time()
    step4_analysis(conf_dict,logfile)
    step4time = time.time() -t
    wlog("running time for clustering: %s"%(step4time),logfile)
    # summary step, integrate all QC figure and expression matrix, generate qC report with latex
    step5_summary(conf_dict,logfile)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)

