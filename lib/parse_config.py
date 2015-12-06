#!/usr/bin/env python
# ------------------------------------
"""
Function declare:

def gen_conf (species, name)
def read_conf(configfile)

"""
# -----------------------------------


import os,sys
import ConfigParser

import Drseqpipe

### config template
CONFIG_TEMPLATE = os.path.join(Drseqpipe.__path__[0], "Config/Drseq_template.conf")
#CONFIG_TEMPLATE = '/mnt/Storage3/home/huse/Dropseq/Result/DrseqPipe/Config/Drseq_template.conf'
#CONFIG_TEMPLATE = '/mnt/Storage3/CR/Dropseq/drseq/Config/Drseq_template.conf'
### generate a config
def gen_conf(CONF_name):
    '''
    Generate a config file
    '''
    inf = open(CONFIG_TEMPLATE)
    if not CONF_name.endswith('.conf'):
        CONF_name += '.conf'
    outf = open(CONF_name,'w')
    for line in inf:
        outf.write(line)
    outf.close()
    inf.close()

### read config
def read_conf(conf_file):
    '''
    Read config file and return a dict containing all infomation
    '''
    conf_dict = {}
    cf = ConfigParser.SafeConfigParser()
    cf.read(conf_file)
    #section_list = sorted( cf.sections() , key=lambda x:int(x.lstrip('Step').split("_")[0]))
    for st in cf.sections():
        conf_dict[st]={}
        for item in cf.items(st):
            conf_dict[st][item[0]]=item[1]
    return conf_dict

### generate a config file in simple mode
#        make_conf(args.barcode,args.reads,args.species,args.name,args.fover,args.CBL,args.UMIL,args.RF,args.P)

def make_conf(barcode_file,reads_file,species,outname,fover,cellbarcodeL,umiL,geneanno,P,mapindex):
    inf = open(CONFIG_TEMPLATE)
    name = outname
    if os.path.isfile(name+'.conf') and  not fover :
        print 'config file "%s" using same name exist , choose other name or add -f to overwrite'%(name+".conf")
        sys.exit(1)
    outf = open(name+".conf",'w')
    if geneanno:
        print 1
    else :
        print 0
    for line in inf:
        if line.startswith('barcode_file ='):
            newline = 'barcode_file = ' + barcode_file + '\n'
        elif line.startswith('cell_barcode_length ='):
            newline = 'cell_barcode_length = ' + str(cellbarcodeL) + '\n'
        elif line.startswith('umi_length ='):
            newline = 'umi_length = ' + str(umiL) + '\n'
        elif line.startswith('reads_file ='):
            newline = 'reads_file = ' + reads_file + '\n'
        elif line.startswith('outputdirectory ='):
            newline = 'outputdirectory = ' + name + '\n'
        elif line.startswith('outname ='):
            newline = 'outname = ' + name + '\n'
        elif line.startswith('gene_annotation ='):
            if geneanno:
                newline = 'gene_annotation = ' + geneanno + '\n'
            else:
                newline = line
        elif line.startswith('mapindex ='):
            if mapindex:
                newline = 'mapindex = ' + mapindex + '\n'
            else:
                newline = line
        elif line.startswith('genome_version ='):
            newline = 'genome_version = '+ species + '\n'
        elif line.startswith('mapping_p ='):
            newline = 'mapping_p = ' + str(P) + '\n'
        else:
            newline = line
        outf.write(newline)
    inf.close()
    outf.close()
    return  
