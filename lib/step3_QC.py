#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time
# --------------------------
# custom package
# --------------------------

### tool function
from Drseqpipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   createDIR)
# --------------------------
# main 
# --------------------------
def step3_QC(conf_dict,logfile):
    '''
    start RseQC
    mapping stat
    single cell level QC
    '''
    # start
    # create section for 
    
    wlog('Step3: bulk and individual cell QC',logfile)
    ### preparing mapping state dict
    wlog('calculate mapping state',logfile)
    conf_dict['Mapping_stat'] = {}
    conf_dict['Mapping_stat']['umi_gene'] = 0
    conf_dict['Mapping_stat']['cdsN'] = 0
    conf_dict['Mapping_stat']['utr3N'] = 0
    conf_dict['Mapping_stat']['utr5N'] = 0
    conf_dict['Mapping_stat']['intronN'] = 0
    conf_dict['Mapping_stat']['intergenicN'] = 0

    ### calculate mapping state based on QC matrix
    if 1:#conf_dict['General']['dryrun'] == 1:
        inf = open(conf_dict['Step2_ExpMat']['qcmat'])
        for line in inf:
            if line.startswith('cellname'):
                continue
            ll = line.split()
            conf_dict['Mapping_stat']['umi_gene'] += int(ll[2])
            conf_dict['Mapping_stat']['cdsN'] += int(ll[3])
            conf_dict['Mapping_stat']['utr3N'] += int(ll[4])
            conf_dict['Mapping_stat']['utr5N'] += int(ll[5])
            conf_dict['Mapping_stat']['intronN'] += int(ll[6])
            conf_dict['Mapping_stat']['intergenicN'] += int(ll[7])
        inf.close()
        conf_dict['Mapping_stat']['totalreads'] = int(sp('wc -l %s'%(conf_dict['General']['barcode_reform']))[0].split()[0])    
        conf_dict['Mapping_stat']['q30reads'] = int(sp('wc -l %s'%(conf_dict['General']['bed']))[0].split()[0])

   
    ### create  QC dir and conduct QC
    wlog('generate bulk QC measurement',logfile)
    qcdir = conf_dict['General']['outputdirectory'] + 'QC/'
    createDIR(qcdir)
    os.chdir(qcdir)
    conf_dict['QCplots'] = {}
    ### if bulk QC is turned on, conduct bulk QC
    if int(conf_dict['Step3_QC']['bulk_qc']) ==1:
        ## reads quality
        t= time.time()
        cmd = "%s -i %s -o %s"%(conf_dict['Step3_QC']['read_qul'],conf_dict['General']['sam'],conf_dict['General']['outname'])
        rwlog(cmd,logfile,conf_dict['General']['dryrun'])
        ## reads nucleotide composition
        cmd = "%s -i %s -o %s"%(conf_dict['Step3_QC']['read_nvc'],conf_dict['General']['sam'],conf_dict['General']['outname'])
        rwlog(cmd,logfile,conf_dict['General']['dryrun'])
        ## reads GC content
        cmd = "%s -i %s -o %s"%(conf_dict['Step3_QC']['read_gc'],conf_dict['General']['sam'],conf_dict['General']['outname'])
        rwlog(cmd,logfile,conf_dict['General']['dryrun'])
        readsqctime = time.time() -t
        wlog("time for readsqc: %s"%(readsqctime),logfile)
        ## reads genebody coverage
        t= time.time()

        cmd = "%s -i %s -o %s -r %s"%(conf_dict['Step3_QC']['gb_cover'],conf_dict['General']['sam'],conf_dict['General']['outname'],conf_dict['General']['outputdirectory'] + 'annotation/'+conf_dict['General']['outname']+'_gene_anno_fullbed.bed')
        rwlog(cmd,logfile,conf_dict['General']['dryrun'])
        bulkqctime = time.time() -t
        wlog("time for bulkqc: %s"%(bulkqctime),logfile)
        mvcmd1 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.qual.heatmap.pdf',qcdir + conf_dict['General']['outname'] + '_quality_heatmap.pdf')
        mvcmd2 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.NVC_plot.pdf',qcdir + conf_dict['General']['outname'] + '_NVC.pdf')
        mvcmd3 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.GC_plot.pdf',qcdir + conf_dict['General']['outname'] + '_GC.pdf')
        mvcmd4 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.geneBodyCoverage.pdf',qcdir + conf_dict['General']['outname'] + '_GBcover.pdf')
        rwlog(mvcmd1,logfile,conf_dict['General']['dryrun'])
        rwlog(mvcmd2,logfile,conf_dict['General']['dryrun'])
        rwlog(mvcmd3,logfile,conf_dict['General']['dryrun'])
        rwlog(mvcmd4,logfile,conf_dict['General']['dryrun'])

        conf_dict['QCplots']['read_qul'] = qcdir + conf_dict['General']['outname'] + '_quality_heatmap.pdf'
        conf_dict['QCplots']['read_nvc'] = qcdir + conf_dict['General']['outname'] + '_NVC.pdf'
        conf_dict['QCplots']['read_gc'] = qcdir + conf_dict['General']['outname'] + '_GC.pdf'
        conf_dict['QCplots']['gb_cover'] = qcdir + conf_dict['General']['outname'] + '_GBcover.pdf'
    else:
        wlog('bulk QC is turned off, skip bulk QC',logfile)
    
    ### individual cell QC
    wlog('generate individual cell QC measurement',logfile)
    t = time.time()
    conf_dict['QCplots']['duprate'] = qcdir + conf_dict['General']['outname'] + '_duprate.pdf'
    conf_dict['QCplots']['covergn'] = qcdir + conf_dict['General']['outname'] + '_coverGN.pdf'
    conf_dict['QCplots']['intronrate'] = qcdir + conf_dict['General']['outname'] + '_intronrate.pdf'

    if conf_dict['General']['png_for_dot'] == 1:
        conf_dict['QCplots']['umicovergn'] = qcdir + conf_dict['General']['outname'] + '_umi_coverGN.png'
        conf_dict['QCplots']['cumumiduprate'] = qcdir + conf_dict['General']['outname'] + '_cumUMI_duprate.png'
    else:
        conf_dict['QCplots']['umicovergn'] = qcdir + conf_dict['General']['outname'] + '_umi_coverGN.pdf'
        conf_dict['QCplots']['cumumiduprate'] = qcdir + conf_dict['General']['outname'] + '_cumUMI_duprate.pdf'        

    conf_dict['results']['qcmatcc'] = qcdir + conf_dict['General']['outname'] + "_qcmat_clustercell.txt" 
    conf_dict['results']['expmatcc'] = qcdir + conf_dict['General']['outname'] + "_expmat_clustercell.txt" 

    if int(conf_dict['Step3_QC']['select_cell_measure']) ==1:
        use_cutoff = conf_dict['Step3_QC']['covergncluster']
    elif int(conf_dict['Step3_QC']['select_cell_measure']) ==2:
        use_cutoff = conf_dict['Step3_QC']['topumicellnumber']
    else:
        ewlog('select_cell_measure value can only be 1 or 2, current value is %s'%(conf_dict['Step4_Analysis']['select_cell_measure']),logfile)

    cmd = "Rscript %s %s %s %s %s %s %s %s %s %s %s %s"%(conf_dict['rscript']+'individual_qc.r',conf_dict['Step2_ExpMat']['qcmat'],conf_dict['Step2_ExpMat']['expmat'],conf_dict['General']['outname'],conf_dict['Step3_QC']['select_cell_measure'],use_cutoff,conf_dict['Step3_QC']['remove_non_dup_cell'],conf_dict['Step3_QC']['non_dup_cutoff'],conf_dict['Mapping_stat']['umi_gene'],conf_dict['results']['qcmatcc'],conf_dict['results']['expmatcc'],conf_dict['General']['png_for_dot'])
    rwlog(cmd,logfile,conf_dict['General']['dryrun'])
    wlog("Step3 bulk and individual cell QC DONE",logfile)
    individualqctime = time.time() -t
    wlog("time for individualqc: %s"%(individualqctime),logfile)
    return conf_dict








