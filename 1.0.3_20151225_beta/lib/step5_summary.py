#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string

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
                                   createDIR,
                                   textformat,
                                   strlatexformat)
# --------------------------
# main 
# --------------------------
def step5_summary(conf_dict,logfile):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''
    # start
    # create section for 
    
    wlog('Step5: summary',logfile)
    wlog('copy results',logfile)
# Rscript analysis.r expmat outname coverGN highvarZ selectPCcutoff rdnumber maxKnum
    summarydir = conf_dict['General']['outputdirectory'] + 'summary/'
    createDIR(summarydir)
    os.chdir(summarydir)
    
    plot_folder = summarydir + "plots/"
    createDIR(plot_folder)
    os.chdir(plot_folder)
    ### collect results 
    for i in conf_dict['QCplots']:
        if os.path.isfile(conf_dict['QCplots'][i]):
            #realname
            cmd = 'cp %s .'%conf_dict['QCplots'][i]
            rwlog(cmd,logfile)

    result_folder = summarydir + "results/"
    createDIR(result_folder)
    os.chdir(result_folder)
    for i in conf_dict['results']:
        if os.path.isfile(conf_dict['results'][i]):
            cmd = 'cp %s .'%conf_dict['results'][i]
            rwlog(cmd,logfile)

    os.chdir(summarydir)

    wlog('generate qc documents',logfile)
    ### initiate 
    QCdoc = """\documentclass[11pt,a4paper]{article}
\usepackage{tabularx}
\usepackage[english]{babel}
\usepackage{array}
\usepackage{graphicx}
\usepackage{color}
\DeclareGraphicsExtensions{.eps,.png,.pdf,.ps}
\\begin{document}
\\title{QC and analysis reports for Drop-seq data : %s}

\\vspace{-1cm}
\maketitle
\\tableofcontents
\\newpage
\\newpage
\section{Data description}
\\begin{quotation}
Table 1 mainly describe the input file and mapping and analysis parameters.
\end{quotation}
\\begin{table}[h]
\caption{Data description}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }

"""%(strlatexformat(conf_dict['General']['outname']))
    ### table1 prepare parameter
    if int(conf_dict['Step1_Mapping']['q30filter']) == 1:
        q30filter = "True"
    else:
        q30filter = "False"
    if int(conf_dict['Step2_ExpMat']['filterttsdistance']) == 1:
        filtertts = "True"
    else: 
        filtertts = "False"
    if int(conf_dict['Step2_ExpMat']['umidis1']) == 1:
        umidis1 = "True"
    else:
        umidis1 = "False"
    if int(conf_dict['Step3_QC']['remove_non_dup_cell']) == 1:
        rmnodup = "True"
    else:
        rmnodup = "False"
          
    QCdoc += """      
\hline
parameter & value  \\\\
\hline
output name & %s \\\\
\hline
barcode file(file name only) & %s \\\\
\hline
reads file(file name only) & %s \\\\
\hline
reads file format & %s  \\\\
\hline
cell barcode length &  %s \\\\
\hline
UMI length & %s \\\\
\hline
mapping software & %s \\\\
\hline
Q30 filter mapped reads & %s \\\\
\hline
remove reads away TTS & %s \\\\
\hline
"""%(strlatexformat(conf_dict['General']['outname']),
     strlatexformat(conf_dict['General']['barcode_file'].split("/")[-1]),
     strlatexformat(conf_dict['General']['reads_file'].split("/")[-1]),
     conf_dict['General']['format'].upper(),
     str(conf_dict['General']['cell_barcode_length']),
     str(conf_dict['General']['umi_length']),
     conf_dict['Step1_Mapping']['mapping_software_main'],
     q30filter,
     filtertts
     )
    ### table1 part2
    if  filtertts == "True":
        QCdoc += """TTS distance (for remove) & %s bp \\\\
\hline
"""%(str(conf_dict['Step2_ExpMat']['ttsdistance'])) 
    if  int(conf_dict['Step2_ExpMat']['duplicate_measure']) == 1:
        QCdoc += """duplicate rate in each cell & UMI $+$ location \\\\"""
    elif int(conf_dict['Step2_ExpMat']['duplicate_measure']) == 2:
        QCdoc += """duplicate rate in each cell & UMI only \\\\"""
    elif int(conf_dict['Step2_ExpMat']['duplicate_measure']) == 3:
        QCdoc += """duplicate rate in each cell & location only \\\\"""
    else:
        QCdoc += """duplicate rate in each cell & keep all reads \\\\"""
    if int(conf_dict['Step2_ExpMat']['duplicate_measure']) in [1,2]:
        QCdoc += """
\hline
merge UMI ED = 1 & %s \\\\ 
\hline"""%(umidis1)
    if  int(conf_dict['Step3_QC']['select_cell_measure']) == 1:
        QCdoc += """
select STAMPs & %s covered gene \\\\
\hline"""%(str(conf_dict['Step3_QC']['covergncluster']))
    elif int(conf_dict['Step3_QC']['select_cell_measure']) == 2:
        QCdoc += """
select STAMPs & top %s UMI count \\\\
\hline"""%(str(conf_dict['Step3_QC']['topumicellnumber']))
    QCdoc += """
remove low duplicate rate cell & %s \\\\ 
\hline """%(rmnodup)
    if  rmnodup == "True":
        QCdoc += """
low duplicate rate cutoff & %s  \\\\
\hline"""%(str(conf_dict['Step3_QC']['non_dup_cutoff']))
    QCdoc += """
z-score for highly variable gene & %s \\\\ 
\hline 
cumulative variance for selecting PC & %s \\\\
\hline """%(str(conf_dict['Step4_Analysis']['highvarz']),
     str(100*float(conf_dict['Step4_Analysis']['selectpccumvar']))+'\\%')
 
    if  int(conf_dict['Step4_Analysis']['clustering_method']) == 1:
        QCdoc += """
cluster method & k-means (Gap statistics, first stable) \\\\"""
    elif int(conf_dict['Step4_Analysis']['clustering_method']) == 2:
        QCdoc += """
cluster method & k-means (Gap statistics, maxSE) \\\\"""
    elif int(conf_dict['Step4_Analysis']['clustering_method']) == 3:
        QCdoc += """
cluster method & k-means (custom, k=%s) \\\\"""%(conf_dict['Step4_Analysis']['custom_k'])
    else:
        QCdoc += """
cluster method & DBScan (eps=%s) \\\\"""%(conf_dict['Step4_Analysis']['custom_d'])
    QCdoc += """
\hline
\end{tabularx}
\end{table}
"""
    ### bulk QC
    QCdoc += """
\\newpage
\\newpage
\section{Reads level QC}
In the reads level QC step we measured the quality of sequencing reads, including nucleotide quality and composition. In the reads level QC step and Bulk-cell level QC step we randomly sampled down total reads to 5 million and used a published package called ``RseQC" for reference.(Wang, L., Wang, S. and Li, W. (2012) )
\subsection{Reads quality}
\\begin{quotation}
Reads quality is one of the basic reads level quality control methods. We plotted the distribution of a widely used Phred Quality Score at every position of sequence. Phred Quality Score was calculate by a python function $ord(Q) - 33$. Color in the heatmap represented frequency of this quality score observed at this position. Red represented higher frequency while blue was lower frequency.     
\end{quotation}
\\begin{figure}[h]
        \caption{Reads quality} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{Reads nucleotide composition}
\\begin{quotation}
We assessed the nucleotide composition bias of a sample. The proportion of four different nucleotides was calculated at each position of reads. Theoretically four nucleotides had similar proportion at each position of reads. For Drop-seq sample we observed higher A count at the 3'end of reads, because of the 3'end polyA tail generated in sequencing cDNA libaray.
\end{quotation}
\\begin{figure}[h]
        \caption{Reads nucleotide composition} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{Reads GC content}
\\begin{quotation}
Distribution of GC content of each read. 
\end{quotation}
\\begin{figure}[h]
        \caption{Reads GC content} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['QCplots']['read_qul'].split("/")[-1]),
     (conf_dict['QCplots']['read_nvc'].split("/")[-1]),
     (conf_dict['QCplots']['read_gc'].split("/")[-1])
    )

    QCdoc += """
\\newpage
\\newpage
\section{Bulk-cell level QC}
In the bulk-cell level QC step we measured the performance of total Drop-seq reads. In this step we did't separate cell or remove ``empty" cell barcodes, just like treated the sample as bulk RNA-seq sample.
\subsection{Reads alignment summary}
\\begin{quotation}
The following table shows mappability and distribution of total Drop-seq reads. Note that UMI number was calculated by removing duplicate reads (which have identical genomic location, cell barcode and UMI sequences). Mappable reads was after Q30 filtering if Q30 filter function was turned on. \\\\
** the percentage was calculated by dividing total reads number \\\\
*** the percentage was calculated by divding total UMI number
\end{quotation}
\\begin{table}[h]
\caption{Reads alignment summary}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|X| }
    
\hline
genomic region(Category) &  reads number \\\\
\hline
total reads & %s \\\\
\hline
mappble reads &  %s (%s\\%%)* \\\\
\hline
total UMI count & %s (%s\\%%)* \\\\
\hline
CDS exon UMI count & %s (%s\\%%)** \\\\
\hline
3'UTR UMI count & %s (%s\\%%)** \\\\
\hline
5'UTR UMI count & %s (%s\\%%)** \\\\
\hline
intron UMI count & %s (%s\\%%)** \\\\
\hline
intergenic UMI count & %s (%s\\%%)** \\\\
\hline

\end{tabularx}
\end{table}
"""%(textformat(str(conf_dict['Mapping_stat']['totalreads'])),
     textformat(str(conf_dict['Mapping_stat']['q30reads'])),
     str( round(100*conf_dict['Mapping_stat']['q30reads']*1.0/conf_dict['Mapping_stat']['totalreads'], 2)),
     textformat(str(conf_dict['Mapping_stat']['umi_gene'])),
     str( round(100*conf_dict['Mapping_stat']['umi_gene']*1.0/conf_dict['Mapping_stat']['totalreads'], 2)),
     textformat(str(conf_dict['Mapping_stat']['cdsN'])),
     str( round(100*conf_dict['Mapping_stat']['cdsN']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['utr3N'])),
     str( round(100*conf_dict['Mapping_stat']['utr3N']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['utr5N'])),
     str( round(100*conf_dict['Mapping_stat']['utr5N']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['intronN'])),
     str( round(100*conf_dict['Mapping_stat']['intronN']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['intergenicN'])),
     str( round(100*conf_dict['Mapping_stat']['intergenicN']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)))
     ### genebody coverage
    QCdoc += """
\\newpage
\\newpage
\subsection{Gene body coverage}
\\begin{quotation}
Aggregate plot of reads coverage on all genes. Theoretically we observe a unimodal(single bell) distribution, but for Drop-seq sample we observed an enrichment at 3'end because of the CEL-seq like protocol used in sequencing cDNA library. (Klein, A.M., et al. (2015) )
\end{quotation}
\\begin{figure}[h]
        \caption{Gene body coverage} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['QCplots']['gb_cover'].split("/")[-1]))

    QCdoc += """

\\newpage
\\newpage
\section{Individual-cell level QC}
In this step we focused on the quality of individual cell and distinguishing cell barcodes from STAMPs (single-cell transcriptomes attached to microparticles)
\subsection{Reads duplicate rate distribution}
\\begin{quotation}
Drop-seq technology has an innate advantage of detect duplicate reads and amplification bias because of the barcode and UMI information. Here we plotted the distribution of duplicate rate in each cell barcode (though most of cell barcodes don't contain cells, they still have RNA) and observed a bimodal distribution of duplicate rate. We set an option for users to discard cell barcodes with low duplicate rate in following steps. The vertical line represented the cutoff (duplicate rate $>=$ 0.1) of discarding cell barcodes with low duplicate rate.
\end{quotation}
\\begin{figure}[h]
        \caption{Reads dupliate rate distribution} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict['QCplots']['duprate'].split("/")[-1])
    if int(conf_dict['Step3_QC']['select_cell_measure']) == 1:
        QCdoc += """
\\newpage
\\newpage
\subsection{Reads duplicate rate vs. cumulative covered gene number}
\\begin{quotation}
Reads duplicate rate versus cumulative covered gene numbers. Cell barcodes were ranked by the number of covered genes. The duplicate rate (y-axis, left side) was plotted as a function of ranked cell barcode. Red curve represented the number of genes covered by top N cell barcodes (y-axis, right side). N was showed by x-axis. 
\end{quotation}
\\begin{figure}[h]
        \caption{Reads duplicate rate vs. cumulative covered gene number} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{UMI vs. covered gene number}
\\begin{quotation}
Covered gene number was plotted as a function of the number of UMI (i.e. unique read). We observed a clearly different pattern for two groups of cell barcodes with different reads duplicate rate (blue dots versus red and purple dots). Purple dots represented the selected STAMPs for the cell-clustering analysis.Note that we use only STAMPs selected in this step for following analysis. The other cell barcodes are discarded. 
\end{quotation}
\\begin{figure}[h]
        \caption{UMI v.s. covered gene number} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict['QCplots']['cumumiduprate'].split("/")[-1],
     conf_dict['QCplots']['umicovergn'].split("/")[-1])
    else:
        QCdoc += """
\\newpage
\\newpage
\subsection{Reads duplicate rate vs. cumulative covered gene number}
\\begin{quotation}
Reads duplicate rate versus cumulative UMI count. Cell barcodes were ranked by the number of UMI count. The duplicate rate (y-axis, left side) was plotted as a function of ranked cell barcode. Red curve represented cumulative UMI number from top N cell barcodes (y-axis, right side). N was showed by x-axis. 
\end{quotation}
\\begin{figure}[h]
        \caption{Reads duplicate rate vs. cumulative covered gene number} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{UMI vs. covered gene number}
\\begin{quotation}
Covered gene number was plotted as a function of the number of UMI (i.e. unique read). We observed a clearly different pattern for two groups of cell barcodes with different reads duplicate rate (blue dots versus red and purple dots). Purple dots represented the selected STAMPs for the cell-clustering analysis.Note that we used only STAMPs selected in this step for following analysis. The other cell barcodes were discarded. 
\end{quotation}
\\begin{figure}[h]
        \caption{UMI v.s. covered gene number} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict['QCplots']['cumumiduprate'].split("/")[-1],
     conf_dict['QCplots']['umicovergn'].split("/")[-1])
     
    QCdoc += """
\\newpage
\\newpage
\subsection{Covered gene number distribution}
\\begin{quotation}
Histogram of covered gene number of selected STAMPs
\end{quotation}
\\begin{figure}[h]
        \caption{Covered gene number} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{Intron rate distribution}
\\begin{quotation}
Intron rate is a effective method to measure the quality of a RNA-seq sample. We plotted a histogram of intron rate of every STAMP barcodes. Intron rate was defined as $\\frac{intron\\ reads\\ number}{intron + exon\\ reads\\ number}$ 
\end{quotation}
\\begin{figure}[h]
        \caption{Intron rate distribution} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
 
"""%(conf_dict['QCplots']['covergn'].split("/")[-1],
     conf_dict['QCplots']['intronrate'].split("/")[-1])
    
    if int(conf_dict['Step4_Analysis']['clustering_method']) in [3,4]:
        pass
    else:
        if int(conf_dict['Step4_Analysis']['clustering_method']) in [1,2]:
            selectM = 'first stable gap'
        else:
            selectM = 'maxSE'
        QCdoc += """
\\newpage
\\newpage
\section{Cell-clustering level QC}
This step composed by k-means clustering based on t-SNE dimentional reduction result and Gap statistics to determine best k.
\subsection{Gap statistics}
\\begin{quotation}
We conducted a k-means clustering based on t-SNE dimentional reduction output. Gap statistics followed by the ``%s" method was performed to determine the best k in k-means clustering (to determine how many groups the data should have). 
\end{quotation}
\\begin{figure}[h]
        \caption{Gap statistics} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

"""%(selectM,conf_dict['QCplots']['gapstat'].split("/")[-1])
    
    QCdoc += """
\\newpage
\\newpage
\subsection{Clustering plot}
\\begin{quotation}
Scatter plot represented visualization of t-SNE dimensional reduction output of selected STAMP barcodes. STAMP barcodes were colored according to the clustering result and cluster numbers were printed in the center of each cluster. 
\end{quotation}
\\begin{figure}[h]
        \caption{Clustering plot} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
 
"""%(conf_dict['QCplots']['cluster'].split("/")[-1])
    if os.path.isfile(conf_dict['QCplots']['silhouette']):
        QCdoc += """
\\newpage
\\newpage
\subsection{Silhouette of clustering}
\\begin{quotation}
Silhouette method was used to interpretate and validate the consistency within clusters defined in previous steps.  
\end{quotation}
\\begin{figure}[h]
        \caption{Clustering plot} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
 
"""%(conf_dict['QCplots']['silhouette'].split("/")[-1])
    
    QCdoc += """
\\newpage
\\newpage
\subsection{Clustering plot}
\\begin{quotation}
STAMPs was by the total number of UMI based on t-SNE visualization. 
\end{quotation}
\\begin{figure}[h]
        \caption{Clustering plot} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
 
"""%(conf_dict['QCplots']['umicolor'].split("/")[-1])
    
   
    QCdoc += """
\\newpage
\\newpage
\section{Output list}
\\begin{quotation}
All output files were described in the following table
\end{quotation}
\\begin{table}[h]
\caption{output list}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }
    
\hline
description & filename \\\\
\hline
expression matrix for selected STAMPs & %s  \\\\
"""%(strlatexformat(conf_dict['results']['expmatcc'].split("/")[-1]))
    if int(conf_dict['Step4_Analysis']['pctable']) == 1:
        QCdoc += """
\hline
top2 components of PCA dimentional reduction result & %s \\\\         
"""%(strlatexformat(conf_dict['results']['pctable'].split("/")[-1]))
    if int(conf_dict['Step4_Analysis']['cortable']) == 1:
        QCdoc += """
\hline
pairwise correlation matrix & %s \\\\
"""%(strlatexformat(conf_dict['results']['cortable'].split("/")[-1]))
    QCdoc += """
\hline
All features of selected STAMPs & %s \\\\
\hline
summary QC report & %s \\\\
\hline

\end{tabularx}
\end{table} 
\end{document} 
"""%(strlatexformat(conf_dict['results']['features'].split("/")[-1]),strlatexformat(conf_dict['General']['outname'])+"\_summary.pdf")

    os.chdir(plot_folder)

    latexfile = conf_dict['General']['outname'] + '_summary.tex'
    outf = open(latexfile,'w')
    outf.write(QCdoc)
    outf.close()
    cmd = "pdflatex %s"%(latexfile)
    cmd2 = 'cp %s %s'%(conf_dict['General']['outname'] + '_summary.pdf',summarydir)
    if conf_dict['General']['latex'] == 1:
        rwlog(cmd,logfile)
        rwlog(cmd,logfile)
        rwlog(cmd2,logfile)
        for files in os.listdir(plot_folder):
            if os.path.isfile(files) and files[-12:-4] == "_summary":
                if not files[-4:] in ['.tex','.pdf',',png','.txt']:
                    cmd = "rm %s"%(files)
                    rwlog(cmd,logfile)
        wlog('pdflatex was detected in default PATH, generate summary report %s'%('summary/'+conf_dict['General']['outname'] + '_summary.pdf'),logfile)
    else:
        wlog('pdflatex was not detected in default PATH, generate summary report .tex file in summary/plots folder, you can move the whole folder to the environment with pdflatex installed and run cmd: "pdflatex %s"'%(conf_dict['General']['outname'] + '_summary.tex'),logfile)
   
        
    if conf_dict['clean']:
        wlog('--clean pararmeter was turned on, remove internal files with large size',logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_symbol.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_cds.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_3utr.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_5utr.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_TTSdis.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_combined.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_barcode_reform.txt'),logfile)

    wlog('Step5 summary DONE, check %s for final outputs'%(summarydir),logfile)


    return conf_dict









