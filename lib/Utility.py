#!/usr/bin/env python
"""

Function declare:

def CMD         (cmd)
def sp          (cmd)
def raise_error ()
def pdf_name    (input_name)
def wlog        (message,logfile)
def ewlog       (message,logfile)
def rwlog       (message,logfile,dryrun)

"""
import subprocess
import sys
import os

def CMD(cmd):
    os.system(cmd)

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
    
def raise_error():
    '''
    Raise an error messgae and exit
    '''
    print 'error occurs, check log file~!'
    sys.exit(1)

def pdf_name(input_name):
    '''
    Change filename to pdf file name
    '''
    outputname = "_".join(input_name.split('.')[:-1])+".pdf"
    return outputname
    
def wlog(message,logfile):
    '''
    print a message and write the message to logfile
    '''
    print message
    os.system('echo "%s " >> %s'%(message,logfile))
    
def ewlog(message,logfile):
    '''
    print an error message and write the error message to logfile
    then exit Dr.seq
    error messages start with [ERROR]
    '''
    print "[ERROR] %s "%(message)
    os.system('echo "[ERROR] %s " >> %s'%(message,logfile))
    raise_error()
    
def rwlog(cmd,logfile,DR) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    print "[CMD] %s "%(cmd)
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    if int(DR) == 1:
        pass
    elif int(DR) == 0:
        CMD(cmd)
    else:
        ewlog('dryrun can only be set 0/1',logfile)
   
    
def readAnnotation(annotation):
    '''
    read full annotation file and output as a dictionary 
    file format is fixed to UCSC full annotation format
    '''
    inf = open(annotation)
    outdict = {}
    for line in inf:
        ll = line.split()
        outdict[ll[1]] = ll[12]
    return outdict
    
     
def textformat(inp):
    '''
    transfer 1000000 to  1,000,000 for better visualization in output report
    '''
    o = ''
    comma = 0
    for i in (inp[::-1]):
        comma += 1
        o += i
        if comma%3 == 0:
            o += ','
    return o[::-1].strip(',')

def createDIR(dirname):
    '''
    check dir name and create new dir
    '''
    if not os.path.isdir(dirname):
        os.system('mkdir %s'%(dirname))

def strlatexformat(instr):
    outstr = instr.replace('_','\_')
    return(outstr)
    
def transform_refgene(refgene,ttsdis,outname):

    '''
    transfer full gene annotation downloaded from UCSC to different meterials used in Dr.seq
    including
    1. CDS exon 
    2. 5' utr exon
    3. 3' utr exon
    4. transcript bed
    5. gene symbol bed(merge transcripts belonging to same gene)
    6. TTS +- 400(default) bed
    7. full bed for RseQC input 
    '''
    ret_lst=[]
    inf = open(refgene)
    
    utr5_list = []
    utr3_list = []    
    CDSexon_list = []
    transcript_list = []
    symbol_dict = {}    
    TTS400_list = []
    refbed_list = []
    for line in inf:
        if line.startswith('#'):
            continue
        f = line.strip().split()
        chrom = f[2]
        strand = f[3]
        symbol = f[12]
        txStart = int(f[4])
        txEnd = int(f[5])
        txName = f[1]
        cdsStart = int(f[6])
        cdsEnd = int(f[7])
        exonCount = int(f[8])
        exonStarts = [ int(i) for i in f[9].strip(',').split(',') ]
        exonEnds = [ int(i) for i in f[10].strip(',').split(',') ]
        exonlength = []
        exondisTss = []
        for i in range(len(exonStarts)):
            exonlength.append(exonEnds[i] - exonStarts[i])
            exondisTss.append(exonStarts[i] - txStart)
	    
        if strand == "+" :
            TSS = txStart
            TTS = txEnd
        else:
            TSS = txEnd
            TTS = txStart
        ### collect transcript info
        refbed_list.append([chrom,txStart,txEnd,txName,'0',strand,cdsStart,cdsEnd,'0',exonCount,",".join(map(str,exonlength))+",",",".join(map(str,exondisTss))+","])
        transcript_list.append([chrom,txStart,txEnd,txName,symbol,strand])
        TTS400_list.append([chrom,max(0,TTS-int(ttsdis)),TTS+int(ttsdis)])
        mergeYes = 0
        if not symbol_dict.has_key(symbol):
            symbol_dict[symbol] = [[chrom,txStart,txEnd]]
        else:
            for mergeTX in symbol_dict[symbol]:

                if mergeTX[0] == chrom and mergeTX[1] < txEnd and mergeTX[2] > txStart:
                    mergeTX[1] = min(mergeTX[1],txStart)
                    mergeTX[2] = max(mergeTX[2],txEnd)
                    mergeYes = 1
                    break
                else:
                    mergeYes = 0
            if mergeYes == 0:
                symbol_dict[symbol].append([chrom,txStart,txEnd])
                    
        for st,end in zip(exonStarts,exonEnds):
            if st < cdsStart:
                utr_st = st
                utr_end = min(end,cdsStart)
                utr5_list.append([chrom,utr_st,utr_end])
                if cdsStart < end:
                    CDSexon_list.append([chrom,cdsStart,end])
            if end > cdsEnd:
                utr_st = max(st, cdsEnd)
                utr_end = end
                utr3_list.append([chrom,utr_st,utr_end])
                if st < cdsEnd:
                    CDSexon_list.append([chrom,st,cdsEnd])
            if st >= cdsStart and end <= cdsEnd:
                CDSexon_list.append([chrom,st,end])
    
    ### output tx information
    outf = open(outname + '_gene_anno_fullbed.bed','w')
    for i in refbed_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()        
    
    outf = open(outname+'_gene_anno_transcript.bed','w')
    for i in transcript_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()
    
    outf = open(outname+'_gene_anno_TTSdis.bed','w')
    for i in TTS400_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()

    
    outf = open(outname+'_gene_anno_5utr.bed','w')
    for i in utr5_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()
    
    outf = open(outname+'_gene_anno_3utr.bed','w')
    for i in utr3_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()

    outf = open(outname+'_gene_anno_cds.bed','w')
    for i in CDSexon_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()    
        
    outf = open(outname+'_gene_anno_symbol.bed','w')
    for sym in symbol_dict:
        count=1
        for region in symbol_dict[sym]:
            newll = region + [sym]
            outf.write("\t".join(map(str,newll))+"\n")
            count += 1 
    outf.close()
           
def reform_barcode_fastq(fq,reformtxt,cbL,umiL):
    '''
    transform barcode fastq to another txt format [name,cell_barcode,umi] for following usage
    '''
    lastL = cbL+umiL
    inf = open(fq)
    outf = open(reformtxt,'w')
    count = 0
    for line in inf:
        count += 1
        if count%4 == 1:
            head = line.split()[0][1:]
        if count%4 == 2:
            seq = line.strip()
        if count%4 == 3:
            pass
        if count%4 == 0:
            newll = [head, seq[:cbL], seq[cbL:lastL]]
            outf.write("\t".join(newll)+"\n")
    
    outf.close()
    inf.close()    


def combine_reads(barcodeF,cdsF,utr3F,utr5F,symbolF,ttsdisF,outF,dup_measure):
    '''
    combine annotation information of all reads which is generate in previous step with bedtools
    '''
    barcode_file = open(barcodeF)
    cds_file = open(cdsF)
    utr3_file = open(utr3F)
    utr5_file = open(utr5F)
    symbol_file = open(symbolF)
    TTS400_file = open(ttsdisF)
    outf = open(outF,'w')

    last_read = ["NA"]*8
    next_read_sym = symbol_file.readline().split()
    for line in cds_file:
        current_read_cds = line.split()
        current_read_utr3 = utr3_file.readline().split()
        current_read_utr5 = utr5_file.readline().split()
        current_read_TTS400 = TTS400_file.readline().split()

        read_info = current_read_cds[:6]
        read_name = current_read_cds[3]
        cdsinfo = current_read_cds[6]
        utr3info = current_read_utr3[6]
        utr5info = current_read_utr5[6]
        TTS400info = current_read_TTS400[6]
        newll = []
        if read_name == last_read[3] :
            newll = last_read
        else:
            while(1):
                current_barcode = barcode_file.readline().split()
                if current_barcode[0] == read_name:
                    newll = read_info  + current_barcode[1:] #+ [cdsinfo,utr3info,utr5info]                 
                    break
        if newll == []:
            print 'error in match barcode'
        last_read = newll[:8]
        if len(next_read_sym) > 3 and newll[3] == next_read_sym[3] :
            addsym_list = [next_read_sym[9]]
            while(1):
                next_read_sym = symbol_file.readline().split()
                if len(next_read_sym)>3 and newll[3] == next_read_sym[3]   :
                    if not next_read_sym[9] in addsym_list:
                        addsym_list.append(next_read_sym[9])
                else:
                    break
            addsym = ",".join(addsym_list)
        else:
            addsym = "NA"
                
        newll += [cdsinfo,utr3info,utr5info,TTS400info,addsym] 
        if int(dup_measure) == 1:
            newll[4] = "_".join([newll[7],newll[0],newll[5],newll[1]])
        elif int(dup_measure) == 2:
            newll[4] = newll[7]
        elif int(dup_measure) == 3:
            newll[4] = "_".join([newll[0],newll[5],newll[1]])
        else:
            newll[4] = "NA"                
            
        outf.write("\t".join(newll)+"\n")

    outf.close()

def strdis(str1,str2):
    if len(str1) != len(str2):
        print 'umi have different length in different reads'
        sys.exit(1)
    diff = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            diff += 1
    return diff

def generate_matrix(refgene,inputbed,ttsdis,qcmatfull,qcmat,expmat,coverGNcutoff,umidis1):
    '''
    generate two matrix
    1. expression matrix whose row/column is corresponded to genes and cell_barcodes
    2. QC matrix whose row/column is corresponded to cell_barcodes and measurement(including total reads, #umi, #cds exon reads ,..., )
    '''

    inf = open(refgene)
    allgenes = []
    for line in inf:
        if line.startswith('#'):
            continue
        ll = line.strip().split()
        if not ll[12] in allgenes:
            allgenes.append(ll[12])
    inf.close()

    inf = open(inputbed)

    QCmat = {}
    Expmat = {}

    last_cell = "NA"
    last_umi = "NA"
    for line in inf:
        ll = line.strip().split()
        cell = ll[6]
        umi = ll[4]
    
        if not Expmat.has_key(cell):
            Expmat[cell] = {}
            QCmat[cell] = [0]*9# allreads,umi,cds,3utr,5utr,intron,intergenic,#greatertts400,covered gene number
    
        QCmat[cell][0] += 1
        
        if cell == last_cell and umi == last_umi and umi != "NA" : 
            pass
        elif cell == last_cell and umi != "NA" and int(umidis1) == 1 and strdis(umi,lastumi) == 1 :
            pass
        else:
            if ttsdis ==1 and int(ll[11]) == 0:
                pass
            else:
                if int(ll[8]) > 0 or int(ll[9]) > 0 or int(ll[10]) > 0:
                    if ll[12] != "NA":
                        for target_gene in ll[12].split(","):
                            if not Expmat[cell].has_key(target_gene):
                                Expmat[cell][target_gene] = 0
                            Expmat[cell][target_gene] += 1           
        
            QCmat[cell][1] += 1
            if ll[12] == "NA":
                QCmat[cell][6] += 1
            elif int(ll[8]) > 0:
                QCmat[cell][2] += 1
            elif int(ll[9]) > 0:
                QCmat[cell][3] += 1
            elif int(ll[10]) > 0:
                QCmat[cell][4] += 1
            else:
                QCmat[cell][5] += 1
            if int(ll[11]) > 0:
                QCmat[cell][7] += 1
        
        last_cell =  cell#ll[6]
        last_umi = umi#ll[4]
    
    inf.close()
    outf0 = open(qcmatfull,'w')
    outf1 = open(qcmat,'w')
    outf2 = open(expmat,'w')

    newll = ['cellname','allreads','umi','cds','utr3','utr5','intron','intergenic','awayTTS','coveredGN']
    outf0.write("\t".join(newll)+"\n")
    outf1.write("\t".join(newll)+"\n")

    newll = ['cellname'] + allgenes
    EXPmat = []
    EXPmat.append(newll)
#    outf2.write("\t".join(newll)+"\n")
    for cell in sorted(Expmat.keys()):
        coverN = len(Expmat[cell].keys())
        QCmat[cell][8] = coverN
        newllqc = [cell] + QCmat[cell]
        outf0.write("\t".join(map(str,newllqc))+"\n")
        if coverN >= int(coverGNcutoff):
            outf1.write("\t".join(map(str,newllqc))+"\n")
            explist = []
            for g in allgenes:
                if Expmat[cell].has_key(g):
                    explist.append(Expmat[cell][g])
                else:
                    explist.append(0)
            newll = [cell] + explist
            EXPmat.append(newll)
#            outf2.write("\t".join(map(str,newll))+"\n")
    for i in range(len(EXPmat[0])):
        newll  =[]
        for j in range(len(EXPmat)):
            newll.append(EXPmat[j][i])
        outf2.write("\t".join(map(str,newll))+"\n")
    outf0.close()  
    outf1.close()
    outf2.close()
        
        
   