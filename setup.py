#!/usr/bin/env python
"""Description
Setup script for Dr.seq  -- QC and analysis pipeline for Drop-seq data
Copyright (c) 2015 Shengen Hu <tarelahu@gmail.com>
This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""
import os
import sys
import subprocess
from distutils.core import setup, Extension


if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: Dr.seq requires Python 2.7"
	sys.exit()

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
	
def compile_bedtools():
    curdir = os.getcwd()
    os.chdir('refpackage/bedtools')
    sp('make > makeoutput.out')
    sp('chmod 755 *')
    os.chdir(curdir)
    
def check_bedtools():
    checkhandle = sp('which bedtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
        
def main():
    if not float(sys.version[:3])>=2.6:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4! python 2.6.1 is recommended!\n")
        sys.exit(1)
    has_bedtools = check_bedtools()
    if has_bedtools == 0:
        print 'bedtools is not detected under default PATH, now installing bedtools'
        print 'The installation of bedtools will take serval minutes'
        compile_bedtools()
        print 'bedtools compiled, now installing Dr.seq'
        setup(name="Drseqpipe",
              version="1.0",
              description="Drseq: Drop-seq QC and analysis pipeline",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://Tarela@bitbucket.org/tarela/drseq',
              package_dir={'Drseqpipe' : 'lib'},
              packages=['Drseqpipe'],
              package_data={'Drseqpipe': ['Config/Drseq_template.conf',
                                      'Rscript/analysis.r',
                                      'Rscript/individual_qc.r',
                                      'Rscript/readsbulkQC.r'
                                         ]},
              scripts=['bin/Drseq.py','refpackage/bedtools/bin/bedtools'],
                        
              classifiers=[
            'Development Status :: version1.0 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=[],
          )
    
    else:
        print 'bedtools detected, now installing Dr.seq' 
        setup(name="Drseqpipe",
              version="1.0",
              description="Drseq: Drop-seq QC and analysis pipeline",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://Tarela@bitbucket.org/tarela/drseq',
              package_dir={'Drseqpipe' : 'lib'},
              packages=['Drseqpipe'],
              package_data={'Drseqpipe': ['Config/Drseq_template.conf',
                                      'Rscript/analysis.r',
                                      'Rscript/individual_qc.r',
                                      'Rscript/readsbulkQC.r'
                                         ]},
              scripts=['bin/Drseq.py'],
                        
              classifiers=[
            'Development Status :: version1.0 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=[],
          )


if __name__ == '__main__':
    main()

