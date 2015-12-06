#!/usr/bin/env python

import os
import sys
import stat
from distutils.core import setup, Extension

def check_pkg_dependencies():
    pass
    #CHECK for dependencies:
#    try:
#        import django
#    except ImportError, e:
#        sys.stderr.write("CRITICAL: DJANGO 1.1.1 or greater must be installed\n")
#        sys.exit(1)

#    try:
#        import numpy
#        NUMPY_PATH = numpy.__path__[0]
#        #print NUMPY_PATH
#    except ImportError, e:
#        sys.stderr.write("CRITICAL: numpy 1.3 or greater must be installed\n")
#        sys.exit(1)

def check_settings_file():
    # Do not include lib/settings.py in distribution only
    # lib/settings.py.sample; so if there is no lib/settings.py, quit
    # with some information.
#    if os.path.isfile("lib/settings.py"):
#        pass
#    else:
#        sys.stderr.write("CRITICAL: Please copy the lib/settings.py.sample to lib/settings.py, and modify the ASSEMBLY_DIR setting!")
#        sys.exit(1)
    pass
def main():
    if not float(sys.version[:3])>=2.6:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4! python 2.6.1 is recommended!\n")
        sys.exit(1)
    check_pkg_dependencies()
    check_settings_file()
    setup(name="Drseqpipe",
          version="1.00",
          description="Drseq: Drop-seq QC and analysis pipeline",
          author='Shengen Hu',
          author_email='Tarelahu@gmail.com',
          url='https://Tarela@bitbucket.org/tarela/drseq',
          package_dir={'Drseqpipe' : 'lib'},
          packages=['Drseqpipe'],
          package_data={'Drseqpipe': ['Config/Drseq_template.conf',
                                  'Rscript/analysis.r',
                                  'Rscript/individual_qc.r'
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

