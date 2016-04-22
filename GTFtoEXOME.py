#! /usr/bin/env python

""" 
This script retrieves the exomes of a given ID
Uses the extension '*.GENES' for input and '*.EXOME' for output
"""

import glob, sys, csv
import subprocess

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

# read gene ID from txt ** uses a '*.GENES' extension
names_ruta = 'gst-test.genes'
geneID_names = list()

# retrieve data from TXT file
with open(names_ruta, 'r') as f:
    geneID_names = [line.strip() for line in f]

#create a txt file from each Scaffold

for geneID_name in geneID_names:
    # create db of missing genome files
    comando = "awk '/{0}/ && !/codon/' Dp_geneset_OGS2.gtf".format(geneID_name)
    #comando = 'grep -e {0} Dp_geneset_OGS2.gtf'.format(geneID_name)
    archivo = geneID_name + '.exome'
    log = open(archivo, 'a')  # so that data written to it will be appended
    c = subprocess.Popen(comando, stdout=log, stderr=log, shell=True)
    #subprocess.call(comando.split())
    print archivo
    print "Exome '{0}' successfully created\n".format(geneID_name)


