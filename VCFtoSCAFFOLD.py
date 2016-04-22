#! /usr/bin/env python

""" 
This script runs the fragmentation of a VCF file by Scaffolds, needs a 
scaffold list called 'names-new.genome'
"""

import glob, sys, csv
import subprocess

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

# read gene ID from txt
names_ruta = 'names-new.genome'
scaffold_names = list()

# retrieve data from TXT file
with open(names_ruta, 'r') as f:
    scaffold_names = [line.strip() for line in f]

#create a txt file from each Scaffold

for scaffold_name in scaffold_names:
    # create db of missing genome files
    comando = 'grep -e {0} Dp-var-final.vcf'.format(scaffold_name)
    #comando = 'makeblastdb -in {0}.faa -title db_{1} -dbtype prot -out db_{1}/db_{1} -parse_seqids'.format(dato[0].replace(' ','_'),dato[1])
    archivo = scaffold_name + '.txt'
    log = open(archivo, 'a')  # so that data written to it will be appended
    c = subprocess.Popen(comando, stdout=log, stderr=log, shell=True)
    #subprocess.call(comando.split())
    print archivo
    print "Scaffold '{0}' successfully created\n".format(scaffold_name)


