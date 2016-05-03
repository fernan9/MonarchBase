#! /usr/bin/env python

"""
---DESCRIPTION
This script calls SNPS from simplified versions of VCF according to a gene ID
-outputs the codon position of the SNP first
-eventually will assign aminoacid change

---OUTPUT format
output data should include
'01genid' gene id
'02iniex' initial exome
'03endex' end of exome
'04strnd' strand
'05frame' frame 0,1,2
'06chrom' chromosome
'07posit' position
'08snpid' id
'09refer' reference base
'10alter' alternative base
'11qualy' quality
'12filtr' filter
'13infor' information
'14formt' format
'15gntyp' genotype
-- CODON	flanquing sequences**
-- ANCAA	ancient aminoacid**
-- VARAA	variable aminoacid**

** (position and strand dependent) UNLESS we could call from FASTA
   but local position of snp (genewise) is needed in that case
   special cases appear when the snp is in any border of exon

---EXTENSIONS
External info:
-VCF
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
-GTF
http://www.ensembl.org/info/website/upload/gff.html

Internal info:
-EXOME
files named with the coordinate data
-SCAF
files termed with the names of Scaffolds
-VAR
files containing the variable data for each gene
"""

import glob, sys, csv
import subprocess
import pandas as pd
from StringIO import StringIO
from subprocess import check_output
import numpy as np
import time

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

# 1. First we read all the files in local directory with the extension .EXOME
exomes = glob.glob('*.exome')
paths = glob.glob('*/')
#out filename
filename = time.strftime("%d.%m.%Y_%H.%M.%S")
#out df
col = ['01genid', '02iniex','03endex','04strnd','05frame','06chrom','07posit','08snpid','09refer','10alter','11qualy','12filtr','13infor','14formt','15gntyp']
salida_df = pd.DataFrame(np.nan, index=[0], columns=col)

# 2. Second for each file (describing one gene by exons)
for exome in exomes:
    geneID = exome.split('.')[0]
    print geneID
#	A. load the file as a DataFrame with pandas
    exome_df = pd.read_table(exome, header = None)
#	B. Read each row to apply
    for exon in exome_df.itertuples():
#		a. drive a masking of the related SCAF file
# 			$ awk '{ if ($2 >= limites(433041) && $2 <= limites(433055)) print }' DPSCF300328.scaf
        scaf = exon[1] + '.scaf'
        limit = (exon[4],exon[5])
        comando = "awk '{{ if ($2 >= {0} && $2 <= {1}) print }}' {2}".format(min(limit),max(limit), scaf)
#       next evaluation catches exceptions to avoid crash
        try:
            salida = subprocess.check_output(comando, stderr=subprocess.STDOUT, shell = True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        if salida == "": #this line catches an empty return to avoid crash
            continue
#		b. create a DataFrame of it      
        temp_df = pd.read_table(StringIO(salida), header = None)
        temp_df.columns = ['06chrom','07posit','08snpid','09refer','10alter','11qualy','12filtr','13infor','14formt','15gntyp']
        tamano = len(temp_df.index)
#		c. repeat data from each Exon Row and merge it to the masked Scaffold
        var_df = pd.DataFrame({'01genid' : np.repeat(geneID, tamano),
                               '02iniex' : np.repeat(exon[4], tamano),
                               '03endex' : np.repeat(exon[5], tamano),
                               '04strnd' : np.repeat(exon[7], tamano),
                               '05frame' : np.repeat(exon[8], tamano)})
#		d. paste data to an output file with the GeneID.VAR
        exon_df = pd.concat([var_df, temp_df], axis=1)
#       this side needs to concatenate column-wise before printing each exon
        salida_df = salida_df.append(exon_df, ignore_index=True)
        #print 'salida'
        #print salida_df.head(5)
#3. Assign number codon position of subtitutions
salida_df = salida_df.drop(salida_df.head(1).index)
salida_df['16codon'] = ((abs(salida_df['07posit']-salida_df['02iniex'])+1)-salida_df['05frame'])%3
#4. Ka/Ks ratio

#   final printing
salida_df.to_csv('codonCaller_'+filename+'.var',mode = 'a', encoding = 'utf-8')
