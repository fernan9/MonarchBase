#! /usr/bin/env python

"""
---DESCRIPTION
This script calls SNPS from simplified versions of VCF according to a gene ID
-outputs the codon position of the SNP first
-eventually will assign aminoacid change

---OUTPUT format
output data should include
-- GENID	gene ID
-- INIEX	begin position of exon
-- ENDEX	end position of exon
-- STRND	strand direction
-- PHASE	phase of exon
-- SCAFF	Scaffold
-- ANSNP	ancient SNP
-- VASNP	variable SNP
-- CDPOS	codon position of snip
-- BPROB	associated base probability
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

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

# 1. First we read all the files in local directory with the extension .EXOME
exomes = glob.glob('*.exome')
paths = glob.glob('*/')

# 2. Second for each file (describing one gene by exons)
for exome in exomes:
    geneID = exome.split('.')[0]
    print geneID
#	A. load the file as a DataFrame with pandas
    exome_df = pd.read_table(exome, header = None)
#	B. Read each row to apply

###   probably we need to define here the df for the complete gene

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
        tamano = len(temp_df.index)
#		c. repeat data from each Exon Row and merge it to the masked Scaffold
        var_df = pd.DataFrame({'01genid' : np.repeat(geneID, tamano),
                               '02iniex' : np.repeat(exon[4], tamano),
                               '03endex' : np.repeat(exon[5], tamano),
                               '04strnd' : np.repeat(exon[7], tamano),
                               '04frame' : np.repeat(exon[8], tamano)})
#		d. paste data to an output file with the GeneID.VAR
        result = pd.concat([var_df, temp_df], axis=1)
#       this side needs to concatenate column-wise before printing each exon
        result.to_csv(geneID+'.var',mode = 'a', encoding = 'utf-8')
        print result.head(10)
#3. Assign number codon position of subtitutions
"""
#4. Ka/Ks ratio
"""
