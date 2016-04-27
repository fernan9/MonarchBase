#! /usr/bin/env python

""" 
This script calls SNPS from simplified versions of VCF according to a gene ID
-outputs the codon position of the SNP first
-eventually will assign aminoacid change

output data should include
-- GENID	gene ID
-- SCAFF	Scaffold
-- INIEX	begin position of exon
-- ENDEX	end position of exon
-- PHASE	phase of exon
-- STRND	strand direction
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
"""

import glob, sys, csv
import subprocess
import pandas as pd

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

"""
files named .EXOME with the coordinate data
files termed .SCAF with the names of Scaffolds
files containing the EXOME and the masked SCAF are named ".VAR"
"""

# 1. First we read all the files in local directory with the extension .EXOME
exomes = glob.glob('*.exome')
paths = glob.glob('*/')

# 2. Second for each file (describing one gene by exons)
for exome in exomes:
    geneID = exome.split('.')[0]
#	A. load the file as a DataFrame with pandas
    exome-df = pd.read_table(exome, header = None)
#	B. Read each row to apply
	for exon in exome-df.itertuples():
#		a. drive a masking of the related SCAF file
# 			$ awk '{ if ($2 >= limites(0) && $2 <= limites(1)) print }' DPSCF300328.scaf
        scaf = exon[1] + '.SCAF'
        limit = (exon[4],exon[5])
        comando = "awk '{{ if ($2 >= {0} && $2 <= {1}) print }}' {2}".format(min(limit), max(limit), scaf)
        # option 1
		archivo = geneID + '.var'
        log = open(archivo, 'a')  # so that data written to it will be appended
        c = subprocess.Popen(comando, stdout=log, stderr=log, shell=True)

""" hace falta saber la salida de el comando para convertilo a un dataframe
"""	
		# option 2
		from subprocess import check_output
		out = check_output(comando.split())

#		b. create a DataFrame of it
#		c. repeat data from each Exon Row and merge it to the masked Scaffold
#		d. append data to an output file with the GeneID.VAR


#2b filter VCF by exons	
limites = {$3, $4}
ordear limites


grep -e 'DPSCF300328' Dp-var-final.vcf | awk '{ if ($2 >= 178059 && $2 <= 178170) print }'

#3. Assign number codon position of subtitutions

#4. Ka/Ks ratio
