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

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code



# 1. First we read all the files in local directory with the extension .EXOME
exomes = glob.glob('*.exome')
paths = glob.glob('*/')

#print type(exomes)

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
#        comando = ["awk", "'{{ if ($2 >= {0} && $2 <= {1}) print }}'".format(min(limit),max(limit)), "{0}".format(scaf)]
        print limit
        try:
            salida = subprocess.check_output(comando, stderr=subprocess.STDOUT, shell = True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

"""
#		b. create a DataFrame of it
        temp_df = pd.read_table(StringIO(salida), header = None)
#		c. repeat data from each Exon Row and merge it to the masked Scaffold
        d = {'GENID' : np.repeat(geneID, len(df)),
            'INIEX' : np.repeat(exon(4), len(df)),
            'ENDEX' : np.repeat(exon(5), len(df)),
            'STRND' : np.repeat(exon(7), len(df)),
            'PHASE' : np.repeat(exon(8), len(df))}
        var_df = pd.DataFrame(d)
#		d. paste data to an output file with the GeneID.VAR
        result = pd.concat([var_df, temp_df], axis=1)
        df.to_csv(geneID+'.var', encoding = 'utf-8')
        print result.head(10)
#3. Assign number codon position of subtitutions

#4. Ka/Ks ratio
"""
