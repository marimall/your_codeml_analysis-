import numpy as np
import os
import pandas as pd
import csv
import scipy.stats as stats
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import Phylo
from io import StringIO
from Bio.Seq import Seq
from itertools import groupby, tee
import mygene
from pyfaidx import Fasta


def parse_args():   #### for all the arguments the user is goin to provide
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-fasta", default="/dev/stdin")
    parser.add_argument("--output-phy", default="/dev/stdout")
    parser.add_argument("--gene")

    return parser.parse_args()

def fasta_iter(fasta_name): 

######reads file line by line instead of loading it all in memory

#"first open the file outside "
    fh = open(fasta_name)

# ditch the boolean (x[0]) and just keep the header or sequence since
# we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
    # drop the ">"
        headerStr = header.__next__()[1:].strip()

    # join all sequence lines to one.
    
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

def gene_finder(genes, my_genes):
    desired_genes = {}
    mg = mygene.MyGeneInfo()
    ids =  mg.querymany(genes, scopes='ensembl.transcript')
    return ids

  ##  if len(ids)>2:
    ##    desired_genes[ids["querry"]] = ids["symbol"]
    
    ##return desired_genes

    
def fasta_to_phyllip(fasta, phyllip) :
    
    with open(fasta) as handle:
        records = AlignIO.parse(handle, "fasta")

        with open(phyllip, "w") as output_handle:
            AlignIO.write(records, output_handle, "phylip-sequential")
            
def manual_fasta_to_phyllip(fasta, phyllip) :
    species = []
    seqs = []

    with open(fasta, "r") as fp:
        lines = fp.readlines()
        
        for l in range(0, len(lines), 2) :
            sp , seq = lines[l], lines[l+1]
            species.append(sp.strip()[1:])
            seqs.append(seq.strip())
    
        
    one =  str(len(species))
    two =  len(seqs[0])
    
   
       
    if two%3 == 0 :
        with open(phyllip, "w") as ff:
            ff.write(" " + one + " " + str(two) + "\n")
            for r in range(len(species)) :
                ff.write(species[r] + "  " + seqs[r] + "\n")
    elif two%3 == 1 :
        with open(phyllip, "w") as ff:
            ff.write(" " + one + " " + str(two -1) + "\n")
            for r in range(len(species)) :
                ff.write(species[r] + "  " + seqs[r] + "\n")
    else:
        with open(phyllip, "w") as ff:
            ff.write(" " + one + " " + str(two - 2) + "\n")
            for r in range(len(species)) :
                ff.write(species[r] + "  " + seqs[r] + "\n")
				
				
def translate(seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein


def write_codeml(file, inputf, tree, results, model ) :
    with open(file, "r") as f :
        content = f.readlines()
    content[0] = '      seqfile = ' + inputf + ' * sequence data filename\n'
    content[1] = '     treefile = ' + tree +'      * tree structure file name\n'
    content[2] = '      outfile = ' + results + '           * main result file name\n'
    content[18] = '        model = {}\n'.format(model)
    with open(file, 'w') as filehandle:
        filehandle.writelines("%s" % place for place in content)


########## -------- main ----------------------------

fiter = fasta_iter("knownCanonical.exonNuc.fa") ## later change to args.input_fasta

fiter, fiter_tee = tee(fiter, 2) ##### get a duplicate generator to regenerate fasta an isolate genes of our prefference
#fiter_tee = list(fiter_tee)
gene_list = ["GRIN"+i for i in ["1","2A", "2B", "2C", "2D", "3A", "3B"]] ####args.gene_list
desired_genes = []

genes = []
for ff in fiter:
    headerStr, seq = ff
    head = headerStr.split(".")[0]
    genes.append(head)

all_genes = list(set(genes))    

trancripts = gene_finder(all_genes, gene_list) #### find all ensembl trancripts of your genes

a_dict = {}
for g in trancripts:    #######creates a dictionary of ensembl names and gene names
    if "symbol" in g.keys() :
        #if g["symbol"] in gene_list:
        a_dict[g['query']] = g["symbol"] 


####find all trancripts of your gene list

all_trans= { d : "" for d in gene_list }

for ff in fiter_tee:
    headerStr, seq = ff
    head = headerStr.split(".")[0]
    #print(head)
    if head in a_dict:
        if a_dict[head] in all_trans.keys() :
            all_trans[a_dict[head]] = all_trans[a_dict[head]] + headerStr + "\n" + seq + "\n"

if os.path.isdir("trancripts_to_fasta") == False:
    os.mkdir( "trancripts_to_fasta" ) 
    #os.system("mkdir trancripts_to_fasta")

for d in all_trans:   ####### creates fasta files with each trancript

    tran = all_trans[d].split("\n")
    number = 1
    species = tran[0].split("_")[1]
    filename = open( d + "_" + str(1) + ".fa" , "w" )
    filename.write(">" + species +"\n" + tran[1] +"\n")
    
    for t in range(2,len(tran)-1,2):
        species = tran[t].split("_")[1]
   
        if species == "hg38" :
            filename.close()
            number +=1 
            filename =open( d + "_" + str(number) + ".fa"  , "w" )
            #tran[t] = "<" + tran[t]
            species = tran[t].split("_")[1]
            filename.write(">" + species+"\n" + tran[t+1] +"\n")
        
        elif species == "odoRosDiv1" :
            filename.write(">" + species+"  \n" + tran[t+1] +"\n")
        else :
            filename.write(">" + species+"\n" + tran[t+1] +"\n")
     
        
    filename.close()
    
all_files = os.listdir()   ###### later replace it by creating specific file that contains 


if os.path.isdir("phyllip_files") == False:
    os.mkdir( "phyllip_files" )  
    #os.system("mkdir phyllip_files")
    
for af in all_files :
    if af.startswith("GRIN") and af.endswith(".fa") :
        #print(af)
        manual_fasta_to_phyllip(af, af[0:-3] + "m.phy") ####later change output file
        
  
if os.path.isdir("codeml_results") == False:
    os.mkdir( "codeml_results" )  
    #os.system("mkdir codeml_results")

#### four runs for each model maybe change so the user can specify the models

for af in all_files :
    if af.startswith("GRIN") and af.endswith(".phy") :
        write_codeml("codeml.ctl", af ,"hg38.H0.nh", "{}_equal_omega.res".format(af[0:-4]), 0 )
        os.system("codeml")
        write_codeml("codeml.ctl", af ,"hg38.H1.nh", "{}_res2_onebranch.res".format(af[0:-4]), 2 )
        os.system("codeml")
        write_codeml("codeml.ctl", af ,"hg38H2.nh", "{}_res2_allbranches.res".format(af[0:-4]), 2 )
        os.system("codeml") 
        write_codeml("codeml.ctl", af ,"hg38.H3.nh", "{}_res2hg.res".format(af[0:-4]), 2 )
        os.system("codeml")


    
all_files = os.listdir()
lnl = {}
trees = {}
for af1 in all_files :
    if af1.startswith("GRIN") and af1.endswith(".res") :
        with open(af1) as f:    
            lines = [line.rstrip() for line in f]
        starts = [n for n, l in enumerate(lines) if l.startswith('tree length')]
        trees[af1] = lines[starts[0]+4]
        lnl[af1] = [q for q in lines if q.startswith("lnL")]


lnlcor = { i : lnl[i][0].split(":")[3].split("     ")[0] for i in lnl }
columns = [i.split("_")[2] for i in lnlcor.keys() ]
rows_names = [i.split("_")[0] +"_" +i.split("_")[1] for i in lnlcor]


df = pd.DataFrame(columns = set(columns), index = set(rows_names))
df.to_csv("all_likelihoods.txt")


