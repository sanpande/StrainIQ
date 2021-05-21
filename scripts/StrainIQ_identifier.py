from Bio import SeqIO
import os
import sys
import datetime
import math
import pandas as pd

def main():
    model = "gi_ref_list_21.tsv"
    sample = "set1.fastq"
    identifier(model, sample)

def identifier(model, sample):
   b_ngrams = set()
   score={} #square score
   score2={} #cube score
   cNgram={}
   #thisFile = os.path.realpath(__file__)
   #cd = os.path.dirname(thisFile)

   #Read config file
   config = model[:-5]+".conf"
   with open(config) as config_file:
       for line in config_file:
           key = line.strip().split("=")[0]
           if key == 'n':
               n = int(line.strip().split("=")[1])
           if key == 'number of bacteria':
               no_genome = int(line.strip().split("=")[1])
           if key == 'cutoff':
               cutoff =  float(line.strip().split("=")[1])
   #initialize score matrix
   for i in range(1,no_genome + 1):
      score[i]=0
      score2[i]=0
      cNgram[i]=0
   # initialize score matrix END
   fasta_sequences = SeqIO.parse(open(sample),'fastq')
   for fasta in fasta_sequences:
       name, sequence = fasta.id, str(fasta.seq)
       num_ngrams = (len(sequence)-n)+1
       for i in range(num_ngrams):
           temp = convert_binary(sequence[i:i+n])      
           if(len(temp) == n*2):
               b_ngrams.add(temp)

   for line in open(model):
      m_ngram, m_genomes, m_score, m_score2 = line.strip().split('\t')
      if m_ngram in b_ngrams:
         for g in m_genomes.split(','):
            score[int(g)]=score[int(g)]+float(m_score)
            score2[int(g)]=score2[int(g)]+float(m_score2)
            cNgram[int(g)]=cNgram[int(g)]+1

   score_file = os.path.splitext(os.path.basename(sample))[0]+".score"
   fOut = open(score_file,"w")
   for g,s in score.items():
      fOut.write("%s\t%s\t%s\n" %(g,s,cNgram[g]))
   fOut.close()

   #Use cutoff value to predict strains
   id_map = model[:-5]+".map"
   taxa_map = model[:-5]+".taxa"
   score_pd = pd.read_csv(score_file, sep='\t', names=['genomeID', 'score', 'cNgram'])
   map_pd = pd.read_csv(id_map, sep='\t')
   complete_pd = pd.merge(left=score_pd, right=map_pd, left_on='genomeID', right_on='genomeID')
   #normalize using nFactor
   complete_pd["adjScore"] = complete_pd["score"]/complete_pd["uNgrams"]*pow(complete_pd["cNgram"]/complete_pd["uNgrams"],4)
   prediction_pd = complete_pd[complete_pd["adjScore"] >= cutoff].sort_values(by=['adjScore'])['genomeID']
   taxa_prediction = pd.read_csv(taxa_map, sep='\t')
   prediction_combined = pd.merge(left=prediction_pd, right=taxa_prediction, left_on='genomeID', right_on='gid')
   prediction_file = os.path.splitext(os.path.basename(sample))[0]+".prediction"
   prediction_combined[['gid','RefSeqassemblyaccession','strain','phylum','class','order','family','genus','species','Organism.Name']].to_csv(prediction_file, index=False)


def convert_binary(str):
   ret_str = ""
   for amino in str.upper():
      if amino == "A":
         ret_str = "".join([ret_str,"00"])
      elif amino == "C":
         ret_str = "".join([ret_str,"01"])
      elif amino == "G":
         ret_str = "".join([ret_str,"10"])
      elif amino == "T":
         ret_str = "".join([ret_str,"11"])
   return ret_str
  
if __name__ == "__main__":
    main()
