#!/usr/bin/python
import csv
from Bio import SeqIO
import os
import sys
import json
import pickle
import datetime
import math

def main():
   model(21,"gi_mock_gutPlusMixture_ref_list.tsv")


def model(n, b_list):
   ngrams = {}
   number_of_bacteria = 0
   model = os.path.splitext(os.path.basename(b_list))[0] +'_'+str(n)+'.dsem'
   conf = os.path.splitext(os.path.basename(b_list))[0] +'_'+str(n)+'.conf'
   with open(b_list) as b:
      for line in b:
         file = line.rstrip().split("\t")
         b_id = file[0]
         print("processing ",file[1])
         b_ngrams = set()
         fasta_sequences = SeqIO.parse(open(file[1]),'fasta')
         for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            num_ngrams = (len(sequence)-n)+1
            for i in range(num_ngrams):
               temp = convert_binary(sequence[i:i+n])
               if(len(temp) == n*2):
                   b_ngrams.add(temp)
         for bin_item in b_ngrams:
            if bin_item in ngrams:
               ngrams[bin_item] = ','.join([ngrams[bin_item],b_id])
            elif (bin_item != "") :
               ngrams[bin_item] = b_id
         b_ngrams.clear()
         number_of_bacteria = number_of_bacteria + 1
   #Write configuration file
   with open(conf,'w') as conf_file:
       conf_file.write("n="+str(n)+"\nnumber of bacteria="+str(number_of_bacteria))
   #Write dsem
   writer = csv.writer(open(model, 'w'), delimiter='\t')
   start_time = datetime.datetime.now()
   for key, value in ngrams.items():
      num_genomes = value.count(",")+1
      score = abs(calculate_score(number_of_bacteria, num_genomes))
      writer.writerow([key, value, math.pow(float(score),2),math.pow(float(score),3)])

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
  
def calculate_score(total_genomes, ngram_genomes):
    return math.log(abs(total_genomes/ngram_genomes))/math.log(abs(total_genomes))

if __name__ == "__main__":
    main()
