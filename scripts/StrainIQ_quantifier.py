from Bio import SeqIO
import os
import sys
from datetime import datetime
import math
import pandas as pd

def main():
    model = "gi_ref_list_21.tsv"
    sample = "set1.fastq"
    quantifier(model, sample)
    
def quantifier(model, sample, prediction):
    ip_fastaq = sample
    ngrams_model = {}
    unique_ngrams = {}
    n=0
    #Read config file
    config = model[:-5]+".conf"
    with open(config) as config_file:
        for line in config_file:
            key = line.strip().split("=")[0]
            if key == 'n':
                n = int(line.strip().split("=")[1])
    if n==0:
        print("n-size not specified in the config file")
        sys.exit()
    
    #read predicted genomes in to a set
    prediction_pd = pd.read_csv(prediction)
    predicted_genomes = prediction_pd['gid'].tolist()
    print(predicted_genomes)

    print("Start reading model...", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    #Extract the genomes belonging to each n-gram from model file. Read model line by line
    for line in open(model):
        m_ngram, m_genomes, sq_score, cu_score = line.strip().split('\t')
        if ',' in m_genomes:
            m_genomes_set = set([e for e in m_genomes.strip().split(',')])
            intersect = m_genomes_set.intersection(predicted_genomes)
            if bool(intersect):
                if len(intersect) > 1:
                    ngrams_model[m_ngram] = ','.join(intersect)
                else:
                    unique_ngrams[m_ngram] = list(intersect)[0]
        else:
            if m_genomes.strip() in predicted_genomes:
                unique_ngrams[m_ngram] = m_genomes.strip()
    print("Completed reading model...", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

    #unique_temp_model
    #non_unique_temp_model
    with open('non_unique_model_21_noCU_predicted.tsv', 'w') as f:
        for key, value in ngrams_model.items():
            f.write('%s\t%s\n' % (key, value))
 
    with open('unique_model_21_noCU_predicted.tsv','w') as f:
        for key, value in unique_ngrams.items():
            f.write('%s\t%s\n' % (key, value))


    #Start processing unique n-grams
    print("Read unique n-grams  started...", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    print("Read fastq file and check agains unique n-grams STARTED...", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    cnt = 0
    read_genome_ngram = {}
    relative_abundance = {}

    fastq_sequences = SeqIO.parse(open(ip_fastaq),'fastq')
    for fasta in fastq_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        num_ngrams = (len(sequence)-n)+1
        read_genome_ngram[sequence] = 0
        cnt = cnt + 1
        for i in range(num_ngrams):
            temp = convert_binary(sequence[i:i+n])
            if temp in unique_ngrams:
                read_genome_ngram[sequence] = temp+"\t"+unique_ngrams[temp]
                break;
    print("Found "+str(cnt)+" reads")
    print("Read fastq file and check against unique n-grams  ENDED...", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

    print("Calculate relative abundance")
    for key, value in read_genome_ngram.items():
        ungram_genomeID = str(value).split('\t')
        if len(ungram_genomeID) == 2:
            if ungram_genomeID[1] not in relative_abundance:
                relative_abundance[ungram_genomeID[1]] = 0
            else:
                relative_abundance[ungram_genomeID[1]] = relative_abundance[ungram_genomeID[1]] + 1
    total = sum(relative_abundance.values())
    for key, value in relative_abundance.items():
        relative_abundance[key] = relative_abundance[key]/total

    relative_abundance_unique = os.path.splitext(os.path.basename(sample))[0]+".rabundance"

    with open(relative_abundance_unique, 'w') as f:
        for key, value in relative_abundance.items():
            f.write('%s\t%s\n' % (key, value))


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
