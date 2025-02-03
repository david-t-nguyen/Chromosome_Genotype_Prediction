#import libraries
import pysam
import os
import pandas as pd
import numpy as np
import math
import collections
import sys

def get_sequences(pos, r):
    '''
    pos: position of gene we are looking for
    r: all reads
    returns sequences in arr containing the position of the gene we are looking at
    '''
    sequences = []

    for read in r:
        if pos in read.get_reference_positions():
            #collect sequence
            sequences.append(read.get_reference_sequence())

    return sequences

def ignore_wrong_reads(sequence, major, minor):
    '''
    sequence: the sequence of alleles
    major: major allele
    minor: minor allele
    returns sequence with only major or minor allele
    '''
    res = ''
    for allele in sequence:
        #only allow correct alleles
        if allele.upper() in [major, minor]:
            res += allele.upper()
    return res



def posterior_probability(sequence, major, minor, maf, e):
    '''
    sequence: sequence of alleles
    major: major allele (let it be represented by A)
    minor: minor allele (let it be represented by B)
    maf: minor allele freqency
    e: error
    returns a dictionary of the allele cominations and posterior probability
    Assumptions:
    1) Sa = 1/2
    2) PAA = (1 - maf) ^ 2
       PBB = maf ^ 2
       PAB = the rest
    '''

    #get the number of observations for major and minor
    c = collections.Counter(sequence)
    nA = c[major]
    nB = c[minor]

    Sa = 1/2

    pAA = (1-maf)**2
    pBB = maf**2
    pAB = 1 - pAA - pBB

    p_AA = (1-e)**nA * e**nB *pAA
    p_AB = ((1-e)*Sa + e*(1-Sa))**nA * (e*Sa+(1-e)*(1-Sa))**nB *pAB
    p_BB = e**nA * (1-e)**nB * pBB


    post_AA = p_AA / (p_AA + p_AB + p_BB)
    post_AB = p_AB / (p_AA + p_AB + p_BB)
    post_BB = p_BB / (p_AA + p_AB + p_BB)

    posterior_probs = {
        major+major: post_AA,
        major+minor: post_AB,
        minor+minor: post_BB
    }


    return posterior_probs

def chromosome_summary(sequences, major, minor, maf, e, name, pos):
    '''
    sequences: arr of sequences
    major: major allele (let it be represented by A)
    minor: minor allele (let it be represented by B)
    maf: minor allele freqency
    e: error
    name: name of chromosome
    pos: position of chromosome
    returns dictionary with information about the most probable genotype
    '''

    post_probs = {
        major+major: [],
        major+minor: [],
        minor+minor: []
    }

    for s in sequences:
        #get sequence with only correct alleles
        new_seq = ignore_wrong_reads(s, major, minor)
        #get posterior probabilies
        probs = posterior_probability(new_seq, major, minor, maf, e)
        for key in probs:
            post_probs[key].append(probs[key])

    combo = ''
    p = -1
    #only collect data for most probable genotype
    for key in post_probs:
        mean = sum(post_probs[key])/len(post_probs[key])
        if mean > p:
            combo = key
            p = mean
    data = {
        'chromosome': name,
        'position': pos,
        'puntative_genotype': combo,
        'posterior_probability': p,
        'n_reads': len(sequences)
    }

    return data

def chromosome_summary2(sequences, major, minor, maf, e, name, pos):
    '''
    sequences: arr of sequences
    major: major allele (let it be represented by A)
    minor: minor allele (let it be represented by B)
    maf: minor allele freqency
    e: error
    name: name of chromosome
    pos: position of chromosome
    returns array of dictionaries with information about each possible genotype
    '''

    post_probs = {
        major+major: [],
        major+minor: [],
        minor+minor: []
    }
    for s in sequences:
        new_seq = ignore_wrong_reads(s, major, minor)
        probs = posterior_probability(new_seq, major, minor, maf, e)
        for key in probs:
            post_probs[key].append(probs[key])

    data_arr = []
    for key in post_probs:
        mean = sum(post_probs[key])/len(post_probs[key])
        data = {
            'chromosome': name,
            'position': pos,
            'puntative_genotype': key,
            'posterior_probability': mean,
            'n_reads': len(sequences)
        }
        data_arr.append(data)
    return data_arr

if __name__ == "__main__":
   '''
   script assumes it will be called in the form
   python3 snps_caller.py <.bam file> <metadata .tsv file>
   '''
   #check for correct files
   if len(sys.argv) == 3:
      input_files_flag = True
      bam_file = sys.argv[1]
      tsv_file = sys.argv[2]

      if bam_file[-4:] != '.bam':
         input_files_flag = False
      if tsv_file[-4:] != '.tsv':
         input_files_flag = False

      if not input_files_flag:
         print('Incorrect input files. Default files will be used instead:')
         print('bam file: align_sort.bam')
         print('tsv file: putatative_snps.tsv')
         bam_file = 'align_sort.bam'
         tsv_file = 'putatative_snps.tsv'

   else:
      print('Wrong number of input files. Default files will be used instead:')
      print('bam file: align_sort.bam')
      print('tsv file: putatative_snps.tsv')
      bam_file = 'align_sort.bam'
      tsv_file = 'putatative_snps.tsv'

   #index bam file and then open
   curr_dir = os.getcwd()
   bam_pathway = os.path.join(curr_dir, bam_file)
   pysam.index(bam_pathway)
   samfile = pysam.AlignmentFile(bam_pathway, "rb")

   #read metadata file
   curr_dir = os.getcwd()
   puntative_snps_path = os.path.join(curr_dir, tsv_file)
   snps = pd.read_csv(puntative_snps_path, sep='\t')

   #get all the reads for the samfile
   reads = []
   for read in samfile.fetch():
       reads.append(read)

   #arr to store data to turn into dataframe
   data = []

   #collect data for each chromosome from metadata
   for i, row in snps.iterrows():
       #info from metadata
       position = row['pos']
       chr_name = row['chr']
       major = row['ref']
       minor = row['alt']
       maf = row['maf']
       e = 0.05

       sequences = get_sequences(position, reads)
       summary = chromosome_summary(sequences, major, minor, maf, e, chr_name, position)
       data.append(summary)

   #turn info into dataframe
   df = pd.DataFrame.from_dict(data)
   #output file 
   df.to_csv('output.tsv', sep="\t")


   #second output file with information about all possible genotypes
   data2 = []
   for i, row in snps.iterrows():
       #infor from metadata
       position = row['pos']
       chr_name = row['chr']
       major = row['ref']
       minor = row['alt']
       maf = row['maf']
       e = 0.05

       sequences = get_sequences(position, reads)
       summary2 = chromosome_summary2(sequences, major, minor, maf, e, chr_name, position)

       for data_genotype in summary2:
           data2.append(data_genotype)

   #turn info into dataframe
   df2 = pd.DataFrame.from_dict(data2)
   #second output file
   df2.to_csv('all_poss_genotypes_output.tsv', sep="\t")
