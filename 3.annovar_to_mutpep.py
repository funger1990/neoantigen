#!/usr/bin/env python
from __future__ import division

__author__ = 'Fan Zhang'


"""
edited from coding_change.pl by Kai Wang
1. ajacent SNVs
GKN2:NM_182536:G451T:G151W, GKN2:NM_182536:G463A:E155K
    1) same allele
    2) different allele

bug!!!!!!!!!!!
ENST00000453130
Exons: 4, Coding exons: 3, Transcript length: 423 bps, Translation length: 49 residues
CDS 5' incomplete
Known nonsense mediated decay

"""


import argparse
import re
from operator import itemgetter
import sys
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


WINDOW_SIZE = 15
# for indel, amino acids after mutation are all mismatches, only keep several for wild type
TRAILING_LEN = 3

genefile = '/lustre1/zeminz_pkuhpc/01.bin/annovar/humandb/hg19_ensGene.txt'
fastafile = '/lustre1/zeminz_pkuhpc/01.bin/annovar/humandb/hg19_ensGeneMrna.fa'


# ----------------------------------------
class MutPep(object):
    def __init__(self, transcript, beg, end ,obs, uniq_id, other_info):
        self.transcript = transcript
        self.beg = beg
        self.end = end
        self.obs = obs
        self.uniq_id = uniq_id
        self.other_info = other_info

    def add_info(self, aa_old, aa_new):
        self.aa_old = aa_old
        self.aa_new = aa_new


# -----------------------------------------
def load_annfile(annfile):
    queue = []
    need_trans = set()

    with open(annfile) as f:
        # ignore header
        f.readline()

        for line in f:
            field = line.strip().split('\t')
            
            genome_chr = field[0]
            genome_beg = field[1]
            gene = field[3]
            vclass = field[4]
            identifier = field[5]
            caller = field[6]
            vaf = round(float(field[7]), 2)
            sr = field[8]

            # TNFRSF9:NM_001561:exon8:c.566_567TT
            if identifier.count(':') < 4:    # synonymous
                continue

            # SLITRK4:NM_001184749:exon2:c.1449_1456del:p.N483fs
            ensgene, transcript, exon, cchange, aachange = identifier.split(':')
            cchange = cchange.split('.')[1]
            aachange = aachange.split('.')[1]

            """
                    cchange         fan-format
            snv     G11C            10,11,C
            mnv     11_12delinsGC   10,12,GC
            del     11delC          10,11,''
            del     11_13del        10,13,''
            ins     11_12insG       11,11,G
            ins     11dupC          11,11,C
                    11_15GG         10,15,GG
            """

            m = re.search(r'^([agctAGCT]?)(\d+)_?(\d*)(del)?(ins)?(dup)?([agctAGCT]*)$', cchange)
            if not m:
                sys.stderr.write('WARNING: {}\n'.format(cchange))
                print cchange
            ref, beg, end, is_del, is_ins, is_dup, obs = m.group(1, 2, 3, 4, 5, 6, 7)

            if end == '':
                end = beg
            beg = int(beg)
            end = int(end)
            if is_del and not is_ins and not is_dup:
                obs = ''

            # FAN: 0-based which convenient for downstream slicing
            if (is_dup or is_ins) and not is_del:    # insertion
                end = beg
            else:   # non-insertion
                beg -= 1
            
            uniq_id = ':'.join((transcript, cchange, aachange))
            other_info = ':'.join(map(str, (genome_chr, genome_beg, gene, vaf, sr, vclass, caller)))
            
            m = MutPep(transcript, beg, end ,obs, uniq_id, other_info)
            
            queue.append(m)
            need_trans.add(transcript)

    return queue, need_trans


def load_genefile(genefile, need_trans):
    mrnastart = {}
    mrnaend = {}

    with open(genefile) as f:
        for line in f:
            field = line.strip().split('\t')
            # first field is bin
            name, strand, txstart, txend, cdsstart, cdsend, exonstart, exonend = itemgetter(1, 3, 4, 5, 6, 7, 9, 10)(field)

            if name not in need_trans:
                continue

            #next we need to make sure that there is no intron between transcription beg and translation beg
            # (this is rare but it happens when cdsstart is not in the first exon)
            exonstart = exonstart.strip(',').split(',')
            exonend = exonend.strip(',').split(',')

            txstart = int(txstart) + 1
            txend = int(txend)
            cdsstart = int(cdsstart) + 1
            cdsend = int(cdsend)
            exonstart = [int(ii) + 1 for ii in exonstart]
            exonend = map(int, exonend)

            if strand == '+':
                intron = 0
                for i in xrange(len(exonstart)):
                    if i > 0:
                        intron += exonstart[i] - exonend[i-1] - 1
                    if cdsstart >= exonstart[i] and cdsstart <= exonend[i]:
                        mrnastart_1 = cdsstart - txstart + 1 - intron
                    if cdsend >= exonstart[i] and cdsend <= exonend[i]:
                        mrnaend_1 = cdsend - txstart + 1 - intron
            elif strand == '-':
                intron = 0
                for i in reversed(xrange(len(exonstart))):
                    if i < len(exonstart) - 1:
                        intron += exonstart[i+1] - exonend[i] - 1
                    if cdsend >= exonstart[i] and cdsend <= exonend[i]:
                        mrnastart_1 = txend - cdsend + 1 - intron
                    if cdsstart >= exonstart[i] and cdsstart <= exonend[i]:
                        mrnaend_1 = txend - cdsstart + 1 - intron

            mrnastart[name] = mrnastart_1
            mrnaend[name] = mrnaend_1

    return mrnastart, mrnaend


def load_fastafile(fastafile, need_trans):
    mrnaseq = {}
    with open(fastafile) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            if seq_record.id in need_trans:
                mrnaseq[seq_record.id] = str(seq_record.seq)

    return mrnaseq


def process_mutation(mutpep, mrnastart, mrnaend, mrnaseq):
    transcript = mutpep.transcript
    beg = mutpep.beg
    end = mutpep.end
    obs = mutpep.obs
    
    if transcript not in mrnaseq or transcript not in mrnastart or transcript not in mrnaend:
        return None

    utr5 = mrnaseq[transcript][: mrnastart[transcript] - 1]
    utr3 = mrnaseq[transcript][mrnaend[transcript] :]

    dna = mrnaseq[transcript][mrnastart[transcript] - 1 : mrnaend[transcript]]

    if end > len(dna):
        return None
    
    # explicitly trim dna to 3 multiples
    dna = dna[: len(dna) // 3 * 3]
    protein1 = str(Seq(dna, generic_dna).translate(to_stop=True))

    # FAN: mutant. customized-defined beg, end
    """
    bug!!!!!!!!!!
    for conservation, don't add utr3

    """
    dna2 = dna[:beg] + obs + dna[end:]
    dna2 = dna2[: len(dna2) // 3 * 3]
    protein2 = str(Seq(dna2, generic_dna).translate(to_stop=True))

    # synonymous
    if protein1 == protein2:
        return None
    
    """
    mutpep format
                wt      mut
    snv     15[F]15     15[Y]15
    mnv     15[FF]15    15[YY]15
    inframe 15[]15      15[Y]15 
            15[F]15     15[]15
    frameshift  15[FFF] 15[YYY]
    stoploss    15[]    15[YYY]          

    """

    # SNV or inframe indel
    if (beg - end + len(obs)) % 3 == 0:
        # python index
        aa_beg = beg // 3
        aa_end = (end - 1) // 3 + 1
        aa_end_new = (beg + len(obs) - 1) // 3 + 1

        """
        bug. a = 'abcdef'; a[-1:1] == '';
        """
        protein1_window = protein1[max(0, aa_beg - WINDOW_SIZE) : aa_beg] + \
                '[' + protein1[aa_beg : aa_end] + ']' + \
                protein1[aa_end : aa_end + WINDOW_SIZE]

        protein2_window = protein2[max(0, aa_beg - WINDOW_SIZE) : aa_beg] + \
                '[' + protein2[aa_beg : aa_end_new] + ']' + \
                protein2[aa_end_new : aa_end_new + WINDOW_SIZE]
        
    else:
        # frameshift or stoploss
        aa_beg = beg // 3

        # max length of WT is the same as mutated peptide
        protein1_window = protein1[max(0, aa_beg - WINDOW_SIZE) : aa_beg] + \
                '[' + protein1[aa_beg : len(protein2)] + ']'

        protein2_window = protein2[max(0, aa_beg - WINDOW_SIZE) : aa_beg] + \
                '[' + protein2[aa_beg :] + ']'

    mutpep.add_info(protein1_window, protein2_window)

    return mutpep


def main(genefile, fastafile, annfile, outfile):
    # load annovar file
    queue, need_trans = load_annfile(annfile)

    # load genefile, fastafile (only part of file)
    mrnastart, mrnaend = load_genefile(genefile, need_trans)
    mrnaseq = load_fastafile(fastafile, need_trans)
    
    # process each element
    queue_valid = []
    for mutpep in queue:
        mutpep = process_mutation(mutpep, mrnastart, mrnaend, mrnaseq)
        if mutpep is not None:
            queue_valid.append(mutpep)

    # output
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    with open(outfile, 'w') as f:
        for m in queue_valid:
            f.write('{}\t{}\t{}\t{}\n'.
                    format(m.uniq_id, m.aa_old, m.aa_new, m.other_info))


# -----------------------------------------------
if __name__ == '__main__':
    annfile = sys.argv[1]
    outfile = sys.argv[2]

    main(genefile, fastafile, annfile, outfile)
