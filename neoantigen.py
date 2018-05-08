#!/usr/bin/env python
from __future__ import division
__author__ = 'Fan Zhang'

import re
import os
import sys
import gzip
import math
from collections import defaultdict
from operator import itemgetter
import numpy as np
import pandas as pd
import subprocess
from tempfile import NamedTemporaryFile
from multiprocessing import Pool
from functools import partial

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


WINDOW_SIZE = 15
TRAILING_LEN = 3
TMPDIR = '/lustre1/zeminz_pkuhpc/fan/test'

file_mrna_anno = '/lustre1/zeminz_pkuhpc/01.bin/annovar/humandb/hg19_ensGene.txt'
file_mrna = '/lustre1/zeminz_pkuhpc/01.bin/annovar/humandb/hg19_ensGeneMrna.fa'
file_pro = '/lustre1/zeminz_pkuhpc/00.database/gencode/release_19/gencode.v19.pc_translations.fa.gz'

pan_version = 4
#dir_pan = '/lustre1/zeminz_pkuhpc/01.bin/netMHC/netMHCpan-3.0'
dir_pan = '/lustre1/zeminz_pkuhpc/01.bin/netMHC/netMHCpan-4.0'


# ==================================================
def load_protein(enst_list):
    d = dict()
    with gzip.open(file_pro) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            enst = seq_record.id.split('.')[0]
            if enst in enst_list:
                d[enst] = str(seq_record.seq)

    return d


def change_aa(row, protein_ref):
    tx = row['enst']
    aachange = row['aachange']
    m = re.search(r'([A-Z])(\d+)([A-Z])', aachange)
    ref, beg, alt = m.group(1, 2, 3)
    beg = int(beg) - 1
    end = beg + 1

    pro1_raw = protein_ref[tx]

    if pro1_raw[beg] == ref:
        pro2_raw = pro1_raw[: beg] + alt + pro1_raw[beg + 1 :]
    else:
        print '[wrong annotation]', tx, aachange
        return None, None

    # left index must be at least 0
    pro1 = pro1_raw[max(0, beg - WINDOW_SIZE) : beg] + \
            '[' + pro1_raw[beg] + ']' + \
            pro1_raw[end : end + WINDOW_SIZE]
    pro2 = pro2_raw[max(0, beg - WINDOW_SIZE) : beg] + \
            '[' + pro2_raw[beg] + ']' + \
            pro2_raw[end : end + WINDOW_SIZE]

    return pro1, pro2


def parse_cchange(cchange):
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

    m = re.search(
        r'^([agctAGCT]?)(\d+)_?(\d*)(del)?(ins)?(dup)?([agctAGCT]*)$',
        cchange)
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

    return beg, end, obs


def load_mrna(enst_list):
    d = defaultdict(dict)

    with open(file_mrna) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            enst = seq_record.id
            if enst in enst_list:
                d[enst]['seq'] = str(seq_record.seq)

    with open(file_mrna_anno) as f:
        for line in f:
            field = line.strip().split('\t')
            # first field is bin
            name, strand, txstart, txend, cdsstart, cdsend, exonstart, exonend = itemgetter(1, 3, 4, 5, 6, 7, 9, 10)(field)

            if name not in enst_list:
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

            d[name]['beg'] = mrnastart_1
            d[name]['end'] = mrnaend_1

    return d


def translate_mut(row, mrna_ref):
    tx = row['enst']
    beg = row['beg']
    end = row['end']
    obs = row['obs']
    tx_seq = mrna_ref[tx]['seq']
    tx_beg = mrna_ref[tx]['beg']
    tx_end = mrna_ref[tx]['end']

    if not tx_seq or not tx_beg or not tx_end:
        print '[No transcript]', tx
        return None, None

    utr5 = tx_seq[: tx_beg - 1]
    utr3 = tx_seq[tx_end :]

    dna = tx_seq[tx_beg - 1 : tx_end]

    if end > len(dna):
        return None, None

    # explicitly trim dna to 3 multiples
    dna = dna[: len(dna) // 3 * 3]
    pro1_raw = str(Seq(dna, generic_dna).translate(to_stop=True))

    # FAN: mutant. customized-defined beg, end
    dna2 = dna[:beg] + obs + dna[end:] + utr3
    dna2 = dna2[: len(dna2) // 3 * 3]
    pro2_raw = str(Seq(dna2, generic_dna).translate(to_stop=True))

    # synonymous indel
    if pro1_raw == pro2_raw:
        print '[total same]', tx
        return None, None

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

    if (beg - end + len(obs)) % 3 == 0:
        # SNV, MNV or inframe indel
        aa_beg = beg // 3
        aa_end = (end - 1) // 3 + 1
        aa_end_new = (beg + len(obs) - 1) // 3 + 1
        # left index must be at least 0
        pro1 = pro1_raw[max(0, aa_beg - WINDOW_SIZE) : aa_beg] + \
                '[' + pro1_raw[aa_beg : aa_end] + ']' + \
                pro1_raw[aa_end : aa_end + WINDOW_SIZE]
        pro2 = pro2_raw[max(0, aa_beg - WINDOW_SIZE) : aa_beg] + \
                '[' + pro2_raw[aa_beg : aa_end_new] + ']' + \
                pro2_raw[aa_end_new : aa_end_new + WINDOW_SIZE]
    else:
        # frameshift or stoploss
        aa_beg = beg // 3
        pro1 = pro1_raw[aa_beg - WINDOW_SIZE : aa_beg] + \
                '[' + pro1_raw[aa_beg : aa_beg + TRAILING_LEN] + ']'
        pro2 = pro2_raw[aa_beg - WINDOW_SIZE : aa_beg] + \
                '[' + pro2_raw[aa_beg :] + ']'

    if pro2.endswith('[]'):
        print '[stopgain]', tx
        return None, None
    elif pro2.startswith('[]'):
        print '[beginning loss?]', tx
        return None, None

    # TODO
    i1 = pro1.index('[') + 1
    j1 = pro1.index(']')
    i2 = pro2.index('[') + 1
    j2 = pro2.index(']')

    def swap(s, i, j):
        t = list(s)
        t[i], t[j] = t[j], t[i]
        s = ''.join(t)
        return s

    while i1 < j1 and i2 < j2:
        if pro1[i1] == pro2[i2]:
            print '[leading same]', tx, pro1, pro2
            pro1 = swap(pro1, i1 - 1, i1)
            pro2 = swap(pro2, i2 - 1, i2)
            i1 += 1
            i2 += 1
        else:
            break

    if i1 == j1 and i2 == j2:
        print '[frameshift same]', tx
        return None, None

    return pro1, pro2


def mutation_to_protein(df):
    # uid, enst, cchange, aachange
    df['substitute'] = df['aachange'].map(
        lambda x: bool(re.search(r'[A-WYZ]\d+[A-WYZ]', x)))

    # 1) single aa change
    df1 = df[df['substitute']]
    if not df1.empty:
        protein_ref = load_protein(df1['enst'].unique())
        df1['old'], df1['new'] = zip(* df1.apply(
            lambda x: change_aa(x, protein_ref), axis=1))
    else:
        df1['old'] = np.nan
        df1['new'] = np.nan

    # 2) complex aa change
    df2 = df[~df['substitute']]
    if not df2.empty:
        mrna_ref = load_mrna(df2['enst'].unique())
        df2['beg'], df2['end'], df2['obs'] = zip(* df2['cchange'].map(parse_cchange))
        df2['old'], df2['new'] = zip(* df2.apply(
            lambda x: translate_mut(x, mrna_ref), axis=1))
    else:
        df2['old'] = np.nan
        df2['new'] = np.nan

    df = pd.concat([df1, df2])
    df = df[['uid', 'enst', 'cchange', 'aachange', 'old', 'new']]
    df.fillna('-', inplace=True)

    return df


# =============================================
def chop_one(uid, old, new, k_list):
    idx_old1 = old.find('[')
    idx_old2 = old.find(']')
    idx1 = new.find('[')
    idx2 = new.find(']')

    if len(old) == len(new) and old[: idx_old1] == new[: idx1] and \
            old[idx_old2 :] == new[idx2 :]:
        has_wt = True
    else:
        has_wt = False

    # remove '[' and ']'
    old = old.replace('[', '').replace(']', '')
    idx_old2 -= 1
    new = new.replace('[', '').replace(']', '')
    idx2 -= 1

    result = []
    for k in k_list:
        for i in xrange(max(0, idx1 - k + 1), min(len(new) - k + 1, idx2)):
            pep_new = new[i : i+k]
            if has_wt:
                pep_old = old[i : i+k]
            else:
                pep_old = '-'
            result.append([uid, pep_old, pep_new])

    return result


def chop_kmer(df, k_list=[9, 10, 11]):
    # uid, old, new
    result = []
    for i, row in df.iterrows():
        result_1 = chop_one(row['uid'], row['old'], row['new'], k_list)
        result += result_1
    df = pd.DataFrame(result, columns=['uid', 'pep_old', 'pep_new'])

    return df


# ==============================================
def reformat_hla(hla):
    hla = hla.upper()
    m = re.search(r'(HLA-)?([ABC]\*\d+:\d+)(:?\d*)[NLSCAQ]?', hla)
    if m:
        hla = m.group(2)
        hla = 'HLA-' + hla.replace('*', '')
    else:
        print '[wrong hla]', hla
        hla = None

    return hla


def predict_netmhc_one(pep_list, hla):
    with open(os.path.join(dir_pan, 'data', 'allelenames')) as f:
        valid_hla = [line.strip().split()[0] for line in f]

    hla_pan = reformat_hla(hla)
    if hla_pan not in valid_hla:
        print '[not in netmhcpan] ', hla
        return None

    with NamedTemporaryFile() as pepfile:
        with open(pepfile.name, 'w') as f:
            f.write('\n'.join(pep_list) + '\n')

        if pan_version == 4:
            p = subprocess.Popen([
                os.path.join(dir_pan, 'netMHCpan'),
                '-a', hla_pan,
                '-p', pepfile.name,
                '-BA'
            ], stdout=subprocess.PIPE)
        elif pan_version == 3:
            p = subprocess.Popen([
                os.path.join(dir_pan, 'netMHCpan'),
                '-a', hla_pan,
                '-p', pepfile.name,
            ], stdout=subprocess.PIPE)

        out = p.communicate()[0]

    result = []
    for line in out.splitlines():
        if 'PEPLIST' in line and 'high binders' not in line:
            field = line.strip().split()
            pep = field[2]
            ic50 = field[12]
            rank = field[13]
            result.append([hla, pep, ic50, rank])
    df = pd.DataFrame(result, columns=['hla', 'pep', 'ic50', 'rank'])

    return df


def predict_netmhc(pep_list, hla_list, thread=1):
    block = int(math.ceil(len(pep_list) / thread))
    pep_list_split = [pep_list[i * block : (i+1) * block] for i in range(thread)]
    list_df = []

    for hla in hla_list:
        if thread == 1:
            list_df.append(predict_netmhc_one(pep_list, hla))
        else:
            p = Pool(thread)
            list_df += p.map(partial(predict_netmhc_one, hla=hla),
                             pep_list_split)

    df = pd.concat(list_df)
    df['rank'] = df['rank'].astype(float)
    df['ic50'] = df['ic50'].astype(float)

    return df


# ===============================================
if __name__ == '__main__':
    #df = pd.DataFrame({
    #    'uid': [1, 2],
    #    'enst': ['ENST00000256078', 'ENST00000357668'],
    #    'cchange': ['-', '2818_2838del'],
    #    'aachange': ['G12D', '940_946del']
    #})
    #df = mutation_to_protein(df.iloc[[1], :])

    #df = pd.DataFrame({
    #    'uid': ['haha'],
    #    'old': ['KKFISQIVDTLDV[SC]DKLAQVGLVQYSSSV'],
    #    'new': ['KKFISQIVDTLDV[PA]DKLAQVGLVQYSSSV']
    #})
    #df = chop_kmer(df, [9, 10, 11])

    print 'netMHCpan', pan_version
    df = predict_netmhc(
        ['KIGDFGLATEK', 'FGLATEKSRW', 'EKSRWSGSHQF', 'LATEKSRWSG', 'IGDFGLATEK'],
        ['A*03:25', 'A*24:02', 'B*15:71', 'B*35:02', 'C*04:01', 'C*07:01'],
        1)

    print df

