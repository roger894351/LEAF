#!/usr/bin/env python
import os
import sys
import re
import zlib
import numpy
import random
import subprocess


def load_mask_bed_file(mask_bed_file):
    mask_dict = {}
    if os.path.exists(mask_bed_file):
        FH = open(mask_bed_file)
        for line in FH:
            line = line.rstrip()
            t = line.split('\t')
            for x in range(int(t[1]), int(t[2])+1):
                if t[0] not in mask_dict:
                    mask_dict[t[0]] = {}
                mask_dict[t[0]][str(x)] = 0
        FH.close()
    return mask_dict


def load_gff3_file(gff3_file):
    mRNA_dict = {}
    mRNA2gene = {}
    FH = open(gff3_file)
    for line in FH:
        line = line.rstrip()
        if not line.startswith("#") and line != '':
            t = line.split('\t')
            if t[2] == 'mRNA':
                gene = None
                mRNA = None
                for x in t[-1].split(';'):
                    if x.startswith('ID='):
                        mRNA = x[3:]
                    if x.startswith('Parent='):
                        gene = x[7:]
                if gene and mRNA:
                    mRNA2gene[mRNA] = gene
            if t[2] == 'CDS':
                mRNA = None
                for x in t[-1].split(';'):
                    if x.startswith('Parent='):
                        mRNA = x[7:]
                if mRNA:
                    if mRNA not in mRNA_dict:
                        mRNA_dict[mRNA] = {'chr': t[0], 'region': {}}
                    mRNA_dict[mRNA]['region'][int(t[3])] = int(t[4])
    FH.close()
    temp_gff3_annotation = job_id+".gff3_annotation.bed"
    OH = open(temp_gff3_annotation, "w")
    for mRNA in mRNA_dict:
        cds_starts = mRNA_dict[mRNA]['region'].keys()
        cds_starts.sort()
        for cds_start in cds_starts:
            text = [mRNA_dict[mRNA]['chr'], str(cds_start), str(
                mRNA_dict[mRNA]['region'][cds_start]), mRNA, mRNA2gene[mRNA], 'CDS']
            OH.write('\t'.join(text)+'\n')
        text = [mRNA_dict[mRNA]['chr'], str(
                mRNA_dict[mRNA]['region'][cds_starts[0]]+1), str(
                mRNA_dict[mRNA]['region'][cds_starts[0]]+2), mRNA, mRNA2gene[mRNA], 'SPLICE_A/D']
        OH.write('\t'.join(text)+'\n')
        for cds_start in cds_starts[1:-1]:
            text = [mRNA_dict[mRNA]['chr'], str(
                cds_start-2), str(cds_start-1), mRNA, mRNA2gene[mRNA], 'SPLICE_A/D']
            OH.write('\t'.join(text)+'\n')
            text = [mRNA_dict[mRNA]['chr'], str(
                mRNA_dict[mRNA]['region'][cds_start]+1), str(
                mRNA_dict[mRNA]['region'][cds_start]+2), mRNA, mRNA2gene[mRNA], 'SPLICE_A/D']
            OH.write('\t'.join(text)+'\n')
        text = [mRNA_dict[mRNA]['chr'], str(
            cds_starts[-1]-2), str(cds_starts[-1]-1), mRNA, mRNA2gene[mRNA], 'SPLICE_A/D']
        OH.write('\t'.join(text)+'\n')
    OH.close()
    return temp_gff3_annotation


def readlines_reverse(input_file, buffer=1024*1024):
    with open(input_file, 'r') as somefile:
        somefile.seek(0, os.SEEK_END)
        size = somefile.tell()
        lines = ['']
        rem = size % buffer
        pos = max(0, (size // buffer - 1) * buffer)
        while pos >= 0:
            somefile.seek(pos, os.SEEK_SET)
            data = somefile.read(rem + buffer) + lines[0]
            rem = 0
            lines = re.findall('[^\n]*\n?', data)
            ix = len(lines) - 2
            while ix > 0:
                yield lines[ix]
                ix -= 1
            pos -= buffer
        else:
            yield lines[0]


def list_text2array(input_text):
    sample_array = []
    for x in input_text.split(','):
        if '-' not in x:
            sample_array.append(int(x))
        else:
            s = x.split('-')
            for y in range(int(s[0]), int(s[1])+1):
                sample_array.append(y)
    return sample_array


def parse_cdr_line(cdr_line, sample_size):
    t = cdr_line.split('\t')
    ref = t[5].split("|")[0]
    loci = t[0:4]
    array_collection = []
    alternative_collection = []
    gt_collection = [x.split('|')[1].split(':') for x in t[6:]]
    for x in gt_collection:
        for xx in x:
            if xx not in [ref, '^', '!']:
                alternative_collection.append(xx)
    alternative_collection = list(set(alternative_collection))
    if len(alternative_collection) == 0:
        alt_i = None
        array = [0] * sample_size
        for x in t[6:]:
            s = x.split('|')
            sample_list = list_text2array(s[0])
            gt = s[1].split(':')
            if '^' in gt:
                for i in sample_list:
                    array[i] = 5
            elif '!' in gt:
                c = 0
                for gt_i in gt:
                    if gt_i != ref:
                        c += 1
                for i in sample_list:
                    array[i] = c
        array_collection.append(loci+[ref, alt_i] + array)
    else:
        for alt_i in alternative_collection:
            array = [0] * sample_size
            for x in t[6:]:
                s = x.split('|')
                sample_list = list_text2array(s[0])
                gt = s[1].split(':')
                if '^' in gt:
                    for i in sample_list:
                        array[i] = 5
                else:
                    c = 0
                    for gt_i in gt:
                        if gt_i == alt_i and gt_i != ref:
                            c += 1
                    for i in sample_list:
                        array[i] = c
            array_collection.append(loci+[ref, alt_i] + array)
    return array_collection


def convert_cdr2matrix(cdr_file, mask_dict={}):
    sample_index = {}
    for line in readlines_reverse(cdr_file):
        line = line.rstrip()
        if line.startswith("#") and 'FILE-INDEX' in line:
            t = line.split()
            sample_index[t[2]] = t[3]
        else:
            if not line.startswith("#"):
                break
    sample_size = len(sample_index)
    matrix = {}
    FH = open(cdr_file)
    variant_file = job_id+".variant.bed"
    OH = open(variant_file, 'w')
    for line in FH:
        line = line.rstrip()
        if line.startswith("#"):
            pass
        else:
            t = line.split()
            out = True
            if t[0] in mask_dict:
                if t[1] in mask_dict[t[0]]:
                    out = False
            if out:
                for variant in parse_cdr_line(line, sample_size):
                    loci = '-'.join(variant[0:3])
                    if loci not in matrix:
                        matrix[loci] = []
                    matrix[loci].append(variant)
                    if t[1] == t[2]:
                        OH.write(
                            '\t'.join(variant[0:3])+'\t'+'-'.join(variant[0:4])+'\n')
                    else:
                        OH.write('\t'.join(
                            [t[0], str(int(t[1])-1), str(int(t[2])+1)])+'\t'+'-'.join(variant[0:4])+'\n')
    FH.close()
    OH.close()
    return matrix, sample_index, variant_file


def collapse_variant(matrix, compress=False):
    def f(x):
        if x > 5:
            return 5
        return x
    for loci in matrix:
        temp = matrix[loci]
        new_temp = None
        if len(temp) == 1:
            new_temp = temp[0]
        else:
            new_temp = temp[0][0:3]
            _variant_type = list(set([x[3] for x in temp]))
            new_temp.append(','.join(_variant_type))
            new_temp.append(temp[0][4])
            _alt_allele = list(set([x[5] for x in temp]))
            new_temp.append(','.join(_alt_allele))
            new_array = numpy.array(temp[0][6:])
            for i in range(1, len(temp)):
                new_array += numpy.array(temp[i][6:])
            new_array = [f(x) for x in new_array]
            new_temp.extend(list(new_array))
        if compress:
            # zlib.decompress()
            _string = ''.join([str(ii) for ii in new_temp[6:]])
            new_temp = new_temp[0:6] + [zlib.compress(_string)]
        matrix[loci] = [new_temp]
    return matrix


def find_variant_to_print(variant_file, region_anno_file):
    temp_variant_to_print_file = job_id+'.variant_to_print.bed'
    cmd = ['bedtools', 'intersect', '-wa', '-wb', '-a', variant_file,
           '-b', region_anno_file, '>', temp_variant_to_print_file]
    cmd = ' '.join(cmd)
    returned_value = subprocess.call(cmd, shell=True)
    var_dict = {}
    FH = open(temp_variant_to_print_file)
    for line in FH:
        line = line.rstrip()
        t = line.split('\t')
        if len(t) > 9:
            gene = t[8]
            mRNA = t[7]
            var_type = t[9]
            loci = t[3]
            if gene not in var_dict:
                var_dict[gene] = {}
            if mRNA not in var_dict[gene]:
                var_dict[gene][mRNA] = {}
            var_dict[gene][mRNA][loci] = var_type
    return var_dict


def print_matrix_for_CARVA(variant_to_print, matrix, sample_index, output):
    OH = open(output, 'w')
    text = ['Gene', 'Transcript', 'Location'] + \
        [sample_index[str(i)] for i in range(len(sample_index))]
    OH.write('\t'.join(text)+'\n')

    for gene in variant_to_print:
        for mRNA in variant_to_print[gene]:
            for loci in variant_to_print[gene][mRNA]:
                var_type = variant_to_print[gene][mRNA][loci]
                new_loci = '-'.join(loci.split('-')[0:3])
                if new_loci in matrix:
                    for variant in matrix[new_loci]:
                        text = [str(_x) for _x in variant]
                        if var_type == 'SPLICE_A/D':
                            xxx = loci.split('-')
                            if xxx[3] != 'insertion' and xxx[3] != 'deletion':
                                text[3] = 'SPDA'
                            else:
                                text[3] = xxx[3]
                        OH.write(
                            '\t'.join(([gene, mRNA, new_loci+'-'+text[3]]+text[6:]))+'\n')
    OH.close()
    return None


# MAIN
if len(sys.argv) < 4:
    print "Need three para: input_cdr mask_file output_cdr."
    sys.exit()
job_id = 'TEMP00' + str(random.randint(10000000, 100000000-1))
input_cdr = sys.argv[1]
input_mask = sys.argv[2]
output_file = sys.argv[3]
mask_dict = load_mask_bed_file(input_mask)
gff3_file = '/rsrch3/home/epi/yyu4/database/Human/ref_UCSC_hg19/refGene_hg19.gff3'
region_anno_file = load_gff3_file(gff3_file)
matrix, sample_index, variant_file = convert_cdr2matrix(input_cdr, mask_dict)
variant_to_print = find_variant_to_print(variant_file, region_anno_file)
print_matrix_for_CARVA(variant_to_print, matrix, sample_index, output_file)
os.system('rm -f '+job_id+'.*')


# How to run
# python /rsrch3/home/epi/yyu4/local/XPAToolkit2/tools/cdr2matrix.py merge.small.cdr nomask out.matrix
# perl ~/local/bsub/bsub_m.pl cdr2matrx python /rsrch3/home/epi/yyu4/local/XPAToolkit2/tools/cdr2matrix.py merge.small.cdr nomask out.matrix
# bsub -J cdr2matrix -W 23:59 -n 6 -q medium -M 128 -R rusage[mem=128] -o cdr2matrix.log -e cdr2matrix.err "python ~/local/XPAToolkit2/tools/cdr2matrix.py merge.anno.edit.cdr nomask merge.anno.edit.matrix"

# Pending part
# matrix = collapse_variant(matrix, compress=False)
