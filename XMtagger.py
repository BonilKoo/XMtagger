import sys
import os
import pysam

def mkdir(dir_name):
    os.system(f'mkdir -p {dir_name}')

def check_index_file(bamfile):
    if not os.path.exists(bamfile+'.bai'):
        print('\nError: The bam file is not indexed.')
        sys.exit()

def check_XM_tag(bamfile):
    for read in bamfile.fetch():
        if read.has_tag('XM'):
            if read.get_tag('XM', with_value_type=True)[1] == 'Z':
                return True
        else:
            return False

def check_paired_end(bamfile):
    for read in bamfile.fetch():
        if read.is_paired:
            return True
        else:
            return False

def check_need_reverse_complement(read):
    if (not read.is_reverse and read.is_read1) or (read.is_reverse and read.is_read2):
        return False
    else:
        return True

def reverse_complement(seq):
    return seq.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()[::-1]

bam_file_name = sys.argv[2]
check_index_file(bam_file_name)
bam_file = pysam.AlignmentFile(bam_file_name)
if check_XM_tag(bam_file):
    print('\nError: XM tag already exists.')
    sys.exit()

ref_file_name = sys.argv[1]
new_bam_file_name = sys.argv[3]
mkdir('/'.join(new_bam_file_name.split('/')[:-1]))

ref_genome = pysam.FastaFile(ref_file_name)
new_bam_file = pysam.AlignmentFile(new_bam_file_name, 'wb', template=bam_file)

flag_paired_end = check_paired_end(bam_file)

print('\nStart making XM tag.\n')
read_num = 0
for read in bam_file.fetch():
    read_num += 1
    read_seq = read.query_alignment_sequence.upper()
    ref_seq = ref_genome.fetch(read.reference_name, read.reference_start-2, read.reference_end+2).upper()

    tmp_read_seq = '--'
    tmp_ref_seq = ref_seq[:2]
    used_read_len = 0
    used_ref_len = 0

    for operation, length in read.cigartuples:
        if operation == 0:
            tmp_read_seq += read_seq[used_read_len:used_read_len+length]
            tmp_ref_seq += ref_seq[2+used_ref_len:2+used_ref_len+length]
            used_read_len += length
            used_ref_len += length
        elif operation == 1:
            tmp_read_seq += read_seq[used_read_len:used_read_len+length]
            tmp_ref_seq += '-' * length
            used_read_len += length
        elif operation == 2:
            tmp_read_seq += '-' * length
            tmp_ref_seq += ref_seq[2+used_ref_len:2+used_ref_len+length]
            used_ref_len += length
    tmp_read_seq += '--'
    tmp_ref_seq += ref_seq[-2:]

    if flag_paired_end:
        flag_reverse_complement = check_need_reverse_complement(read)
    else:
        flag_reverse_complement = read.is_reverse

    if flag_reverse_complement:
        tmp_read_seq = reverse_complement(tmp_read_seq[:-2])
        tmp_ref_seq = reverse_complement(tmp_ref_seq[:-2])
    else:
        tmp_read_seq = tmp_read_seq[2:]
        tmp_ref_seq = tmp_ref_seq[2:]

    XM_tag = ''
    for idx in range(len(tmp_read_seq)-2):
        if tmp_read_seq[idx] == '-':
            continue
        elif tmp_read_seq[idx] == 'N':
            XM_tag += '.'
        elif tmp_ref_seq[idx] == 'C':
            if (tmp_read_seq[idx+1] == '-' or tmp_read_seq[idx+2] == '-') and (idx not in [len(tmp_read_seq)-3, len(tmp_read_seq)-4]):
                tmp_tmp_read_seq = tmp_read_seq[idx]
                tmp_tmp_ref_seq = tmp_ref_seq[idx]
                flag_tmp = 0
                tmp_count = 1
                while flag_tmp != 2:
                    if tmp_read_seq[idx+tmp_count] != '-':
                        tmp_tmp_read_seq += tmp_read_seq[idx+tmp_count]
                        tmp_tmp_ref_seq += tmp_ref_seq[idx+tmp_count]
                        flag_tmp += 1
                    tmp_count += 1
                if tmp_tmp_ref_seq[:2] == 'CG':
                    if tmp_tmp_read_seq[0] == 'C':
                        XM_tag += 'Z'
                    elif tmp_tmp_read_seq[0] == 'T':
                        XM_tag += 'z'
                    else:
                        XM_tag += '.'
                elif tmp_tmp_ref_seq in ['CAG', 'CTG', 'CCG']:
                    if tmp_tmp_read_seq[0] == 'C':
                        XM_tag += 'X'
                    elif tmp_tmp_read_seq[0] == 'T':
                        XM_tag += 'x'
                    else:
                        XM_tag += '.'
                elif tmp_tmp_ref_seq in ['CAA', 'CAT', 'CAC', 'CTA', 'CTT', 'CTC', 'CCA', 'CCT', 'CCC']:
                    if tmp_tmp_read_seq[0] == 'C':
                        XM_tag += 'H'
                    elif tmp_tmp_read_seq[0] == 'T':
                        XM_tag += 'h'
                    else:
                        XM_tag += '.'
                elif '-' in tmp_tmp_ref_seq or 'N' in tmp_tmp_ref_seq:
                    if tmp_tmp_read_seq[0] == 'C':
                        XM_tag += 'U'
                    elif tmp_tmp_read_seq[0] == 'T':
                        XM_tag += 'u'
                    else:
                        XM_tag += '.'

            elif tmp_ref_seq[idx:idx+2] == 'CG':
                if tmp_read_seq[idx] == 'C':
                    XM_tag += 'Z'
                elif tmp_read_seq[idx] == 'T':
                    XM_tag += 'z'
                else:
                    XM_tag += '.'
            elif tmp_ref_seq[idx:idx+3] in ['CAG', 'CTG', 'CCG']:
                if tmp_read_seq[idx] == 'C':
                    XM_tag += 'X'
                elif tmp_read_seq[idx] == 'T':
                    XM_tag += 'x'
                else:
                    XM_tag += '.'
            elif tmp_ref_seq[idx:idx+3] in ['CAA', 'CAT', 'CAC', 'CTA', 'CTT', 'CTC', 'CCA', 'CCT', 'CCC']:
                if tmp_read_seq[idx] == 'C':
                    XM_tag += 'H'
                elif tmp_read_seq[idx] == 'T':
                    XM_tag += 'h'
                else:
                    XM_tag += '.'
            elif '-' in tmp_ref_seq[idx:idx+3] or 'N' in tmp_ref_seq[idx:idx+3]:
                if tmp_read_seq[idx] == 'C':
                    XM_tag += 'U'
                elif tmp_read_seq[idx] == 'T':
                    XM_tag += 'u'
                else:
                    XM_tag += '.'
        else:
            XM_tag += '.'
    if flag_reverse_complement:
        XM_tag = XM_tag[::-1]

    read.set_tag(tag='XM', value=XM_tag, value_type='Z')
    new_bam_file.write(read)
    if read_num % 100000 == 0:
        print(f'{read_num} reads are tagged.')
ref_genome.close()
bam_file.close()
new_bam_file.close()
print('\nXM tagging completed!!!\n')
