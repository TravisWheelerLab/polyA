import sys
import re


def read_file(filename_hmm_align):
    """
    opens file and returns contents of file as one string
    """
    file_contents = ''
    with open(filename_hmm_align, 'r') as f_hmm_align:
        file_contents = f_hmm_align.read()

    return file_contents


def get_info(align_parts):
    """
    takes info line split by white space and assigns all info fields to variables
    """
    start: int
    stop: int
    subfam = ''
    consensus_start = ''
    consensus_stop = ''
    strand = '+'
    flank = ''
    chrom_seq = ''
    subfam_seq = ''
    strand_seq = '-1'

    model_len = align_parts[3].split()[13]
    align_parts = align_parts[6:]

    score_line = align_parts[0]
    score = score_line.split()[1]
    # on subfam
    subfam_line = align_parts[2]
    subfam_parts = subfam_line.split()
    subfam = subfam_parts[0]
    consensus_start = subfam_parts[1]
    subfam_seq = subfam_parts[2].upper().replace('.', '-')
    consensus_stop = subfam_parts[3]
    # name start subfam_seq stop
    # on chrom
    chrom_line = align_parts[4]
    chrom_parts = chrom_line.split()
    chrom = chrom_parts[0]
    start = chrom_parts[1]
    chrom_seq = chrom_parts[2].upper()
    stop = chrom_parts[3]
    # doesn't have chrom info because it's an artificial seq
    if 'chr' not in chrom:
        chrom = 'chr0:0000-0000'

    flank = str(int(model_len) - int(consensus_stop))
    if int(stop) < int(start):
        strand = '-'
        strand_seq = 't'
        start, stop = stop, start
    return subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, chrom_seq, subfam_seq, strand_seq


def print_info(subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, strand_seq, chrom_seq, subfam_seq, f_out):
    """
    prints all info in table format
    """
    f_out.write(f'{subfam} {chrom} {score} {strand} {strand_seq} {start} {stop} {consensus_start} {consensus_stop} {flank} {chrom_seq} {subfam_seq}\n')


if __name__ == "__main__":
    filename_hmm_align = sys.argv[1]
    file_contents = read_file(filename_hmm_align)

    filename_out = filename_hmm_align + ".tbl.txt"
    f_out = open(filename_out, 'w')
    f_out.write("# ID TR SC SD TG ST SP CST CSP FL chrom subfam\n")

    # =GF ID  L1MA3_3end
    # =GF TR  chr2:128345832-128354866
    # =GF SC  455.1
    # =GF SD  -
    # =GF TQ  t
    # =GF ST  3262
    # =GF SP  3714
    # =GF CST 587
    # =GF CSP 1049
    # =GF FL  3

    # alignments are between >> and PP chars in hmm alignment file
    alignments = re.findall(r'>>([\S\s]*?)PP', file_contents)
    for region in alignments:
        align_parts = region.splitlines()
        subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, chrom_seq, subfam_seq, strand_seq = get_info(align_parts)
        print_info(subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, strand_seq, chrom_seq, subfam_seq, f_out)
    f_out.close()