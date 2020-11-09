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
    strand = '+' # how to get this
    flank = '' # ? How much is left in subfam
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

    if int(consensus_stop) < int(consensus_start): # strings
        strand = '-'
        strand_seq = 't'
        flank = str(int(model_len) - int(consensus_start))
    else:
        flank = str(int(model_len) - int(consensus_stop))
    if int(stop) < int(start):
        strand = '-'
        strand_seq = 'q'
    return subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, chrom_seq, subfam_seq, strand_seq


def print_info(subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, strand_seq, f_out):
    """
    prints all info in correct format for stockholm
    """
    f_out.write(f'#=GF ID  {subfam}\n')
    f_out.write(f'#=GF TR  {chrom}\n')
    f_out.write(f'#=GF SC  {score}\n')
    f_out.write(f'#=GF SD  {strand}\n')
    f_out.write(f'#=GF TQ  {strand_seq}\n')
    f_out.write(f'#=GF ST  {start}\n')
    f_out.write(f'#=GF SP  {stop}\n')
    f_out.write(f'#=GF CST {consensus_start}\n')
    f_out.write(f'#=GF CSP {consensus_stop}\n')
    f_out.write(f'#=GF FL  {flank}\n')


def print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out):
    """
    takes chrom and subfam names and chrom and subfam seqs and prints in
    stockholm format
    """
    f_out.write(f"{chrom}   {chrom_seq}\n")
    f_out.write(f"{subfam}   {subfam_seq}\n")

    f_out.write("//\n")


if __name__ == "__main__":
    filename_hmm_align = sys.argv[1]
    file_contents = read_file(filename_hmm_align)

    filename_out = filename_hmm_align + ".sto"
    f_out = open(filename_out, 'w')
    f_out.write("# STOCKHOLM 1.0\n")

    alignments = re.findall(r'>>([\S\s]*?)PP', file_contents)
    for region in alignments:
        align_parts = region.splitlines()
        subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, chrom_seq, subfam_seq, strand_seq = get_info(align_parts)
        print_info(subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, strand_seq, f_out)
        print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out)
    f_out.close()