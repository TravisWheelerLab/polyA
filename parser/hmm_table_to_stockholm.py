# read in table files -> stockholm format
# basically just print it differently
import sys
import re


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
    filename_hmm_alignment = sys.argv[1]
    filename_out = filename_hmm_alignment + ".sto"
    f_out = open(filename_out, 'w')
    f_out.write("# STOCKHOLM 1.0\n")

    # skip first line with table header
    file_content = open(filename_hmm_alignment).readlines()[1:]
    for alignment in file_content:
        subfam, chrom, score, strand, strand_seq, start, stop, consensus_start, consensus_stop, flank, chrom_seq, subfam_seq = alignment.split()
        print_info(subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, strand_seq, f_out)
        print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out)

    f_out.close()
