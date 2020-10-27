import sys
import re
from sys import stdout


def read_file(filename_cm):
    """
    opens file and returns contents of file as one string
    """
    file_contents = ''
    with open(filename_cm, 'r') as f_cm:
        file_contents = f_cm.read()

    return file_contents



def get_info(info_array):
    """
    takes info line split by white space and assigns all info fields to variables
    """
    score = info_array[0]
    chrom = info_array[4]
    chrom = info_array[4]

    # doesn't have chrom info because it's an artificial seq
    if 'chr' not in chrom:
        chrom = 'chr0:0000-0000'

    start = info_array[5]
    stop = info_array[6]
    subfam = ''
    consensus_start = ''
    consensus_stop = ''
    strand = ''
    flank = ''

    # if 6 contains (), then start/stop = 7/8
    if "(" in info_array[6]:  # FIXME - shouldn't this be 5?
        start = info_array[6]
        stop = info_array[7]

    if len(info_array) == 12:  # not compliment
        strand = '+'
        subfam = info_array[8]
        consensus_start = info_array[9]
        consensus_stop = info_array[10]
        flank = info_array[11]

        # if 10 contains (), then start/stop = 7/8
        if "(" in info_array[9]:
            flank = info_array[9]
            start = info_array[10]
            stop = info_array[11]

    else:  # compliment
        strand = '-'
        subfam = info_array[9]
        consensus_start = info_array[10]
        consensus_stop = info_array[11]
        flank = info_array[12]

        if "(" in info_array[10]:
            flank = info_array[10]
            consensus_start = info_array[11]
            consensus_stop = info_array[12]

    flank = flank.strip(")")
    flank = flank.strip("(")

    return subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank


def print_info(C, subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, f_out):
    """
    prints all info in correct format for stockholm
    """
    f_out.write(f'#=GF ID  {subfam}\n')
    f_out.write(f'#=GF TR  {chrom}\n')
    f_out.write(f'#=GF SC  {score}\n')
    f_out.write(f'#=GF ST  {strand}\n')

    if strand == '-':
        if C:
            f_out.write(f'#=GF TQ  t\n')
        else:
            f_out.write(f'#=GF TQ  q\n')
    else:
        f_out.write(f'#=GF TQ  -1\n')

    f_out.write(f'#=GF ST  {start}\n')
    f_out.write(f'#=GF SP  {stop}\n')
    f_out.write(f'#=GF CST {consensus_start}\n')
    f_out.write(f'#=GF CSP {consensus_stop}\n')
    f_out.write(f'#=GF FL  {flank}\n')


def get_alignment(alignment_array):
    """
    takes cm alignment split on newline and return string of chrom seq and subfam seq
    """
    chrom_seq = ""
    subfam_seq = ""

    # for (my $i = 0; $i < @ Alignment; $i = $i + 4){
    for i in range(0, len(alignment_array), 4):
        m_chrom = re.search(r'.+?\s\d+\s(.+?)\s\d+', alignment_array[i])
        if m_chrom:
            chrom_seq += m_chrom.group(1)

        m_subfam = re.search(r'.+?\s\d+\s(.+?)\s\d+', alignment_array[i+2])
        if m_chrom:
            subfam_seq += m_subfam.group(1)

    return chrom_seq, subfam_seq


def print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out):
    """
    takes chrom and subfam names and chrom and subfam seqs and prints in
    stockholm format
    """
    f_out.write(f"{chrom}   {chrom_seq}\n")
    f_out.write(f"{subfam}   {subfam_seq}\n")

    f_out.write("//\n")


if __name__ == "__main__":
    filename_cm = sys.argv[1]
    file_contents = read_file(filename_cm)

    filename_out = filename_cm + ".sto"
    f_out = open(filename_out, 'w')

    f_out.write("# STOCKHOLM 1.0\n")

    alignments = re.findall(r'\s*?\d+\s+[0-9]+\.[0-9]+\s+[0-9.]+\s+[0-9.]+\s+.+?\n\n[\s\S]+?Transitions', file_contents)

    for region in alignments:
        info_line: str = ''
        alignment: str = ''
        region = region.strip()
        m = re.match(r'(.+?)\n([\s\S]+)\n\nTransitions', region)
        if m:
            info_line = m.group(1)
            alignment = m.group(2)

        alignment = alignment.strip("\n")

        info_array = info_line.split()

        subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank = get_info(info_array)

        alignment_array = alignment.split("\n")

        # if C is true: on reverse strand and target is reversed
        C = 0
        if "C" in alignment_array[0]:
            C = 1

        chrom_seq, subfam_seq = get_alignment(alignment_array)

        print_info(C, subfam, chrom, score, strand, start, stop, consensus_start, consensus_stop, flank, f_out)
        print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out)

    f_out.close()