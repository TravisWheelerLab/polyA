import re


def read_file(filename_cm):
    """
    opens file and returns contents of file as one string
    """
    file_contents = ""
    with open(filename_cm, "r") as f_cm:
        file_contents = f_cm.read()

    return file_contents


def get_info(info_array):
    """
    takes info line split by white space and assigns all info fields to variables
    """
    score = info_array[0]
    chrom = info_array[4]

    # format chrom name to match "(.+):(\d+)-(\d+)"
    if "chr" not in chrom.lower():
        # doesn't have chrom info because it's an artificial seq
        chrom = "chr0:0000-0000"
    chrom = chrom.replace("Chr", "chr")
    chrom_start, chrom_end = chrom.split("chr")
    chrom_end_values = chrom_end.split("_")
    if len(chrom_end_values) == 3:
        # reformat chrom_end
        chrom_end = (
            chrom_end_values[0]
            + ":"
            + chrom_end_values[1]
            + "-"
            + chrom_end_values[2]
        )
    chrom = chrom_start + "chr" + chrom_end

    start = info_array[5]
    stop = info_array[6]
    subfam = ""
    consensus_start = ""
    consensus_stop = ""
    strand = ""
    flank = ""

    # if 6 contains (), then start/stop = 7/8
    if "(" in info_array[6]:  # FIXME - shouldn't this be 5?
        start = info_array[6]
        stop = info_array[7]

    # 14 for cm .align files, might need to change
    if len(info_array) == 14:  # not compliment
        strand = "+"
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
        strand = "-"
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

    return (
        subfam,
        chrom,
        score,
        strand,
        start,
        stop,
        consensus_start,
        consensus_stop,
        flank,
    )


def get_score_matrix(matrix_name):
    """
    grabs score matrix info from original matrix file, puts it in correct format, returns string
    """
    matrix_file_contents = read_file("../fixtures/matrices/" + matrix_name)
    matrix_file_contents = matrix_file_contents.split("\n")
    background_freqs_line = matrix_file_contents[0]
    # reformat background freqs
    background_freqs = background_freqs_line.split()[1:]
    background_freqs_dict = {}
    for i in range(0, len(background_freqs), 2):
        char = background_freqs[i][0]
        freq = float(background_freqs[i + 1])
        background_freqs_dict[char] = freq
    freqs_string = "BACKGROUND FREQS: " + str(background_freqs_dict)
    return "\n".join(([freqs_string] + matrix_file_contents[1:]))


def print_info(
    C,
    subfam,
    chrom,
    score,
    strand,
    start,
    stop,
    consensus_start,
    consensus_stop,
    flank,
    matrix_name,
    f_out_sto,
):
    """
    prints all info in correct format for stockholm
    """
    f_out_sto.write(f"#=GF ID  {subfam}\n")
    f_out_sto.write(f"#=GF TR  {chrom}\n")
    f_out_sto.write(f"#=GF SC  {score}\n")
    f_out_sto.write(f"#=GF SD  {strand}\n")

    if strand == "-":
        if C:
            f_out_sto.write(f"#=GF TQ  t\n")
        else:
            f_out_sto.write(f"#=GF TQ  q\n")
    else:
        f_out_sto.write(f"#=GF TQ  -1\n")

    f_out_sto.write(f"#=GF ST  {start}\n")
    f_out_sto.write(f"#=GF SP  {stop}\n")
    f_out_sto.write(f"#=GF CST {consensus_start}\n")
    f_out_sto.write(f"#=GF CSP {consensus_stop}\n")
    f_out_sto.write(f"#=GF FL  {flank}\n")
    f_out_sto.write(f"#=GF MX  {matrix_name}\n")


def get_alignment(alignment_array):
    """
    takes cm alignment split on newline and return string of chrom seq and subfam seq
    """
    chrom_seq = ""
    subfam_seq = ""

    # for (my $i = 0; $i < @ Alignment; $i = $i + 4){
    for i in range(0, len(alignment_array), 4):
        m_chrom = re.search(r".+?\s\d+\s(.+?)\s\d+", alignment_array[i])
        if m_chrom:
            chrom_seq += m_chrom.group(1)

        m_subfam = re.search(r".+?\s\d+\s(.+?)\s\d+", alignment_array[i + 2])
        if m_chrom:
            subfam_seq += m_subfam.group(1)

    return chrom_seq, subfam_seq


def print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out_sto):
    """
    takes chrom and subfam names and chrom and subfam seqs and prints in
    stockholm format
    """
    f_out_sto.write(f"{chrom}   {chrom_seq}\n")
    f_out_sto.write(f"{subfam}   {subfam_seq}\n")

    f_out_sto.write("//\n")


def print_score_matrix(f_out_matrix, score_matrix, matrix_name):
    """
    print score matrix to its output file with extension ".matrix"
    """
    f_out_matrix.write(matrix_name)
    f_out_matrix.write("\n")
    f_out_matrix.write(score_matrix)
    f_out_matrix.write("//\n")


def convert(filename_rm: str, filename_out_sto: str, filename_out_matrix: str):
    matrices = {}
    file_contents = read_file(filename_rm)

    if not filename_out_sto:
        filename_out_sto = filename_rm + ".sto"
    f_out_sto = open(filename_out_sto, "w")

    if not filename_out_matrix:
        filename_out_matrix = filename_rm + ".matrix"
    f_out_matrix = open(filename_out_matrix, "w")

    f_out_sto.write("# STOCKHOLM 1.0\n")
    f_out_sto.write(f"# ALIGNMENT TOOL RepeatMasker\n")

    alignments = re.findall(
        r"\s*?\d+\s+[0-9]+\.[0-9]+\s+[0-9.]+\s+[0-9.]+\s+.+?\n\n[\s\S]+?Transitions",
        file_contents,
    )

    for region in alignments:
        info_line: str = ""
        alignment: str = ""
        matrix_name: str = ""
        region = region.strip()

        m_matrix = re.search(r"Matrix = (.+?)\n", region)
        if m_matrix is None:
            raise RuntimeError("did not find Matrix stanza")
        matrix_name = m_matrix[1]

        if (
            matrix_name not in matrices
        ):  # do not put duplicate of matrices in output file
            matrices[matrix_name] = 0
            score_matrix = get_score_matrix(matrix_name)
            print_score_matrix(f_out_matrix, score_matrix, matrix_name)

        m = re.match(r"(.+?)\n([\s\S]+)\n\nMatrix", region)
        if m:
            info_line = m.group(1)
            alignment = m.group(2)

        alignment = alignment.strip("\n")

        info_array = info_line.split()

        (
            subfam,
            chrom,
            score,
            strand,
            start,
            stop,
            consensus_start,
            consensus_stop,
            flank,
        ) = get_info(info_array)

        alignment_array = alignment.split("\n")

        # if C is true: on reverse strand and target is reversed
        C = 0
        if "C" in alignment_array[0]:
            C = 1

        chrom_seq, subfam_seq = get_alignment(alignment_array)

        print_info(
            C,
            subfam,
            chrom,
            score,
            strand,
            start,
            stop,
            consensus_start,
            consensus_stop,
            flank,
            matrix_name,
            f_out_sto,
        )
        print_alignment(chrom_seq, subfam_seq, chrom, subfam, f_out_sto)

    f_out_sto.close()
    f_out_matrix.close()
