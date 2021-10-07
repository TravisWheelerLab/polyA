from typing import Tuple
from subprocess import Popen, PIPE


def pick_merged_subfam_name(subfam_pair: Tuple[str, str]):
    # input: ('AluSg4#SINE/Alu', 'AluSx1#SINE/Alu')
    # output: AluSg4, AluSx1, AluSg4_merged
    subfam_A = subfam_pair[0].split("#")[0]
    subfam_B = subfam_pair[1].split("#")[0]
    # pick a merged subfam name
    merged_name = ""
    if subfam_A.endswith("_merged"):
        # subfam A has been merged
        merged_name = subfam_A
    elif subfam_B.endswith("_merged"):
        # subfam B has been merged
        merged_name = subfam_B
        # switch names so subfam A is always the merged one
        subfam_A, subfam_B = subfam_B, subfam_A
    else:
        # neither have been merged yet
        # default merged name on subfam A
        merged_name = subfam_A + "_merged"
    return subfam_A, subfam_B, merged_name


def merge_subfams(subfam_pair: Tuple[str, str]):
    msa_dir_path = "/Users/audrey/Desktop/notebook/2021/subfamily-reliability/alu-data/alu-subfams-msa/"
    subfam_A, subfam_B, merged_name = pick_merged_subfam_name(subfam_pair)
    subfam_A_instances = msa_dir_path + subfam_A + ".fas"
    subfam_B_instances = msa_dir_path + subfam_B + ".fas"
    merged_subfam_instances = msa_dir_path + merged_name + ".fas"

    # create file with all instances
    outfile = open(merged_subfam_instances, "w")
    with open(subfam_A_instances, "r") as infile:
        outfile.write(infile.read())
    with open(subfam_B_instances, "r") as infile:
        outfile.write(infile.read())
    outfile.close()

    # create MSA from the merged instances file
    # TODO: eventually use refiner for this
    merged_subfam_msa = msa_dir_path + merged_name + ".afa"
    with open(merged_subfam_msa, "w") as f_out_merged_msa:
        process = Popen(['mafft', merged_subfam_instances], stdout=f_out_merged_msa, stderr=PIPE)
        process.communicate()
    # should now have the msa for the new thing
    # get a consensus from this?
    hmm_out = msa_dir_path + "merged_subfam.hmm"
    process = Popen(['hmmbuild', hmm_out, merged_subfam_msa], stdout=PIPE, stderr=PIPE)
    process.communicate()
    process = Popen(['hmmemit', '-c', hmm_out], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    merged_consensus_seq = str(stdout).split("consensus")[1].replace("\\n", "\n")
    merged_consensus_seq = merged_consensus_seq[:len(merged_consensus_seq) - 1]
    return merged_consensus_seq





