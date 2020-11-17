from typing import TextIO, List, Dict


def load_hmm(file: TextIO, subfams: List[str]):
    hmm_dict = {}
    char_pos: List[str] = []
    line = file.readline()
    while line:
        # new subfam HMM
        if line.strip().startswith("NAME"):
            parts = line.strip().split()
            name = parts[1]
            if name in subfams:
                hmm_sub_dict = {}
                while line and not line.strip().startswith("//"):
                    if line.strip().startswith("HMM"):
                        char_pos = line.strip().split()[1:]
                        line = file.readline()
                        transition_pos = line.strip().split()
                    if line.strip().startswith("COMPO"):
                        line = file.readline()
                        unrelated_vals = line.strip().split()
                        unrelated_dict = {}
                        for i in range(len(char_pos)):
                            unrelated_dict[char_pos[i]] = unrelated_vals[i]
                        hmm_sub_dict["unrelated"] = unrelated_dict
                    if line.strip().split()[0].isdigit():
                        # get per position scores
                        # read in next 3 lines
                        pos_dict = {}
                        emmission_dict = {}
                        transition_dict = {}
                        pos = line.strip().split()[0]
                        emmission = line.strip().split()[1:]
                        for i in range(len(char_pos)):
                            emmission_dict[char_pos[i]] = emmission[i]
                        file.readline()
                        transition = file.readline().strip().split()
                        for i in range(len(char_pos)):
                            transition_dict[transition_pos[i]] = transition[i]
                        pos_dict["emission"] = emmission_dict
                        pos_dict["transition"] = transition_dict
                        hmm_sub_dict[pos] = pos_dict
                    line = file.readline()
                hmm_dict[name] = hmm_sub_dict
        line = file.readline()
    return hmm_dict
