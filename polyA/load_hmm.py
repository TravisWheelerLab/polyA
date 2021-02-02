from typing import TextIO, List, Dict


# TODO(george, audrey): Figure out the true return type
def load_hmm(file: List[str], subfams: List[str]) -> Dict[str, dict]:
    hmm_dict = {}
    char_pos: List[str] = []
    i = 0
    while i < len(file):
        line = file[i]
        # new subfam HMM
        if line.strip().startswith("NAME"):
            parts = line.strip().split()
            name = parts[1]
            if name in subfams:
                hmm_sub_dict = {}
                while line and not line.strip().startswith("//"):
                    if line.strip().startswith("HMM"):
                        chars = line.strip().split()[1:]
                        i += 1
                        line = file[i]
                        transitions = line.strip().split()
                    if line.strip().startswith("COMPO"):
                        # first row of COMPO
                        line = file[i]
                        unrelated_vals = line.strip().split()[1::]
                        unrelated_dict = {}
                        for j in range(len(chars)):
                            unrelated_dict[chars[j]] = unrelated_vals[j]
                        hmm_sub_dict["random"] = unrelated_dict
                    if line.strip().split()[0].isdigit():
                        # get per position scores
                        # read in next 3 lines
                        pos_dict = {}
                        emmission_dict = {}
                        transition_dict = {}
                        pos = line.strip().split()[0]
                        emmission = line.strip().split()[1:]
                        for j in range(len(chars)):
                            emmission_dict[chars[j]] = emmission[j]
                        i += 2  # 2 read lines
                        line = file[i]
                        transition = line.strip().split()
                        for j in range(len(transition)):
                            transition_dict[transitions[j]] = transition[j]
                        pos_dict["emission"] = emmission_dict
                        pos_dict["transition"] = transition_dict
                        hmm_sub_dict[int(pos)] = pos_dict
                    i += 1
                    line = file[i]
                hmm_dict[name] = hmm_sub_dict
        i += 1
    return hmm_dict
