from typing import TextIO, List, Dict


# TODO(george, audrey): Figure out the true return type
def load_hmm(file: List[str], subfams: List[str]) -> Dict[str, dict]:
    hmm_dict = {}
    char_pos: List[str] = []
    i = 0
    line = file[i]
    while i < len(file):
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
                        i += 1
                        line = file[i]
                        unrelated_vals = line.strip().split()
                        unrelated_dict = {}
                        for i in range(len(chars)):
                            unrelated_dict[chars[i]] = unrelated_vals[i]
                        hmm_sub_dict["unrelated"] = unrelated_dict
                    if line.strip().split()[0].isdigit():
                        # get per position scores
                        # read in next 3 lines
                        pos_dict = {}
                        emmission_dict = {}
                        transition_dict = {}
                        pos = line.strip().split()[0]
                        emmission = line.strip().split()[1:]
                        for i in range(len(chars)):
                            emmission_dict[chars[i]] = emmission[i]
                        i += 1
                        line = file[i]
                        transition = file.readline().strip().split()
                        for i in range(len(transition)):
                            transition_dict[transitions[i]] = transition[i]
                        pos_dict["emission"] = emmission_dict
                        pos_dict["transition"] = transition_dict
                        hmm_sub_dict[int(pos)] = pos_dict
                    i += 1
                    line = file[i]
                hmm_dict[name] = hmm_sub_dict
        i += 1
        line = file[i]
    return hmm_dict
