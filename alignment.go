package polyA

var skipState = &Alignment{
	Subfamily:         "skip",
	ChromName:         "skip_chrom",
	ChromStart:        0,
	ChromStop:         0,
	Score:             0,
	Start:             0,
	Stop:              0,
	ConsensusStart:    0,
	ConsensusStop:     0,
	Sequence:          "",
	SubfamilySequence: "",
	Strand:            "",
	Flank:             0,
	SubMatrixName:     "",
	GapInit:           0,
	GapExt:            0,
}

func SkipState() *Alignment {
	return skipState
}

type Alignment struct {
	Subfamily         string
	ChromName         string
	ChromStart        int
	ChromStop         int
	Score             int
	Start             int
	Stop              int
	ConsensusStart    int
	ConsensusStop     int
	Sequence          string
	SubfamilySequence string
	Strand            string
	Flank             int
	SubMatrixName     string
	GapInit           float64
	GapExt            float64
}

func (a *Alignment) ChromLength() int {
	return a.ChromStop - a.ChromStart + 1
}
