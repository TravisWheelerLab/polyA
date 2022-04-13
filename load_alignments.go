package polyA

import (
	"bufio"
	"io"
	"regexp"
	"strconv"
	"strings"
)

var chromMetaRegexp *regexp.Regexp

func init() {
	chromMetaRegexp = regexp.MustCompile("(.+):(\\d+)-(\\d+)")
}

func parseMetaLine(line string) (key string, value string, found bool) {
	trimmedLine := strings.Trim(line, " \t")
	if strings.HasPrefix(trimmedLine, "#=") {
		parts := strings.Fields(trimmedLine)
		if len(parts) != 3 {
			return "", "", false
		}

		return parts[1], parts[2], true
	}

	return "", "", false
}

func parseAlignmentLine(line string) (key string, value string, found bool) {
	trimmedLine := strings.Trim(line, " \t")
	parts := strings.Fields(trimmedLine)
	if len(parts) != 2 {
		return "", "", false
	}

	return parts[0], parts[1], true
}

func parsePreambleLine(line string) bool {
	trimmedLine := strings.Trim(line, " \t")
	upperLine := strings.ToUpper(trimmedLine)
	return strings.HasPrefix(upperLine, "# STOCKHOLM")
}

func parseToolLine(line string) string {
	trimmedLine := strings.Trim(line, " \t")
	upperLine := strings.ToUpper(trimmedLine)
	if strings.HasPrefix(upperLine, "# ALIGNMENT TOOL") {
		parts := strings.Fields(trimmedLine)
		return strings.Join(parts[3:], " ")
	}

	return ""
}

func parseTerminatorLine(line string) bool {
	return strings.Trim(line, " \t") == "//"
}

func parseChromMeta(line string) (name string, start int, stop int, found bool) {
	trimmedLine := strings.Trim(line, " \t")
	match := chromMetaRegexp.FindStringSubmatch(trimmedLine)
	if match == nil {
		return "", 0, 0, false
	}

	chromName := match[1]
	chromStart, err := strconv.Atoi(match[2])
	if err != nil {
		panic("bad match - failed to convert to int")
	}
	chromStop, err := strconv.Atoi(match[3])
	if err != nil {
		panic("bad match - failed to convert to int")
	}

	return chromName, chromStart, chromStop, true
}

type Loader struct {
	addSkipState bool
	done         bool
	loadedCount  int
	scanner      *bufio.Scanner
}

func NewLoader(file io.Reader, addSkipState bool) *Loader {
	return &Loader{
		addSkipState: addSkipState,
		scanner:      bufio.NewScanner(file),
	}
}

func (l *Loader) Load(a *Alignment) bool {
	if l.done {
		return false
	}

	if l.loadedCount == 0 && l.addSkipState {
		a.Subfamily = skipState.Subfamily
		a.ChromName = skipState.ChromName
		a.ChromStart = skipState.ChromStart
		a.ChromStop = skipState.ChromStop
		a.Score = skipState.Score
		a.Start = skipState.Start
		a.Stop = skipState.Stop
		a.ConsensusStart = skipState.ConsensusStart
		a.ConsensusStop = skipState.ConsensusStop
		a.Sequence = skipState.Sequence
		a.SubfamilySequence = skipState.SubfamilySequence
		a.Strand = skipState.Strand
		a.Flank = skipState.Flank
		a.SubMatrixName = skipState.SubMatrixName
		a.GapInit = skipState.GapInit
		a.GapExt = skipState.GapExt

		l.loadedCount += 1
		return true
	}

	meta := make(map[string]string)
	var seqs []string

	for {
		more := l.scanner.Scan()
		line := l.scanner.Text()

		if parsePreambleLine(line) {
			continue
		}

		if parseToolLine(line) != "" {
			continue
		}

		key, value, found := parseMetaLine(line)
		if found {
			meta[key] = value
			continue
		}

		_, seq, found := parseAlignmentLine(line)
		if found {
			seqs = append(seqs, seq)
			continue
		}

		if parseTerminatorLine(line) {
			if len(seqs) != 2 {
				panic("file format error - wrong number of sequences")
			}

			chromName, chromStart, chromStop, found := parseChromMeta(meta["TR"])
			if !found {
				panic("file format error - missing TR field")
			}

			// TODO: Write helper(s) to check for presence of key and then get the value

			score, err := strconv.Atoi(meta["SC"])
			if err != nil {
				panic("file format error - invalid SC field")
			}

			start, err := strconv.Atoi(meta["ST"])
			if err != nil {
				panic("file format error - invalid ST field")
			}

			stop, err := strconv.Atoi(meta["SP"])
			if err != nil {
				panic("file format error - invalid SP field")
			}

			consensusStart, err := strconv.Atoi(meta["CST"])
			if err != nil {
				panic("file format error - invalid CST field")
			}

			consensusStop, err := strconv.Atoi(meta["CSP"])
			if err != nil {
				panic("file format error - invalid CSP field")
			}

			var sequence string
			var subfamilySequence string
			if meta["TQ"] == "t" {
				sequence = reverse(seqs[0])
				subfamilySequence = reverse(seqs[1])
			} else {
				sequence = seqs[0]
				subfamilySequence = seqs[1]
			}

			flank, err := strconv.Atoi(meta["FL"])
			if err != nil {
				panic("file format error - invalid FL field")
			}

			gapInit, err := strconv.ParseFloat(meta["GI"], 64)
			if err != nil {
				panic("file format error - invalid GI field")
			}

			gapExt, err := strconv.ParseFloat(meta["GE"], 64)
			if err != nil {
				panic("file format error - invalid GE field")
			}

			a.Subfamily = meta["ID"]
			a.ChromName = chromName
			a.ChromStart = chromStart
			a.ChromStop = chromStop
			a.Score = score
			a.Start = start
			a.Stop = stop
			a.ConsensusStart = consensusStart
			a.ConsensusStop = consensusStop
			a.Sequence = sequence
			a.SubfamilySequence = subfamilySequence
			a.Strand = meta["SD"]
			a.Flank = flank
			a.SubMatrixName = meta["MX"]
			a.GapInit = gapInit
			a.GapExt = gapExt

			l.loadedCount += 1

			if !more {
				l.done = true
				l.scanner = nil
			}

			return true
		}

		if !more {
			l.done = true
			l.scanner = nil

			return false
		}
	}
}

func LoadAlignments(file io.Reader, addSkipState bool) <-chan *Alignment {
	loaded := make(chan *Alignment)

	go func() {
		if addSkipState {
			loaded <- SkipState()
		}

		meta := make(map[string]string)
		var seqs []string

		scanner := bufio.NewScanner(file)
		for {
			more := scanner.Scan()
			line := scanner.Text()

			if parsePreambleLine(line) {
				continue
			}

			if parseToolLine(line) != "" {
				continue
			}

			key, value, found := parseMetaLine(line)
			if found {
				meta[key] = value
				continue
			}

			_, seq, found := parseAlignmentLine(line)
			if found {
				seqs = append(seqs, seq)
				continue
			}

			if parseTerminatorLine(line) {
				if len(seqs) != 2 {
					panic("file format error - wrong number of sequences")
				}

				chromName, chromStart, chromStop, found := parseChromMeta(meta["TR"])
				if !found {
					panic("file format error - missing TR field")
				}

				// TODO: Write helper(s) to check for presence of key and then get the value

				score, err := strconv.Atoi(meta["SC"])
				if err != nil {
					panic("file format error - invalid SC field")
				}

				start, err := strconv.Atoi(meta["ST"])
				if err != nil {
					panic("file format error - invalid ST field")
				}

				stop, err := strconv.Atoi(meta["SP"])
				if err != nil {
					panic("file format error - invalid SP field")
				}

				consensusStart, err := strconv.Atoi(meta["CST"])
				if err != nil {
					panic("file format error - invalid CST field")
				}

				consensusStop, err := strconv.Atoi(meta["CSP"])
				if err != nil {
					panic("file format error - invalid CSP field")
				}

				var sequence string
				var subfamilySequence string
				if meta["TQ"] == "t" {
					sequence = reverse(seqs[0])
					subfamilySequence = reverse(seqs[1])
				} else {
					sequence = seqs[0]
					subfamilySequence = seqs[1]
				}

				flank, err := strconv.Atoi(meta["FL"])
				if err != nil {
					panic("file format error - invalid FL field")
				}

				gapInit, err := strconv.ParseFloat(meta["GI"], 64)
				if err != nil {
					panic("file format error - invalid GI field")
				}

				gapExt, err := strconv.ParseFloat(meta["GE"], 64)
				if err != nil {
					panic("file format error - invalid GE field")
				}

				loaded <- &Alignment{
					Subfamily:         meta["ID"],
					ChromName:         chromName,
					ChromStart:        chromStart,
					ChromStop:         chromStop,
					Score:             score,
					Start:             start,
					Stop:              stop,
					ConsensusStart:    consensusStart,
					ConsensusStop:     consensusStop,
					Sequence:          sequence,
					SubfamilySequence: subfamilySequence,
					Strand:            meta["SD"],
					Flank:             flank,
					SubMatrixName:     meta["MX"],
					GapInit:           gapInit,
					GapExt:            gapExt,
				}

				meta = make(map[string]string)
				seqs = nil
			}

			if !more {
				break
			}
		}

		close(loaded)
	}()

	return loaded
}

func reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}
