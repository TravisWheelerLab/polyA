package main

import (
	"github.com/traviswheelerlab/polyA"
	"os"
)

func main() {
	stoPath := os.Args[1]
	stoFile, err := os.Open(stoPath)
	if err != nil {
		panic("failed to open file")
	}

	count := 0
	//for _ = range polyA.LoadAlignments(stoFile, true) {
	//	count += 1
	//}

	loader := polyA.NewLoader(stoFile, true)
	a := &polyA.Alignment{}
	for loader.Load(a) {
		count += 1
	}

	println(count)
}
