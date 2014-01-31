package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime"
	"strings"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
)

const (
	DEBUG = false
)

func calgc(seq *linear.Seq) (gc float64, err error) {
	seqStr := string(seq.Seq.String())
	gccount := strings.Count(seqStr, "C") + strings.Count(seqStr, "G")
	gcratio := float64(gccount) / float64(seq.Len())
	return gcratio, nil
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

var flagVerbose bool
var flagNumProcs int

func init() {
	flag.BoolVar(&flagVerbose, "v", false, "Print verbose output")
	flag.IntVar(&flagNumProcs, "p", 1, "Use up to this many processors")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: %s [options] scaffolds.fasta\n\n", os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()

	if flag.NArg() < 1 {
		flag.Usage()
		fmt.Fprintf(os.Stderr, "Error: Incorrect number of arguments.\n")
		os.Exit(1)
	}

	runtime.GOMAXPROCS(flagNumProcs)

	if flagVerbose {
		fmt.Fprintf(os.Stderr, "verbose: %t\n", flagVerbose)
		fmt.Fprintf(os.Stderr, "using %d/%d available procs \n", flagNumProcs, runtime.NumCPU())
		fmt.Fprintf(os.Stderr, "trailing args: %s\n", flag.Args())
	}
}

func main() {
	var ctgsFileName string
	if flag.Arg(0) != "" {
		ctgsFileName = flag.Arg(0)
	} else {
		ctgsFileName = "test.fasta"
	}

	ctgsFile, err := os.Open(ctgsFileName)
	check(err)
	defer ctgsFile.Close()
	ctgs := seqio.NewScanner(
		fasta.NewReader(bufio.NewReader(ctgsFile), linear.NewSeq("", nil, alphabet.DNA)))

	for ctgs.Next() {
		gc, err := calgc(ctgs.Seq().(*linear.Seq))
		if err != nil {
			panic(err)
		}

		if DEBUG {
			fmt.Fprintf(os.Stderr, "Acc:%s   Seq:%s   GC:%f\n", ctgs.Seq().Name(), ctgs.Seq().(*linear.Seq).String(), gc)
		}
		
		fmt.Printf("%s\t%f\n", ctgs.Seq().Name(), gc)
	}

}
