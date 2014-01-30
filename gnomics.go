package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
)

func gcContent()

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
		fmt.Fprintf(os.Stderr,"Error: Incorrect number of arguments.\n")
		os.Exit(1)
	}

	runtime.GOMAXPROCS(flagNumProcs)

	if flagVerbose {
		fmt.Fprintf(os.Stderr,"verbose: %t\n", flagVerbose)
		fmt.Fprintf(os.Stderr,"using %d/%d available procs \n", flagNumProcs, runtime.NumCPU())
		fmt.Fprintf(os.Stderr,"trailing args: %s\n", flag.Args())
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
		if flagVerbose {
			fmt.Fprintf(os.Stderr, "%s\n", ctgs.Seq().Name())
		}
		
		
	}

}
