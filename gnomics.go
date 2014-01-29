package main

import "os"
import "fmt"
import "flag"


//import "code.google.com/p/biogo"

func main() {

	flag.Usage = func() {
	    fmt.Fprintf(os.Stderr, "usage: %s [options] scaffolds.fasta\n\n", os.Args[0])
		flag.PrintDefaults()
	}
	
	flagVerbosePtr := flag.Bool("v", false, "verbose")
	flag.Parse()
	
	if (*flagVerbosePtr) {
		fmt.Println("verbose:", *flagVerbosePtr)
	}
//	fmt.Println("Args from flag.Args:", flag.Args())
}