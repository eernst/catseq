package cmd

import (
	"bufio"
	"fmt"
	//"log"
	"os"
	"runtime"
	"sync"
	"time"

	"github.com/eernst/catseq/seqmath"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

var minLength int
var maxLength int

func init() {
	RootCmd.AddCommand(filterCmd)

	filterCmd.Flags().BoolP("fasta", "", false, "Input is in FASTA format.")
	filterCmd.Flags().BoolP("fastq", "", false, "Input is in FASTQ format.")
	filterCmd.Flags().IntP("length_min", "", -1, "Minimum sequence length to keep. [0]")
	filterCmd.Flags().IntP("length_max", "", -1, "Maximum sequence length to keep. [∞]")
	filterCmd.Flags().Float64P("error_rate_avg_min", "", 0, "Keep reads with a mean error rate equal to or greater than this. [0.00]")
	filterCmd.Flags().Float64P("error_rate_avg_max", "", 1, "Keep reads with a mean error rate equal to or less than this. [1.00]")
	filterCmd.Flags().Float64P("qual_avg_min", "", 0, "Keep reads with a mean phred base quality equal to or greater than this. [0.00]")
	filterCmd.Flags().Float64P("qual_avg_max", "", -1, "Keep reads with a mean phred base quality equal to or greater than this. [∞]")
}

func passesFilters(s seq.Sequence, flags *pflag.FlagSet) bool {
	minLength, _ := flags.GetInt("length_min")
	maxLength, _ := flags.GetInt("length_max")
	minMeanError, _ := flags.GetFloat64("error_rate_avg_min")
	maxMeanError, _ := flags.GetFloat64("error_rate_avg_max")
	minMeanQ, _ := flags.GetFloat64("qual_avg_min")
	maxMeanQ, _ := flags.GetFloat64("qual_avg_max")

	// Validate argument values

	switch {
	case minLength >= 0 && s.Len() < minLength:
		return false
	case maxLength >= 0 && s.Len() > maxLength:
		return false
	}

	switch s.(type) {
	case *linear.QSeq:
		var seqBytes []byte = make([]byte, s.Len())
		var qualScores int = 0
		var errorProbs float64 = 0

		for i, ql := range s.(*linear.QSeq).Seq {
			seqBytes[i] = byte(ql.L)
			qualScores += int(ql.Q)
			errorProbs += seqmath.ErrorProbForQ(int(ql.Q))
		}

		meanQ := float64(qualScores) / float64(s.Len())
		meanErrorProb := float64(errorProbs) / float64(s.Len())

		switch {
		case meanErrorProb < minMeanError:
			return false
		case meanErrorProb > maxMeanError:
			return false
		case minMeanQ >= 0 && meanQ < minMeanQ:
			return false
		case maxMeanQ >= 0 && meanQ > maxMeanQ:
			return false
		}
	}

	return true
}

func channelSeq(seqsIn *seqio.Scanner) <-chan seq.Sequence {
	out := make(chan seq.Sequence)
	go func() {
		for seqsIn.Next() {
			out <- seqsIn.Seq()
		}
		close(out)
	}()
	return out
}

func filterSeq(in <-chan seq.Sequence, flags *pflag.FlagSet) <-chan seq.Sequence {
	out := make(chan seq.Sequence)
	go func() {
		for seq := range in {
			if passesFilters(seq, flags) {
				if DEBUG {
					fmt.Fprintf(os.Stderr, "PASSED FILTER   Acc: %s		Length: %d\n", seq.Name(), seq.Len())
				}
				out <- seq
			}
			out <- nil
		}
		close(out)
	}()
	return out
}

// Copied from https://blog.golang.org/pipelines
func merge(cs ...<-chan seq.Sequence) chan seq.Sequence {
	var wg sync.WaitGroup
	out := make(chan seq.Sequence)

	// Start an output goroutine for each input channel in cs.  output
	// copies values from c to out until c is closed, then calls wg.Done.
	output := func(c <-chan seq.Sequence) {
		for n := range c {
			out <- n
		}
		wg.Done()
	}
	wg.Add(len(cs))
	for _, c := range cs {
		go output(c)
	}

	// Start a goroutine to close out once all the output goroutines are
	// done.  This must start after the wg.Add call.
	go func() {
		wg.Wait()
		close(out)
	}()
	return out
}

var filterCmd = &cobra.Command{
	Use:   "filter SEQUENCE_FILE",
	Short: "Filter sequences from (multi-)sequence files.",
	Long: `
	
filter input sequences by applying combinations of simple criteria.

FASTQ and FASTA formats are currently supported and guessed based on file 
extension. Seqeunce can be piped in on STDIN, in which case the format must be
specified.`,
	Run: func(cmd *cobra.Command, args []string) {
		StartProfiling()

		flags := cmd.Flags()

		// seqsFile is the multi-fast(a/q) over which we will iterate
		var seqsInFileName string

		if len(args) < 1 {
			// Anything coming in on stdin?
			seqsInFileName = "/dev/stdin"
			fmt.Fprintf(os.Stderr, "No input sequence file given. Reading from STDIN.\n")
		} else {
			seqsInFileName = args[0]
		}
		if DEBUG && len(args) > 0 {
			fmt.Fprintf(os.Stderr, "args[0] is %q\n", args[0])
		}

		var seqsInFormat string
		isFasta, err := flags.GetBool("fasta")
		check(err)
		isFastq, err := flags.GetBool("fastq")
		check(err)
		if isFasta {
			seqsInFormat = FastaFormat
		} else if isFastq {
			seqsInFormat = FastqFormat
		} else {
			seqsInFormat, err = GuessFileFormat(seqsInFileName)
		}
		check(err)

		seqsInFile, err := os.Open(seqsInFileName)
		check(err)
		defer seqsInFile.Close()

		var seqsIn *seqio.Scanner
		switch seqsInFormat {
		case FastaFormat:
			seqsIn = seqio.NewScanner(
				fasta.NewReader(bufio.NewReader(seqsInFile), linear.NewSeq("", nil, alphabet.DNA)))
		case FastqFormat:
			seqsIn = seqio.NewScanner(
				fastq.NewReader(bufio.NewReader(seqsInFile), linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger)))
		case UnknownFormat:
			fmt.Fprintf(os.Stderr, "Unknown input sequence file format.")
			os.Exit(1)
		}

		var underWriter *bufio.Writer = bufio.NewWriter(os.Stdout)
		var writer seqio.Writer
		switch seqsInFormat {
		case FastaFormat:
			writer = fasta.NewWriter(underWriter, MaxInt)
		case FastqFormat:
			writer = fastq.NewWriter(underWriter)
		}
		// TODO: request addition of Flush() to biogo's writers
		defer underWriter.Flush()

		inStream := channelSeq(seqsIn)

		processors := make([]<-chan seq.Sequence, runtime.GOMAXPROCS(0))
		for p := range processors {
			processors[p] = filterSeq(inStream, flags)
		}

		for seq := range merge(processors...) {
			if seq != nil {
				writer.Write(seq)
			}
		}

		time.Sleep(0 * time.Millisecond)

		StopProfiling()
	},
}
