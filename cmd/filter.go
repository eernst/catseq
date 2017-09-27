package cmd

import (
	"fmt"
	"os"
	"runtime"
	"time"

	"github.com/eernst/catseq/pipeline"
	"github.com/eernst/catseq/seqmath"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

var minLength int
var maxLength int

func init() {
	RootCmd.AddCommand(filterCmd)
	filterCmd.Flags().IntP("length_min", "", -1, "Minimum sequence length to keep. [0]")
	filterCmd.Flags().IntP("length_max", "", -1, "Maximum sequence length to keep. [∞]")
	filterCmd.Flags().Float64P("error_rate_avg_min", "", 0, "Keep reads with a mean error rate equal to or greater than this. [0.00]")
	filterCmd.Flags().Float64P("error_rate_avg_max", "", 1, "Keep reads with a mean error rate equal to or less than this. [1.00]")
	filterCmd.Flags().Float64P("qual_avg_min", "", 0, "Keep reads with a mean phred base quality equal to or greater than this. [0.00]")
	filterCmd.Flags().Float64P("qual_avg_max", "", -1, "Keep reads with a mean phred base quality equal to or greater than this. [∞]")
}

func passesFilters(s *seq.Seq, flags *pflag.FlagSet) bool {
	minLength, _ := flags.GetInt("length_min")
	maxLength, _ := flags.GetInt("length_max")
	minMeanError, _ := flags.GetFloat64("error_rate_avg_min")
	maxMeanError, _ := flags.GetFloat64("error_rate_avg_max")
	minMeanQ, _ := flags.GetFloat64("qual_avg_min")
	maxMeanQ, _ := flags.GetFloat64("qual_avg_max")

	// Validate argument values

	switch {
	case minLength >= 0 && s.Length() < minLength:
		return false
	case maxLength >= 0 && s.Length() > maxLength:
		return false
	}

	if len(s.Qual) > 0 {
		if len(s.QualValue) <= 0 {
			vals, err := seq.QualityValue(seq.Sanger, s.Qual)
			s.QualValue = vals
			check(err)
		}
		var qualScores int = 0
		var errorProbs float64 = 0

		for _, score := range s.QualValue {
			qualScores += score
			errorProbs += seqmath.ErrorProbForQ(score)
		}

		meanQ := float64(qualScores) / float64(s.Length())
		meanErrorProb := float64(errorProbs) / float64(s.Length())

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

func filterSeq(in <-chan *fastx.Record, flags *pflag.FlagSet) <-chan *fastx.Record {
	out := make(chan *fastx.Record)
	go func() {
		for rec := range in {
			if passesFilters(rec.Seq, flags) {
				if DEBUG {
					fmt.Fprintf(os.Stderr, "PASSED FILTER   Acc: %s		Length: %d\n", rec.Name, rec.Seq.Length())
				}
				out <- rec
			}
			out <- nil
		}
		close(out)
	}()
	return out
}

var filterCmd = &cobra.Command{
	Use:   "filter SEQUENCE_FILE",
	Short: "Filter sequences from (multi-)sequence files.",
	Long: `
	
Filter input sequences by applying combinations of simple criteria.

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
			seqsInFileName = "-"
			fmt.Fprintf(os.Stderr, "No input sequence file given. Reading from STDIN.\n")
		} else {
			seqsInFileName = args[0]
		}
		if DEBUG && len(args) > 0 {
			fmt.Fprintf(os.Stderr, "args[0] is %q\n", args[0])
		}

		seq.ValidateSeq = false
		reader, err := fastx.NewDefaultReader(seqsInFileName)
		check(err)

		writer, err := xopen.Wopen("-") // "-" for STDOUT
		check(err)
		defer writer.Close()

		// Using the pipeline pattern
		inStream := pipeline.ChannelRec(reader)
		processors := make([]<-chan *fastx.Record, runtime.GOMAXPROCS(0))
		for p := range processors {
			processors[p] = filterSeq(inStream, flags)
		}
		for record := range pipeline.MergeRec(processors...) {
			if record != nil {
				record.FormatToWriter(writer, 0)
			}
		}

		time.Sleep(0 * time.Millisecond)

		StopProfiling()
	},
}
