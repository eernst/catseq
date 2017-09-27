package cmd

import (
	"fmt"
	"os"
	"regexp"
	"runtime"
	"time"

	"github.com/eernst/catseq/pipeline"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"

	"github.com/spf13/cobra"
)

const (
	HeaderField string = "header"
	SeqField           = "seq"
	BothFields         = "both"
)

var grepField string
var invert bool
var ignoreCase bool

func init() {
	RootCmd.AddCommand(grepCmd)
	grepCmd.Flags().BoolP("fasta", "", false, "Input is in FASTA format.")
	grepCmd.Flags().BoolP("fastq", "", false, "Input is in FASTQ format.")
	grepCmd.Flags().StringP("field", "f", "header", "Which field to match the pattern against. One of \"header\",\"seq\", or \"both\".")
	grepCmd.Flags().BoolP("invert-match", "v", false, "Selected lines are those not matching any of the specified patterns.")
	grepCmd.Flags().BoolP("ignore-case", "i", false, "Perform case insensitive matching.")
}

func grepRecs(in <-chan *fastx.Record, regex *regexp.Regexp, grepField string) <-chan *fastx.Record {
	out := make(chan *fastx.Record)
	go func() {
		for rec := range in {
			if DEBUG {
				fmt.Fprintf(os.Stderr, "Matching %v against field %v containing %v ... ", regex.String(), grepField, rec.Name)
			}

			switch grepField {
			case HeaderField:
				if regex.Match(rec.Name) != invert {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "Matched!\n")
					}
					out <- rec
				} else {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "No match!\n")
					}
					out <- nil
				}
			case SeqField:
				if regex.Match(rec.Seq.Seq) != invert {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "Matched!\n")
					}
					out <- rec
				} else {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "No match!\n")
					}
					out <- nil
				}
			case BothFields:
				if (regex.Match(rec.Name) || regex.Match(rec.Seq.Seq)) != invert {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "Matched!\n")
					}
					out <- rec
				} else {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "No match!\n")
					}
					out <- nil
				}
			default:
				fmt.Fprintf(os.Stderr, "Error: Unknown grep field.")
				os.Exit(1)
			}

		}
		close(out)
	}()
	return out
}

var grepCmd = &cobra.Command{
	Use:   "grep PATTERN [SEQUENCE_FILE]",
	Short: "Match a regular expression in sequences from (multi-)sequence files.",
	Long: `
	
grep scans through input sequences and outputs only those that match the 
provided pattern. By default, all data is considered for a match, but the
search can be limited to the header, sequence, quality, or quality header.

The pattern language is the same as the regular expression syntax used by Perl,
Python, etc. Reference: https://golang.org/s/re2syntax

FASTQ and FASTA formats are currently supported and guessed based on file 
extension. Seqeunce can be piped in on STDIN, in which case the format must be
specified.
`,
	Run: func(cmd *cobra.Command, args []string) {
		StartProfiling()

		flags := cmd.Flags()
		var err error

		invert, err = flags.GetBool("invert-match")
		check(err)
		ignoreCase, err = flags.GetBool("ignore-case")
		check(err)

		// seqsFile is the multi-fast(a/q) over which we will iterate
		var seqsInFileName string

		var pattern string
		switch {
		case len(args) == 0:
			fmt.Fprintf(os.Stderr, "Error: Can't grep without a pattern.\n")
			cmd.Usage()
			os.Exit(1)
		case len(args) == 1:
			pattern = args[0]
			// TODO: Check here for valid sequence on stdin
			seqsInFileName = "-"
			fmt.Fprintf(os.Stderr, "No input sequence file given. Reading from STDIN.\n")
		case len(args) == 2:
			pattern = args[0]
			seqsInFileName = args[1]
		default:
			fmt.Fprintf(os.Stderr, "Error: Wrong number of positional arguments given.\n")
			cmd.Usage()
			os.Exit(1)
		}
		if DEBUG && len(args) > 0 {
			fmt.Fprintf(os.Stderr, "grep called with args[]:\n %v\n", args)
		}

		if ignoreCase {
			// perhaps faster to just UC the string
			// http://stackoverflow.com/questions/15326421/how-do-i-do-a-case-insensitive-regular-expression-in-go
			pattern = "(?i)" + pattern
		}
		regex, err := regexp.Compile(pattern)
		check(err)

		seq.ValidateSeq = false
		reader, err := fastx.NewDefaultReader(seqsInFileName)
		check(err)

		writer, err := xopen.Wopen("-") // "-" for STDOUT
		check(err)
		defer writer.Close()

		grepField, err = flags.GetString("field")
		check(err)

		// Using the pipeline pattern
		inStream := pipeline.ChannelRec(reader)
		processors := make([]<-chan *fastx.Record, runtime.GOMAXPROCS(0))
		for p := range processors {
			processors[p] = grepRecs(inStream, regex, grepField)
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
