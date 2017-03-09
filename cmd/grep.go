package cmd

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"time"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq/linear"

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
			// Anything coming in on stdin?
			seqsInFileName = "/dev/stdin"
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

		var seqsInFormat string
		flagFasta, err := flags.GetBool("fasta")
		check(err)
		flagFastq, err := flags.GetBool("fastq")
		check(err)
		if flagFasta {
			seqsInFormat = FastaFormat
		} else if flagFastq {
			seqsInFormat = FastqFormat
		} else {
			seqsInFormat, err = GuessFileFormat(seqsInFileName)
		}
		check(err)
		if DEBUG {
			fmt.Fprintf(os.Stderr, "seqsInFormat: %v\n", seqsInFormat)
		}

		seqsInFile, err := os.Open(seqsInFileName)
		check(err)
		defer seqsInFile.Close()

		var seqsIn *seqio.Scanner
		var underWriter *bufio.Writer
		var writer seqio.Writer
		switch seqsInFormat {
		case FastaFormat:
			seqsIn = seqio.NewScanner(
				fasta.NewReader(bufio.NewReader(seqsInFile), linear.NewSeq("", nil, alphabet.DNA)))
			underWriter = bufio.NewWriter(os.Stdout)
			writer = fasta.NewWriter(underWriter, MaxInt)
			defer underWriter.Flush()
		case FastqFormat:
			seqsIn = seqio.NewScanner(
				fastq.NewReader(bufio.NewReader(seqsInFile), linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger)))
			writer = fastq.NewWriter(bufio.NewWriter(os.Stdout))
		case UnknownFormat:
			fmt.Fprintf(os.Stderr, "Unknown input sequence file format.\n")
			os.Exit(1)
		}

		grepField, err = flags.GetString("field")
		check(err)
		switch grepField {
		case HeaderField:
			for seqsIn.Next() {
				seq := seqsIn.Seq()

				if DEBUG {
					fmt.Fprintf(os.Stderr, "Matching %v against field %v containing %v%v ... ", pattern, grepField, seq.Name(), seq.Description())
				}

				if (regex.MatchString(seq.Name()) || regex.MatchString(seq.Description())) != invert {
					bytesWritten, err := writer.Write(seq)
					check(err)

					if DEBUG {
						fmt.Fprintf(os.Stderr, "Matched!\n")
						fmt.Fprintf(os.Stderr, "Wrote %v bytes.\n", bytesWritten)
					}

				} else {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "No match!\n")
					}
				}
			}
		case SeqField:
			for seqsIn.Next() {
				seq := seqsIn.Seq()

				if DEBUG {
					fmt.Fprintf(os.Stderr, "Matching %v against field %v containing %v%v ... ", pattern, grepField, seq.Name(), seq.Description())
				}

				if (regex.Match(alphabet.LettersToBytes(seq.(*linear.Seq).Seq))) != invert {
					bytesWritten, err := writer.Write(seq)
					check(err)

					if DEBUG {
						fmt.Fprintf(os.Stderr, "Matched!\n")
						fmt.Fprintf(os.Stderr, "Wrote %v bytes.\n", bytesWritten)
					}

				} else {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "No match!\n")
					}
				}
			}
		case BothFields:
			for seqsIn.Next() {
				seq := seqsIn.Seq()

				if DEBUG {
					fmt.Fprintf(os.Stderr, "Matching %v against field %v containing %v%v ... ", pattern, grepField, seq.Name(), seq.Description())
				}

				if (regex.MatchString(seq.Name()) || regex.MatchString(seq.Description())) != invert ||
					regex.Match(alphabet.LettersToBytes(seq.(*linear.Seq).Seq)) {
					bytesWritten, err := writer.Write(seq)
					check(err)

					if DEBUG {
						fmt.Fprintf(os.Stderr, "Matched!\n")
						fmt.Fprintf(os.Stderr, "Wrote %v bytes.\n", bytesWritten)
					}

				} else {
					if DEBUG {
						fmt.Fprintf(os.Stderr, "No match!\n")
					}
				}
			}
		default:
			fmt.Fprintf(os.Stderr, "Error: Unknown grep field.")
			os.Exit(1)
		}

		time.Sleep(0 * time.Millisecond)

		StopProfiling()
	},
}
