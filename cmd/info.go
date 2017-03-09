package cmd

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strings"
	"time"

	"github.com/eernst/catseq/seqmath"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
)

func init() {
	RootCmd.AddCommand(infoCmd)

	infoCmd.Flags().BoolP("fasta", "", false, "Input is in FASTA format.")
	infoCmd.Flags().BoolP("fastq", "", false, "Input is in FASTQ format.")
	infoCmd.Flags().BoolP("summary", "s", false, "Only output summary info for all sequences.")
}

var infoCmd = &cobra.Command{
	Use:   "info SEQUENCE_FILE",
	Short: "Show basic sequence info.",
	Long: `
	
Print basic sequence info including name, length, GC content, average quality,
etc. in a tabular format, one input sequence per row.

FASTQ and FASTA formats are currently supported and guessed based on file 
extension. Seqeunce can be piped in on STDIN, in which case the format must be
specified.`,
	Run: func(cmd *cobra.Command, args []string) {
		StartProfiling()

		flags := cmd.Flags()
		summaryOnly, err := flags.GetBool("summary")
		check(err)
		summaryOut := os.Stderr
		if summaryOnly {
			summaryOut = os.Stdout
		}

		// seqsFile is the multi-fast(a/q) over which we will iterate
		var seqsInFileName string

		if len(args) < 1 {
			// TODO: Can we check if there's something actually coming in on STDIN?
			seqsInFileName = "/dev/stdin"
			fmt.Fprintf(os.Stderr, "No input sequence file given. Reading from STDIN.\n")
		} else {
			seqsInFileName = args[0]
		}
		if DEBUG && len(args) > 1 {
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
				fasta.NewReader(bufio.NewReader(seqsInFile), linear.NewSeq("", nil, alphabet.DNAredundant)))
		case FastqFormat:
			seqsIn = seqio.NewScanner(
				fastq.NewReader(bufio.NewReader(seqsInFile), linear.NewQSeq("", nil, alphabet.DNAredundant, alphabet.Sanger)))
		case UnknownFormat:
			fmt.Fprintf(os.Stderr, "Unknown input sequence file format.")
			os.Exit(1)
		}

		if PrintHeader {
			fmt.Fprintf(os.Stdout, "accession\tlength\tgc-content\tmean quality\tmean P(error)\t")
		}

		totalSeqs := 0
		totalGcCount := 0
		totalSeqLength := 0
		var sumBaseQualityScores int
		var sumMeanQualityScores float64
		var sumBaseErrorProbs float64
		var sumMeanErrorProbs float64
		var seqLens []int

		for seqsIn.Next() {
			s := seqsIn.Seq()

			var seqBytes []byte = make([]byte, s.Len())

			qualScores := 0
			errorProbs := float64(0)
			switch s.(type) {
			case *linear.QSeq:
				for i, ql := range s.(*linear.QSeq).Seq {
					seqBytes[i] = byte(ql.L)
					qualScores += int(ql.Q)
					errorProbs += seqmath.ErrorProbForQ(int(ql.Q))
				}
			case *linear.Seq:
				for i, l := range s.(*linear.Seq).Seq {
					seqBytes[i] = byte(l)
				}
			}

			seqStr := string(seqBytes)
			gcCount := strings.Count(seqStr, "C") + strings.Count(seqStr, "G") + strings.Count(seqStr, "c") + strings.Count(seqStr, "g")
			gcRatio := float64(gcCount) / float64(s.Len())
			gcPercent := gcRatio * 100

			totalSeqs += 1
			totalGcCount += gcCount
			totalSeqLength += s.Len()

			meanBaseQual := float64(qualScores) / float64(s.Len())
			meanErrorProb := float64(errorProbs) / float64(s.Len())
			sumMeanQualityScores += meanBaseQual
			sumMeanErrorProbs += meanErrorProb

			sumBaseErrorProbs += errorProbs
			sumBaseQualityScores += qualScores

			seqLens = append(seqLens, s.Len())

			// Print per-read info
			if !summaryOnly {
				fmt.Fprintf(os.Stdout, "%s\t%d\t%.2f", s.Name(), s.Len(), gcPercent)

				if seqsInFormat == "fastq" {
					fmt.Fprintf(os.Stdout, "\t%.2f\t%.4f", meanBaseQual, meanErrorProb)
				}

				fmt.Fprintf(os.Stdout, "\n")
			}
		}

		totalGcRatio := float64(totalGcCount) / float64(totalSeqLength)
		totalGcPercent := totalGcRatio * 100

		meanQualityPerSeq := float64(sumMeanQualityScores) / float64(totalSeqs)
		meanErrorProbPerSeq := float64(sumMeanErrorProbs) / float64(totalSeqs)

		meanQualityPerBase := float64(sumBaseQualityScores) / float64(totalSeqLength)
		meanErrorProbPerBase := float64(sumBaseErrorProbs) / float64(totalSeqLength)

		// N50 Calc
		sort.Ints(seqLens)
		//seqLens = sort.Reverse(seqLens)
		var seqLensCum []int
		var cumLen int

		cumLen = 0
		var n50 int
		for i, l := range seqLens {
			if i == 0 {
				seqLensCum = append(seqLensCum, seqLens[i])
			} else {
				seqLensCum = append(seqLensCum, seqLens[i]+seqLens[i-1])

			}
			//fmt.Fprintf(os.Stderr, "seqLensCum[%v] = %v\n", i, seqLensCum[i])
			newLen := cumLen + l
			if newLen > totalSeqLength/2 {
				break
			} else {
				n50 = l
				cumLen = newLen
			}
		}

		cumLen = 0
		var n75 int
		for i, l := range seqLens {
			if i == 0 {
				seqLensCum = append(seqLensCum, seqLens[i])
			} else {
				seqLensCum = append(seqLensCum, seqLens[i]+seqLens[i-1])

			}
			//fmt.Fprintf(os.Stderr, "seqLensCum[%v] = %v\n", i, seqLensCum[i])
			newLen := cumLen + l
			if newLen > totalSeqLength/4 {
				break
			} else {
				n75 = l
				cumLen = newLen
			}
		}

		const sep string = "--------------------\n"
		fmt.Fprintf(summaryOut, "\nSUMMARY\n"+sep)
		fmt.Fprintf(summaryOut, "Total Seqs (#): %23d\n", totalSeqs)
		fmt.Fprintf(summaryOut, "Total Length (bp): %20d\n", totalSeqLength)
		fmt.Fprintf(summaryOut, "Overall GC Content (%%): %15.2f\n", totalGcPercent)
		fmt.Fprintf(summaryOut, "Shortest (bp): %24d\n", seqLens[0])
		fmt.Fprintf(summaryOut, "Longest (bp): %25d\n", seqLens[len(seqLens)-1])
		fmt.Fprintf(summaryOut, "N50 (bp): %29d\n", n50)
		fmt.Fprintf(summaryOut, "N75 (bp): %29d\n", n75)
		// Would be better to test the Sequence for quality info
		if seqsInFormat == "fastq" {
			fmt.Fprintf(summaryOut, "\nPER-SEQ\n"+sep)
			fmt.Fprintf(summaryOut, "Mean Phred quality score: %13.2f\n", meanQualityPerSeq)
			fmt.Fprintf(summaryOut, "Mean error rate: %22.4f\n", meanErrorProbPerSeq)
			fmt.Fprintf(summaryOut, "\nPER-BASE\n"+sep)
			fmt.Fprintf(summaryOut, "Mean Phred quality score: %13.2f\n", meanQualityPerBase)
			fmt.Fprintf(summaryOut, "Mean error rate: %22.4f\n", meanErrorProbPerBase)
		}

		time.Sleep(0 * time.Millisecond)

		StopProfiling()
	},
}
