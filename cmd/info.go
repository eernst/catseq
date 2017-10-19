package cmd

import (
	"fmt"
	"os"
	"runtime"
	"sort"
	"sync"
	"time"

	//"github.com/eernst/catseq/pipeline"
	"github.com/eernst/catseq/seqmath"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"

	"github.com/spf13/cobra"
)

type InfoRecord struct {
	Record        *fastx.Record
	GcBases       int
	AtBases       int
	NonATGCNBases int
	NBases        int
	GcRatio       float64
	MeanBaseQual  float64
	MeanErrorProb float64
	SumQ          int
	SumErrorProbs float64
}

func init() {
	RootCmd.AddCommand(infoCmd)
	infoCmd.Flags().BoolP("fasta", "", false, "Input is in FASTA format.")
	infoCmd.Flags().BoolP("fastq", "", false, "Input is in FASTQ format.")
	infoCmd.Flags().BoolP("summary", "s", false, "Only output summary info for all sequences.")
}

func infoSeq(in <-chan fastx.RecordChunk) <-chan *InfoRecord {
	out := make(chan *InfoRecord)
	go func() {
		for chunk := range in {
			for _, rec := range chunk.Data {
				s := rec.Seq
				length := s.Length()

				var gcBases int = 0
				var atBases int = 0
				var nonATGCNBases int = 0
				var nBases int = 0

				for _, char := range s.Seq {
					switch char {
					case 'C', 'c', 'G', 'g', 'S', 's':
						gcBases++
					case 'A', 'a', 'T', 't', 'W', 'w':
						atBases++
					case 'N':
						nBases++
					default:
						nonATGCNBases++
					}
				}

				gcRatio := float64(gcBases) / float64(length-(nonATGCNBases+nBases))

				var qualScores int = 0
				var errorProbs float64 = 0
				var meanBaseQual float64
				var meanErrorProb float64

				if len(s.Qual) > 0 {
					if len(s.QualValue) <= 0 {
						vals, err := seq.QualityValue(seq.Sanger, s.Qual)
						s.QualValue = vals
						check(err)
					}

					for _, score := range s.QualValue {
						qualScores += score
						errorProbs += seqmath.ErrorProbForQ(score)
					}

					meanBaseQual = float64(qualScores) / float64(length)
					meanErrorProb = float64(errorProbs) / float64(length)
				}

				infoRec := &InfoRecord{
					Record:        rec,
					GcBases:       gcBases,
					AtBases:       atBases,
					NonATGCNBases: nonATGCNBases,
					NBases:        nBases,
					GcRatio:       gcRatio,
					MeanBaseQual:  meanBaseQual,
					MeanErrorProb: meanErrorProb,
					SumQ:          qualScores,
					SumErrorProbs: errorProbs}

				out <- infoRec
			}
		}
		close(out)
	}()
	return out
}

func mergeInfoRec(chans ...<-chan *InfoRecord) chan *InfoRecord {
	var wg sync.WaitGroup
	out := make(chan *InfoRecord)
	output := func(c <-chan *InfoRecord) {
		for n := range c {
			out <- n
		}
		wg.Done()
	}
	wg.Add(len(chans))
	for _, c := range chans {
		go output(c)
	}
	go func() {
		wg.Wait()
		close(out)
	}()
	return out
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
			// TODO: Check here for valid sequence on stdin
			seqsInFileName = "-"
			fmt.Fprintf(os.Stderr, "No input sequence file given. Reading from STDIN.\n")
		} else {
			seqsInFileName = args[0]
		}
		if DEBUG && len(args) > 1 {
			fmt.Fprintf(os.Stderr, "args[0] is %q\n", args[0])
		}

		seq.ValidateSeq = false
		reader, err := fastx.NewDefaultReader(seqsInFileName)
		check(err)

		if PrintHeader {
			fmt.Fprintf(os.Stdout, "accession\tlength\tgc-content\tmean quality\tmean P(error)\t\n")
		}

		chunkStream := reader.ChunkChan(runtime.GOMAXPROCS(0), 1<<8)
		processors := make([]<-chan *InfoRecord, runtime.GOMAXPROCS(0))
		for p := range processors {
			processors[p] = infoSeq(chunkStream)
		}

		var totalSeqs int
		var totalSeqLength int
		var totalGcCount int
		var totalNonATGCNBases int
		var totalNBases int
		var sumBaseQualityScores int
		var sumMeanQualityScores float64
		var sumBaseErrorProbs float64
		var sumMeanErrorProbs float64
		var seqLens []int

		for infoRec := range mergeInfoRec(processors...) {
			if infoRec != nil {

				rec := infoRec.Record
				s := rec.Seq
				length := s.Length()

				// Print per-read info
				if !summaryOnly {
					fmt.Fprintf(os.Stdout, "%s\t%d\t%.2f", rec.Name, length, infoRec.GcRatio*100)

					if reader.IsFastq {
						fmt.Fprintf(os.Stdout, "\t%.2f\t%.4f", infoRec.MeanBaseQual, infoRec.MeanErrorProb)
					}

					fmt.Fprintf(os.Stdout, "\n")
				}

				totalSeqs += 1
				totalGcCount += infoRec.GcBases
				totalSeqLength += length
				totalNonATGCNBases += infoRec.NonATGCNBases
				totalNBases += infoRec.NBases
				sumMeanQualityScores += infoRec.MeanBaseQual
				sumMeanErrorProbs += infoRec.MeanErrorProb

				sumBaseQualityScores += infoRec.SumQ
				sumBaseErrorProbs += infoRec.SumErrorProbs

				seqLens = append(seqLens, length)
			}
		}

		totalGcRatio := float64(totalGcCount) / float64(totalSeqLength)
		totalGcPercent := totalGcRatio * 100

		totalGcRatioNoAmbig := float64(totalGcCount) / float64(totalSeqLength-(totalNonATGCNBases+totalNBases))
		totalGcPercentNoAmbig := totalGcRatioNoAmbig * 100

		meanQualityPerSeq := float64(sumMeanQualityScores) / float64(totalSeqs)
		meanErrorProbPerSeq := float64(sumMeanErrorProbs) / float64(totalSeqs)

		meanQualityPerBase := float64(sumBaseQualityScores) / float64(totalSeqLength)
		meanErrorProbPerBase := float64(sumBaseErrorProbs) / float64(totalSeqLength)

		// NXX Calc
		sort.Ints(seqLens)
		var nxx []int = seqmath.Nxx(seqLens, totalSeqLength)

		const sep string = "--------------------\n"
		fmt.Fprintf(summaryOut, "\nSUMMARY\n"+sep)
		fmt.Fprintf(summaryOut, "Total Seqs (#): %23d\n", totalSeqs)
		fmt.Fprintf(summaryOut, "Total Length (bp): %20d\n", totalSeqLength)
		fmt.Fprintf(summaryOut, "GC Content (%%): %23.2f\n", totalGcPercent)
		fmt.Fprintf(summaryOut, "GC Content (%%, no ambig): %13.2f\n", totalGcPercentNoAmbig)
		fmt.Fprintf(summaryOut, "N bases (#): %26d\n", totalNBases)
		fmt.Fprintf(summaryOut, "Non-ATGCN bases (#): %18d\n", totalNonATGCNBases)
		fmt.Fprintf(summaryOut, "Shortest (bp): %24d\n", seqLens[0])
		fmt.Fprintf(summaryOut, "Longest (bp): %25d\n", seqLens[len(seqLens)-1])
		fmt.Fprintf(summaryOut, "Mean (bp): %28d\n", totalSeqLength/totalSeqs)
		fmt.Fprintf(summaryOut, "Median (bp): %26d\n", seqmath.Median(seqLens))
		fmt.Fprintf(summaryOut, "N80 (bp): %29d\n", nxx[80])
		fmt.Fprintf(summaryOut, "N50 (bp): %29d\n", nxx[50])
		fmt.Fprintf(summaryOut, "N20 (bp): %29d\n", nxx[20])
		// Would be better to test the Sequence for quality info
		if reader.IsFastq {
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
