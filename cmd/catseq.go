package cmd

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var RootCmd = &cobra.Command{
	Use:   "catseq",
	Short: "catseq is a toolbox for performing common operations on sequence data.",
	Long:  `catseq is a toolbox for working with genomes, annotations, sequencing reads and the like.`,
	Run: func(cmd *cobra.Command, args []string) {
		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "Please provide a command to run.\n\n")
			flag.Usage()
			os.Exit(1)
		}
		//utils.StopOnErr(RootCmd.Execute())
	},
}

const (
	DEBUG = false
)

const MaxInt = math.MaxInt32

type GoseqCommand string

const (
	// The various commands catseq can perform
	INFO    GoseqCommand = "info"
	FILTER               = "filter"
	VERSION              = "version"
)

var cfgFile string

var Verbose bool
var PrintHeader bool
var NumProcs int

const (
	// Supported sequence file formats
	FastaFormat   string = "fasta"
	FastqFormat          = "fastq"
	UnknownFormat        = "unknown"
)

func GuessFileFormat(filename string) (format string, err error) {
	switch filepath.Ext(strings.ToLower(filename)) {
	case ".fastq", ".fq":
		return FastqFormat, nil
	case ".fasta", ".fa", ".fna", ".faa":
		return FastaFormat, nil
	}
	fmt.Fprintf(os.Stderr, "Unknown file format: %s\n", filepath.Ext(strings.ToLower(filename)))
	return UnknownFormat, nil
}

var MemProfileFileName string
var MemProfileFile *os.File
var CpuProfileFileName string
var CpuProfileFile *os.File

func StartProfiling() {
	if DEBUG && (MemProfileFileName != "" || CpuProfileFileName != "") {
		fmt.Fprintf(os.Stderr, "Starting Profiling.\nLogging CPU Profile to:  %s\nLogging Mem Profile to:  %s\n", CpuProfileFileName, MemProfileFileName)
	}
	if MemProfileFileName != "" {
		var err error
		MemProfileFile, err = os.Create(MemProfileFileName)
		if err != nil {
			log.Fatal(err)
		}
	}

	if CpuProfileFileName != "" {
		CpuProfileFile, err := os.Create(CpuProfileFileName)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(CpuProfileFile)
		//defer pprof.StopCPUProfile()
	}
}

func StopProfiling() {
	if DEBUG && (MemProfileFileName != "" || CpuProfileFileName != "") {
		fmt.Fprintf(os.Stderr, "Stopping Profiling.\n")
	}
	if MemProfileFileName != "" {
		pprof.WriteHeapProfile(MemProfileFile)
		MemProfileFile.Close()
	}
	if CpuProfileFileName != "" {
		pprof.StopCPUProfile()
		CpuProfileFile.Close()
	}
}

func main() {
}

// Execute adds all child commands to the root command sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

func init() {
	cobra.OnInitialize(initConfig)

	// Here you will define your flags and configuration settings.
	// Cobra supports Persistent Flags, which, if defined here,
	// will be global for your application.

	//RootCmd.PersistentFlags().StringVar(&cfgFile, "config", "", "config file (default is $HOME/.catseq.yaml)")
	RootCmd.PersistentFlags().BoolVarP(&Verbose, "verbose", "", false, "Enable verbose output.")
	RootCmd.PersistentFlags().BoolVarP(&PrintHeader, "print-header", "", false, "Include column header in output.")
	//RootCmd.PersistentFlags().IntVarP(&NumProcs, "procs", "p", 1, "Use up to this many processors/cores in parallel.")
	RootCmd.PersistentFlags().StringVarP(&MemProfileFileName, "memprofile", "", "", "Write a memory profile to this file.")
	RootCmd.PersistentFlags().StringVarP(&CpuProfileFileName, "cpuprofile", "", "", "Write a CPU profile to this file.")

	// Cobra also supports local flags, which will only run
	// when this action is called directly.
	RootCmd.Flags().BoolP("help", "h", false, "Show this help message.")

	flag.Usage = func() {
		switch GoseqCommand(flag.Arg(0)) {
		case INFO, FILTER:
			fmt.Fprintf(os.Stderr, "usage: catseq %v [command options...] [sequence files...]\n\n", flag.Arg(0))
			flag.PrintDefaults()
		case VERSION:
			versionCmd.Execute()
		default:
			fmt.Fprintf(os.Stderr, "usage: catseq [options] command [command options...] [sequence files...]\n\n")
			flag.PrintDefaults()
		}
	}

	runtime.GOMAXPROCS(NumProcs)

	if Verbose {
		fmt.Fprintf(os.Stderr, "verbose: %t\n", Verbose)
		fmt.Fprintf(os.Stderr, "using %d/%d available procs \n", NumProcs, runtime.NumCPU())
		fmt.Fprintf(os.Stderr, "trailing args: %s\n", flag.Args())
	}
}

// initConfig reads in config file and ENV variables if set.
func initConfig() {
	if cfgFile != "" { // enable ability to specify config file via flag
		viper.SetConfigFile(cfgFile)
	}

	viper.SetConfigName(".catseq") // name of config file (without extension)
	viper.AddConfigPath("$HOME")   // adding home directory as first search path
	viper.AutomaticEnv()           // read in environment variables that match

	// If a config file is found, read it in.
	if err := viper.ReadInConfig(); err == nil {
		fmt.Println("Using config file:", viper.ConfigFileUsed())
	}
}

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}
