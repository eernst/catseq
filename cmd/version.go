package cmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

const (
	MAJOR    = 1
	MINOR    = 0
	REVISION = 0
)

func init() {
	RootCmd.AddCommand(versionCmd)
}

var timeLayout string // the layout for time.Time

var (
	commitHash string
	buildDate  string
)

var versionCmd = &cobra.Command{
	Use:   "version",
	Short: "Print the version number.",
	Long:  `Output the version number of this binary. What more can be said?`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Printf("goseq version %v.%v.%v\n\n", MAJOR, MINOR, REVISION)
	},
}
