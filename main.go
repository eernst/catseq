package main

import (
	"fmt"
	"github.com/eernst/catseq/cmd"
	"os"
)

func main() {
	err := cmd.RootCmd.Execute()
	if err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}
