package pipeline

import (
	"io"
	"log"
	"sync"

	"github.com/shenwei356/bio/seqio/fastx"
)

func ChannelRec(reader *fastx.Reader) <-chan *fastx.Record {
	out := make(chan *fastx.Record)
	go func() {
		for {
			record, err := reader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				check(err)
				break
			}
			out <- record.Clone()
		}
		close(out)
	}()
	return out
}

// Copied from https://blog.golang.org/pipelines
func MergeRec(chans ...<-chan *fastx.Record) chan *fastx.Record {
	var wg sync.WaitGroup
	out := make(chan *fastx.Record)
	output := func(c <-chan *fastx.Record) {
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

func check(e error) {
	if e != nil {
		log.Fatal(e)
	}
}
