package seqmath

import (
	"fmt"
	"math"
	"sort"
)

var errorProbForQ = make([]float64, 201)

func init() {
	// Pre-compute expensive P(e) up to Q=200
	for q := range errorProbForQ {
		errorProbForQ[q] = math.Pow(10, -1*float64(q)/10)
	}
}

func ErrorProbForQ(qual int) (prob float64) {
	if qual >= 0 && qual < 201 {
		return errorProbForQ[qual]

	} else {
		panic(fmt.Sprintf("Quality value %v out of bounds [0,200]", qual))
	}
}

// Nxx returns an int slice with all values N1..N50..N99 calculated for the input slice
// of sequence lengths. The input slice does not need to be sorted, but the total length
// must also be passed to avoid a second pass.
func Nxx(seqLens []int, totalSeqLength int) (nxx []int) {
	nxx = make([]int, 100)
	var sls = seqLens
	if !sort.IntsAreSorted(sls) {
		sort.Ints(sls)
	}
	var cumLen int = 0
	var n = 1
	for i := range sls {
		l := sls[len(sls)-1-i]
		cumLen += l
		for n < 100 && float64(cumLen) >= float64(n)*0.01*float64(totalSeqLength) {
			nxx[n] = l
			n++
		}
	}
	return nxx
}

func Median(seqLens []int) (median int) {
	var n = len(seqLens)

	switch {
	case n == 1:
		return seqLens[0]
	case n&1 == 0:
		return seqLens[n/2-1]
	default:
		return ((seqLens[n/2] + seqLens[n/2-1]) / 2)

	}
}
