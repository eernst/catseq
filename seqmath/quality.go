package seqmath

import (
	"fmt"
	"math"
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
