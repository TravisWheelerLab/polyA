package polyA

import "math"

func ConfidenceOnly(region []int, lambdas []float64) []float64 {
	confidences := make([]float64, len(region))
	scoreTotal := 0.0

	for index := range region {
		power := float64(region[index]) * lambdas[index]
		convertedScore := math.Pow(2, power)
		confidences[index] = convertedScore
		scoreTotal += convertedScore
	}

	for index := range region {
		confidences[index] = confidences[index] / scoreTotal
	}

	return confidences
}
