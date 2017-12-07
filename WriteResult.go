package main

import (
	"fmt"
	"os"
	"strconv"
)

var filename string
MakeFolder("Analysis/Results")
func WriteResult(allDrug []Drug) string {
	filename = "/Analysis/Results/DrugResistance.txt"
	// Create a new empty sam file.
	resultFile, err := os.Create(filename)
	if err != nil {
		fmt.Println("Sorry: couldn’t create the file!")
		os.Exit(1)
	}
	defer resultFile.Close()
	// Write data to files.

	//fmt.Fprintln(resultFile, headerLine) // Print the order of the chain.
	// Write the result.
	var drugName string
	var geneName string
	var mutationLine string
	var description string
	for i := range allDrug {
		// Write the drug name
		drugName = allDrug[i].name
		fmt.Fprintln(resultFile, drugName)
		for j := range allDrug[i].resistance {
			// Write the resistant gene name
			geneName = allDrug[i].resistance[j].name
			fmt.Fprintln(resultFile, geneName)
			//Write the header line.
			headerLine := "Mutation Genes Position\tAmino Acid Changes\tCodon Changes\tConfidence\tdescription"
			fmt.Fprintln(resultFile, headerLine)
			for k := range allDrug[i].resistance[j].mutations {
				// write each mutation and its properties
				description = allDrug[i].resistance[j].mutations[k].refAA.name + "to" + allDrug[i].resistance[j].mutations[k].altAA.name
				mutationLine = strconv.Itoa(allDrug[i].resistance[j].mutations[k].Pos) + "\t" +
					allDrug[i].resistance[j].mutations[k].refAA.name + "->" + allDrug[i].resistance[j].mutations[k].altAA.name + "\t" +
					allDrug[i].resistance[j].mutations[k].refCodon + "->" + allDrug[i].resistance[j].mutations[k].altCodon + "\t" +
					strconv.FormatFloat(allDrug[i].resistance[j].mutations[k].confidence, 'E', -1, 32) + "\t" + description
				fmt.Fprintln(resultFile, mutationLine)
			}
		}
	}
  return cwd+"/Analysis/Results/DrugResistance.txt"
}
