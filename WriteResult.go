package main

import (
	"fmt"
	"os"
	"strconv"
)

var filename string
func WriteResult(allDrug []Drug) string {
	filename = cwd+"/Analysis/Results/DrugResistance.txt"

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
			if len(allDrug[i].resistance[j].mutations) != 0 {
				for k := range allDrug[i].resistance[j].mutations {
					if allDrug[i].resistance[j].mutations[k].resistance == true {
						geneName = allDrug[i].resistance[j].name
						fmt.Fprint(resultFile, geneName+"\t")
						// write each mutation and its properties
						description = allDrug[i].resistance[j].mutations[k].refAA.name + "(" + allDrug[i].resistance[j].mutations[k].refAA.polarity + ", " +
							strconv.Itoa(allDrug[i].resistance[j].mutations[k].refAA.charge) + ") " +
							"to " + allDrug[i].resistance[j].mutations[k].altAA.name + "(" + allDrug[i].resistance[j].mutations[k].altAA.polarity + ", " +
							strconv.Itoa(allDrug[i].resistance[j].mutations[k].altAA.charge) + ") "
						mutationLine = strconv.Itoa(allDrug[i].resistance[j].mutations[k].Pos) + "\t" +
							allDrug[i].resistance[j].mutations[k].refAA.name + "->" + allDrug[i].resistance[j].mutations[k].altAA.name + "\t" +
							allDrug[i].resistance[j].mutations[k].refCodon + "->" + allDrug[i].resistance[j].mutations[k].altCodon + "\t" +
							strconv.FormatFloat(allDrug[i].resistance[j].mutations[k].confidence, 'E', -1, 32) + "\t" + description
						fmt.Fprintln(resultFile, mutationLine)
					}
				}
			}
		}
	}
	//return cwd+"/Analysis/Results/DrugResistance.txt"
	return "Hi"
	}
