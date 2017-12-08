package main

import (
	"fmt"
	"os"
	"strconv"
	"math"
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
	var geneName string
	var mutationLine string
	fmt.Fprintf(resultFile,"DRUG\tGENE\tPOS\tAA\tCOD\tCONF\tDESC\n")
	for i := range allDrug {
		fmt.Println(allDrug[i].name)
		for j := range allDrug[i].resistance { //for each gene associated with the drug
			for k := range allDrug[i].resistance[j].mutations { //for each mutation
				if allDrug[i].resistance[j].mutations[k].resistance == true {// if the mutation codes for resistance
					mutationLine = ""
					DRUG := GetFullName(allDrug[i].name)
					geneName = allDrug[i].resistance[j].name
					POS := strconv.Itoa(allDrug[i].resistance[j].mutations[k].Pos)
					AA := allDrug[i].resistance[j].mutations[k].refAA.name + " -> " + allDrug[i].resistance[j].mutations[k].altAA.name
					refPol :=allDrug[i].resistance[j].mutations[k].refAA.polarity
					altPol := allDrug[i].resistance[j].mutations[k].altAA.polarity
					refCharge :=ConvertCharge(allDrug[i].resistance[j].mutations[k].refAA.charge)
					altCharge :=ConvertCharge(allDrug[i].resistance[j].mutations[k].altAA.charge)
					COD := allDrug[i].resistance[j].mutations[k].refCodon + "->" + allDrug[i].resistance[j].mutations[k].altCodon + "\t"
					CONF := strconv.FormatFloat(math.Pow(10,(allDrug[i].resistance[j].mutations[k].confidence/10)),'f',2,64)
					DESC := "The change from a "+refPol+" "+refCharge+ " charged amino acid to a "+ altPol + " " + altCharge + " charged amino acid may affect this proteins susceptibility to "+DRUG+
					". We are "+CONF+"'%' sure this mutation is real."
					mutationLine = mutationLine + DRUG+"\t"+geneName+"\t"+POS+"\t"+AA+"\t"+COD+"\t"+CONF+"\t"+DESC
					fmt.Fprintf(resultFile,mutationLine+"\n")
				}
			}
		}
	}
	//return cwd+"/Analysis/Results/DrugResistance.txt"
	return "Hi"
}


func GetFullName(abv string) string {
	switch abv{
	case "RIF":
		return "Rifampin"
	case "INH":
		return "Isoniazid"
	case "PZA":
		return "Pyrazinamide"
	case "SM":
		return "Streptomycin"
	case "AMI":
		return "Aminoglycosides"
	case "EMB":
		return "Ethambutol"
	case "ETH":
		return "Ethionamide"
	case "FLQ":
		return "Fluoroquinolones"
	case "PAS":
		return "Para-Aminosalicylic Acid"
	}
	return ""
}

func ConvertCharge(i int)string{
	if i>0{
		return "positively"
	} else if i<0{
		return "negatively"
	} else {
		return "nuetrally"
	}
}
