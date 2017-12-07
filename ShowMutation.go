package main

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
	"math"
	"math/rand"
)

//store mutations as a struct
type Mutation struct {
	refSeq     string
	altSeq     []string
	Pos        int
	confidence float64
	resistance bool // whether the mutation causes drug resistance
	refAA AA
	altAA AA
	refCodon string
	altCodon string
}

//store gene as a struct contains all the mutations about it
type Gene struct {
	name        string
	boundary    [2]int
	orientation bool
	mutations   []Mutation
	drug        []string
}

func (gene Gene) PrintMutation() {
	for i := range gene.mutations {
		fmt.Println("Mutation:", gene.mutations[i])
	}
}

//store all the mutations in a file
type MutationFile []Mutation

func (all MutationFile) Len() int { // rewrite Len() method
	return len(all)
}
func (all MutationFile) Swap(i, j int) { // rewrite Swap() method
	all[i], all[j] = all[j], all[i]
}
func (all MutationFile) Less(i, j int) bool { // rewrite Less() method, sort from large to small.
	return all[j].Pos < all[i].Pos
}

//store all the drugs
type Drug struct {
	name       string
	resistance []Gene
}


//main
func ShowMutation(vcfName string) {
	vcfLines := ReadVcf(vcfName)
	allMutant = SplitVcf(vcfLines)
	InitializeGene()
	InitializeDrug()
	allMutant = FilterMutaions(allMutant)
	AssignMutant(allGenes, allMutant)
	UpdateDrug(allDrug, allGenes)
	//PrintAll()
	allDrug = FilterResult(allDrug)
}

// Delete those mutations whose confidence is relatively low.
func FilterMutaions(allMutant MutationFile) MutationFile {
	for i := len(allMutant)-1; i >=0; i-- {
		confPCT := ConfTransfer(allMutant[i].confidence)
	 	// Further check whether we keep this mutation
		if !FurtherCheck(confPCT) { // Delete this mutation
			allMutant = append(allMutant[0:i], allMutant[i+1:]...)
		}
	}
	return allMutant
}

//  Further check whether keep the mutation with low confidence
func FurtherCheck(confPCT float64) bool {
	checkNum := float64(rand.Intn(100))/100.0
	if checkNum <= confPCT {
		return true
	} else {
		return false
	}
}

// Transfer the confidence
func ConfTransfer(confidence float64) float64 {
	return 1.0 - math.Pow(10, confidence/(-10.0))
}

//Print drugs and resistances
func PrintAll() {
	for i := range allDrug {
		fmt.Println("Drug Name:", allDrug[i].name)
		for j := range allDrug[i].resistance {
			fmt.Println("Gene Name:", allDrug[i].resistance[j].name)
			allDrug[i].resistance[j].PrintMutation()
		}
	}
}

func UpdateDrug(allDrug []Drug, allGenes []Gene) {
	for i := range allDrug {
		for j := range allGenes {
			for k := range allGenes[j].drug {
				if allDrug[i].name == allGenes[j].drug[k] {
					allDrug[i].resistance = append(allDrug[i].resistance, allGenes[j])
					break
				}
			}
		}

	}
}

//Assign mutant to Genes
func AssignMutant(allGenes []Gene, allMutant MutationFile) {
	for i := range allGenes {
		allGenes[i] = GenesMutation(allGenes[i], allMutant)
	}
}

//read the input file line by line
func ReadVcf(vcfName string) []string {
	// open the file and make sure all went well
	file, err := os.Open(vcfName)
	if err != nil {
		fmt.Printf("An error occurred on opening the inputfile\n" +
			"Does the file exist?\n" +
			"Have you got access to it?\n")
		os.Exit(1)
	}
	// create the variable to hold the lines
	var vcfLines []string = make([]string, 0)
	// for every line in the file
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// append it to the lines slice
		vcfLines = append(vcfLines, scanner.Text())
	}

	// check that all went ok
	if scanner.Err() != nil {
		fmt.Println("Sorry: there was some kind of error during the file reading")
		os.Exit(1)
	}
	// close the file and return the lines
	defer file.Close()
	return vcfLines
}

//spilt the lines using \t and trun into Mutation
func SplitVcf(vcfLines []string) MutationFile {

	for i := range vcfLines {
		if vcfLines[i][0] == '#' {
			continue
		} else {
			//ignore chrom because the bacteria only have 1
			mutant := strings.Split(vcfLines[i], "\t")
			//maybe add this when filter is ready
			/*if mutant[6] != "PASS" {
			  continue*/

			if mutant[4] != "." {
				var thisMutant Mutation
				thisMutant.Pos, _ = strconv.Atoi(mutant[1])
				thisMutant.refSeq = mutant[3]
				thisMutant.altSeq = append(thisMutant.altSeq, mutant[4])
				thisMutant.confidence, _ = strconv.ParseFloat(mutant[5], 64)
				allMutant = append(allMutant, thisMutant)
			}

		}
	}
	sort.Sort(sort.Reverse(MutationFile(allMutant)))
	return allMutant
}

func GenesMutation(gene Gene, mutant MutationFile) Gene {
	//Search gene location
	startMutant, movedRight := SearchMutant(mutant, gene.boundary)
	if startMutant>=len(mutant){
		return gene
	}
	if movedRight { //keep moving right, or larger
		for i := startMutant; mutant[i].Pos > gene.boundary[1]; i++{
			gene.mutations = append(gene.mutations, mutant[i])
		}
	} else {
		for i := startMutant; mutant[i].Pos > gene.boundary[0]; i-- {
			gene.mutations = append(gene.mutations, mutant[i])
		}
	}
	return gene
}

func inBounds(i int, boundary [2]int) int{
	if i>boundary[1]{
		return 2
	} else if i<boundary[0]{
		return 0
	} else {
		return 1
	}
}

func MoveRightToBoundary(mutant MutationFile, i int, boundary [2]int) int{
	for mutant[i].Pos < boundary[1]{
		i += 1
	}
	return i-1
}

//maybe exist problem. NOT ANYMORE! >:()
func SearchMutant(mutant MutationFile, boundary [2]int) (int,bool) {
	now := 0
	low, high := 0, len(mutant)-1
	first := true
	var movedRight bool
	for low <= high {
		now = (low + high) >> 1
		side := inBounds(mutant[now].Pos,boundary)
		switch side {
			case 0: //left
				low = now+1
				movedRight =true
			case 1: //inside!
				if first{
					now = MoveRightToBoundary(mutant,now, boundary)
					movedRight = true
				}
				return now, movedRight
			case 2:
				high = now-1
				movedRight = false
			default:
				fmt.Println("ERROR: somehow, i is not greater, larger, or equal to the boundaries")
		}
		first = false
	}
	return now, movedRight
}
