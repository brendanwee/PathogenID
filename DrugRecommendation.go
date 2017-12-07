package main

import (
	"fmt"
)

var refCDNA string

// main
func DrugRecom() {
	InitializeAA()
	InitializeCodonTable()
	refCDNA = refSeq
	//Transcription(refSeq)
	//fmt.Println(aminoAcid)
	//fmt.Println(CodonTable)
	CheckMutation()
}

// Filter the final result
func FilterResult(allDrug []Drug) []Drug {
	allDrug = FilterDrugMutation(allDrug)
	allDrug = FilterDrug(allDrug)
	return allDrug
}

// Filter genes without mutations
func FilterDrugMutation(allDrug []Drug) []Drug {
	for i := len(allDrug) - 1; i >= 0; i-- {
		for j := len(allDrug[i].resistance) - 1; j >= 0; j-- {
			if allDrug[i].resistance[j].mutations == nil { // Delete this resistance gene
				allDrug[i].resistance = append(allDrug[i].resistance[0:j], allDrug[i].resistance[j+1:]...)
			}
		}
	}
	return allDrug
}

// Filter drug without resistant genes
func FilterDrug(allDrug []Drug) []Drug {
	for i := len(allDrug) - 1; i >= 0; i-- {
		if len(allDrug[i].resistance) == 0 {
			allDrug = append(allDrug[0:i], allDrug[i+1:]...)
			continue
		}
	}
	return allDrug
}

func CheckMutation() {
	for i := range allGenes {
		for j := range allGenes[i].mutations {
			if CheckResistance(&allGenes[i].mutations[j], allGenes[i]) {
				allGenes[i].mutations[j].resistance = true
			}
		}
	}
}

func CheckResistance(m *Mutation, gene Gene) bool {
	if CheckFrameShift(m) {
		return true
	} else if CheckSNPResis(m, gene) {
		return true
	}
	return false
}

// Check whether the mutation causes frame shift.
func CheckFrameShift(m *Mutation) bool {
	if len(m.refSeq) > 1 { // frame shift
		return true
	} else {
		for i := range m.altSeq {
			if len(m.altSeq[i]) > 1 {
				return true
			}
		}
		return false
	}
}

// Check whether the mutation causes SNP
func CheckSNPResis(m *Mutation, gene Gene) bool {
	refAA, altAA := FindAA(m, gene)
	if refAA.charge != altAA.charge {
		return true
	}
	if refAA.polarity != altAA.polarity {
		return true
	}
	return false
}

func FindAA(m *Mutation, gene Gene) (thisRefAA AA, thisAltAA AA) {
	position := (m.Pos - gene.boundary[0]) % 3
	var thisRefcodon, thisAltcodon string
	switch position {
	// first base in a codon
	case 0:
		thisRefcodon = refCDNA[m.Pos-1 : m.Pos+2]
		thisAltcodon = m.altSeq[0] + refCDNA[m.Pos:m.Pos+2]
	// second base in a codon
	case 1:
		thisRefcodon = refCDNA[m.Pos-2 : m.Pos+1]
		thisAltcodon = refCDNA[m.Pos-2:m.Pos-1] + m.altSeq[0] + refCDNA[m.Pos:m.Pos+1]
	// third base in a codon
	case 2:
		thisRefcodon = refCDNA[m.Pos-3 : m.Pos]
		thisAltcodon = refCDNA[m.Pos-3:m.Pos-1] + m.altSeq[0]
		fmt.Println(2, thisRefcodon, thisAltcodon)
	default:
		fmt.Println("An error happened!")
	}
	m.refCodon = thisRefcodon
	m.altCodon = thisAltcodon
	thisRefAA = CodonToAA(thisRefcodon)
	thisAltAA = CodonToAA(thisAltcodon)
	m.refAA = thisRefAA
	m.altAA = thisAltAA
	return thisRefAA, thisAltAA
}

// Find correspond amino acid
func CodonToAA(codon string) AA {
	return CodonTable[codon]
}

// // Transcrip the DNA to mRNA
// func Transcription(refSeq string) {
//   for i := range refSeq {
//     switch refSeq[i] {
//     case 'A':
//       refmRNA = refmRNA + "U"
//     case 'T':
//       refmRNA = refmRNA + "A"
//     case 'C':
//       refmRNA = refmRNA + "G"
//     case 'G':
//       refmRNA = refmRNA + "C"
//     }
//   }
//   fmt.Println(refmRNA)
// }
