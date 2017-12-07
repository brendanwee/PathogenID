package main

import (
  "fmt"
  "bufio"
  "os"
  "strings"
)

var refSeq string //32Mb...

// Read the DNA strand from file and store it into a slice.
func ReadFromFile(filename string) (patient []string, quality []string, seqName []string) {
  // Open the file.
  inFile, err := os.Open(filename)
  if err != nil {
		fmt.Println("Sorry: couldn't open the file!")
		os.Exit(1)
	}
	defer inFile.Close()
  // For every DNA stand line in the file,
	scanner := bufio.NewScanner(inFile)
  count := 0 // count the lines
	for scanner.Scan() {
		// append it to the  slice.
    count ++
    if count%4 == 2 { //sequence
      patient = append(patient, scanner.Text())
    }
    if count%4 == 0 { //quality
      quality = append(quality, scanner.Text())
    }
    if count%4 == 1 { //sequence name
      header := strings.Split(scanner.Text(), " ")
      seqName = append(seqName, header[0])
    }
	}
  return patient, quality, seqName
}

func ReadReference(filename string) string {
  // Open the file.
  inFile, err := os.Open(filename)
  if err != nil {
    fmt.Println("Sorry: couldn't open the file!")
    os.Exit(1)
  }
  defer inFile.Close()
  // For every DNA stand line in the file,
  scanner := bufio.NewScanner(inFile)
  scanner.Scan() // Skip the first line.
  for scanner.Scan() {
    // append it to the  slice.
    refSeq = refSeq + scanner.Text()
  }
  return refSeq
}
