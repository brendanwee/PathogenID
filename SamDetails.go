package main

import(
  "bufio"
  "strconv"
  "os"
  "strings"
)

func SamDetails(samFile string)string{
  prefix := MakeFolder("resources/SamPlots")
  reads, qualities, baseLocations := ReadSamFile(samFile)
  total := 0
  for i:=0;i<len(qualities);i++{
    total += qualities[i]
  }
  MakeHistogram(prefix, "Coverage", baseLocations)
  averageQuality := float64(total)/float64(len(qualities))
  MakeHistogram(prefix, "Quality of Aligned Reads", qualities)
  return "Reads Aligned: " +strconv.Itoa(reads)+ "&#13;&#10;" + "Average Quality For Aligned Reads: " + strconv.FormatFloat(averageQuality,'f',2,64)
}

func ReadSamFile(samFile string) (int, []int, []int){
  file,err := os.Open(samFile)
  qualities := make([]int,0)
  baseLocations := make([]int,0)
  reads := 0
  CheckError(err)
  scanner := bufio.NewScanner(file)
  scanner.Scan()
  for scanner.Scan(){
    reads += 1
    alignedRead := strings.Split(scanner.Text(),"\t")
    bases := make([]int,len(alignedRead[9]))
    pos, err := strconv.Atoi(alignedRead[9])
    CheckError(err)
    for i := range(bases){
      bases[i] = pos+i
    }
    baseLocations = append(baseLocations,bases...)
    qScore,err := strconv.Atoi(alignedRead[4])
    CheckError(err)
    qualities = append(qualities,qScore)
  }
  return reads, qualities, baseLocations
}
