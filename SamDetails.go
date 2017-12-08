package main

import(
  "bufio"
  "strconv"
  "os"
  "strings"
)

func SamDetails(samFile string)string{
  prefix := cwd+"/resources/SamPlots"
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
    reads += 1 // count the reads
    alignedRead := strings.Split(scanner.Text(),"\t")
    bases := make([]int,len(alignedRead[9])) // save each bases location on the genome to generate coverage plot
    pos, err := strconv.Atoi(alignedRead[3]) //the position in the genome
    CheckError(err)
    for i := range(bases){
      bases[i] = pos+i
    }
    baseLocations = append(baseLocations,bases...) //append to master slice
    qScore,err := strconv.Atoi(alignedRead[4]) //grab the qScore of the individual read
    CheckError(err)
    qualities = append(qualities,qScore) // append to master quality
  }
  file.Close()
  return reads, qualities, baseLocations
}
