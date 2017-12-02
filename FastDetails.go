import(
  "fmt"
  "bufio"
  "strings"
  "os"
)

type adapter [2]string //name then sequence

func MakeAdapterMap() map[adapter]int {
  adapterFile,err := os.Open("adapter.fasta")
  CheckError(err)
  scanner := bufio.NewScanner(adapterFile)
  adapterMap := make(map[adapter]int)
  for scanner.Scan(){
    if scanner.Text()!= "" && scanner.Text()[0] == '>'{
      var a adapter
      a[0] = scanner.Text()
      scanner.Scan()
      a[1] = scanner.Text()
      adapterMap[a]=0
    }
  }
  return adapterMap
}

func CountBases(G,A,C,T *int, sequence string){
  for i:=0; i<len(sequence);i++{
    switch sequence[i]{
      case 'G': *G++
      case 'g': *G++
      case 'C': *C++
      case 'c': *C++
      case 'T': *T++
      case 't': *T++
      case 'A': *A++
      case 'a': *A++
    }
  }
}

func SumQScore(sum *int, qscore string){
  for i:=0;i<len(qscore);i++{
    *sum += int(qscore[i])-32
  }
}

func CheckForAdapters(adapterMap map[adapter]int, sequence string){
  for adapter, count := range(adapterMap){
      if strings.Contains(sequence, adapter[1]){
        adapterMap[adapter] = count + 1
      }
  }
}

func PrintReadLengthsGraph(readLengths []int){
  /*categories := make(map[int]int)
  categories[50]=0
  categories[100]=0
  categories[150]=0
  for i:=0;i<23;i++{
    if n<8 {
      cat := 200+100*n
    } else if n<18 {
      cat := 1000+1000*(n-8)
    } else if n<20 {
      cat := 15000+5000*(n-18)
    } else {
      cat := 30000 + 10000*(n-20)
    }
    categories[cat] = 0
  }
  AssignCategory(categories, readLengths)
  prevKey := 0
  for key, value := range(categories){
    fmt.Println(prevKey,"-", key,"=", value)
    prevKey = key
  }*/

}
/*
func AssignCategory(categories map[int]int, readlengths){
  for i:=0; i<len(readLengths);i++{
    if len(line) < key :
      categories[key] = value + 1
      break
  }
}
*/
func FastDetails(readFiles... string) {
  MakeAdapterMap()
  adapterMap := MakeAdapterMap()
  var G int
  var C int
  var T int
  var A int
  var qualitySum int
  var numReads int
  var readLengths []int

  for i:=0; i<len(readFiles);i++{
    reads,err:= os.Open(readFiles[i])
    CheckError(err)
    scanner := bufio.NewScanner(reads)
    for scanner.Scan(){

      if scanner.Text()==""{
        continue
      } else if scanner.Text()[0]=='@' { //fastq
        numReads++
        scanner.Scan()
        readLengths = append(readLengths,len(scanner.Text()))
        CheckForAdapters(adapterMap, scanner.Text())
        CountBases(&G,&A,&C,&T, scanner.Text())
        scanner.Scan()//+
        scanner.Scan()//qscore
        SumQScore(&qualitySum, scanner.Text())
      } else if scanner.Text()[0]=='>'{
        numReads++
        scanner.Scan()
        readLengths = append(readLengths,len(scanner.Text()))
        CountBases(&G,&A,&C,&T, scanner.Text())
      }
    }

  }
  totalBases := G+C+A+T
  gContent := float64(G)/float64(totalBases)*100
  cContent := float64(C)/float64(totalBases)*100
  aContent := float64(A)/float64(totalBases)*100
  tContent := float64(T)/float64(totalBases)*100
  averageReadQuality := float64(qualitySum)/float64(totalBases)
  averageReadLengths := float64(totalBases)/float64(numReads)
  fmt.Println("G/C Content: G:", gContent, "C:", cContent)
  fmt.Println("Per Base Sequnce Content: G:", gContent, "C:", cContent, "A:", aContent, "T:", tContent)
  fmt.Println("Average Read Quality:", averageReadQuality)
  fmt.Println("Average Read Length:", averageReadLengths)
  fmt.Println( "Analyzed", len(readLengths), "reads")
  //PrintReadLengthsGraph(readLengths)
}
