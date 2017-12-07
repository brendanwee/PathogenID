package main

import(
  "fmt"
  "bufio"
  "strings"
  "os"
  "gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
  "gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
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

func PrintReadLengthsGraph(prefix string, readLengths []int){
  lengths := make(plotter.Values, len(readLengths))
	for i :=0; i< len(readLengths); i++ {
		lengths[i] = float64(readLengths[i])
	}
  graph,err := plot.New()
  CheckError(err)
  graph.Title.Text= "Read Lengths"
  histogram, err := plotter.NewHist(lengths,16)
  CheckError(err)
  graph.Add(histogram)
  graph.Save(4*vg.Inch, 4*vg.Inch, prefix+"ReadLengths.png")

}

func DrawBaseContentGraph(gContent, cContent, tContent, aContent float64, prefix string){
  plot,err := plot.New()
  CheckError(err)
  plot.Title.Text = "Base Content"
	plot.Y.Label.Text = "Percent"
  width := vg.Points(40)
  gData := plotter.Values{gContent}
  cData := plotter.Values{cContent}
  tData := plotter.Values{tContent}
  aData := plotter.Values{aContent}
  gBar, err := plotter.NewBarChart(gData, width)
  CheckError(err)
  cBar, err := plotter.NewBarChart(cData, width)
  CheckError(err)
  tBar, err := plotter.NewBarChart(tData, width)
  CheckError(err)
  aBar,err := plotter.NewBarChart(aData, width)
  CheckError(err)
  gBar.Offset= -1.5*width-4
  cBar.Offset= -.5*width-2
  tBar.Offset= .5*width+2
  aBar.Offset= 1.5*width+4
  gBar.Color= plotutil.Color(0)
  cBar.Color= plotutil.Color(1)
  tBar.Color= plotutil.Color(2)
  aBar.Color= plotutil.Color(3)
  plot.Add(gBar,cBar,tBar, aBar)
  plot.Legend.Add("G", gBar)
  plot.Legend.Add("C", cBar)
  plot.Legend.Add("T", tBar)
  plot.Legend.Add("A", aBar)
  plot.Legend.Top = true
  plot.Save(4*vg.Inch, 4*vg.Inch, prefix+"BaseContent.png")
}

func DrawAdapterContent(adapterMap map[adapter]int, prefix string){
  plot,err := plot.New()
  CheckError(err)
  plot.Title.Text = "Adapter Content"
	plot.Y.Label.Text = "Reads"
  width := vg.Points(15)
  i :=0
  for adpt, count := range(adapterMap){
    numReads := plotter.Values{float64(count)}
    adapterBar, err := plotter.NewBarChart(numReads, width)
    CheckError(err)
    offset :=vg.Points(-(23./2.)+.5+float64(i))
    adapterBar.Offset= offset*width
    adapterBar.Color= plotutil.Color(i)
    plot.Add(adapterBar)
    plot.Legend.Add(adpt[0], adapterBar)
    i++
  }
  plot.Save(8*vg.Inch, 4*vg.Inch, prefix+"AdapterContent.png")

}


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
  PrintReadLengthsGraph(strings.Split(readFiles[0],".")[0],readLengths)
  DrawBaseContentGraph(gContent,cContent,tContent,aContent,strings.Split(readFiles[0],".")[0])
  DrawAdapterContent(adapterMap,strings.Split(readFiles[0],".")[0])
}
