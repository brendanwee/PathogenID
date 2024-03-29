package main

import(
  "fmt"
  "bufio"
  "strings"
  "strconv"
  "os"
  "gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
  "gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

type adapter [2]string //name then sequence

/*
MakeAdapterMap makes a map where the keys are the sequences of the adapters
in adapter.fasta and the keys are ints
*/
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
  adapterFile.Close()
  return adapterMap
}

/*
counts all the nucleotide bases in a sequence of bytes
*/
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

/*
adds the q score for each base to the same and appends it to the q slice as well. the q slice is then returned.
*/
func SumQScore(sum *int, qscore string) []int{
  var q []int
  for i:=0;i<len(qscore);i++{
    q = append(q,int(qscore[i])-32)
    *sum += int(qscore[i])-32
  }
  return q
}

/*
for every sequence in the adapter map, a string is queried to determine if the adapter falls within it.
*/
func CheckForAdapters(adapterMap map[adapter]int, sequence string){
  for adapter, count := range(adapterMap){
      if strings.Contains(sequence, adapter[1]){
        adapterMap[adapter] = count + 1
      }
  }
}
/*
Given a path, filename, and some data, MakeHistogram draws a histogram of the data
and saves is at path/filename
*/
func MakeHistogram(prefix, title string, data []int){
  lengths := make(plotter.Values, len(data))
	for i :=0; i< len(data); i++ {
		lengths[i] = float64(data[i])
	}
  graph,err := plot.New()
  CheckError(err)
  graph.Title.Text= title
  histogram, err := plotter.NewHist(lengths,16)
  CheckError(err)
  graph.Add(histogram)
  fmt.Println("Saved ", "/"+prefix+"/"+title+".png")
  graph.Save(4*vg.Inch, 4*vg.Inch, "/"+prefix+"/"+title+".png")
}

/*
Draws a bar graph of a the bases
*/
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

/*
draws of bar chart where each adapter is given a bar in the graph.
*/
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


/*
FastDetails iterates over read files and reads them line by line. If its a
fastq its sequence line and quality line are mined for data. If just a fasta then the sequence line is mined
*/
func FastDetails(readFiles []string) string{
  prefix := cwd+"/resources/Fastq_Plots/"
  MakeAdapterMap()
  adapterMap := MakeAdapterMap()
  var G int
  var C int
  var T int
  var A int
  var qualitySum int
  var numReads int
  var readLengths []int
  var qScores []int

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

        qScores = append(qScores,SumQScore(&qualitySum, scanner.Text())...)
      } else if scanner.Text()[0]=='>'{
        numReads++
        scanner.Scan()
        readLengths = append(readLengths,len(scanner.Text()))
        CountBases(&G,&A,&C,&T, scanner.Text())
      }
    }
    reads.Close()
  }
  totalBases := G+C+A+T
  gContent:= float64(G)/float64(totalBases)*100
  cContent:= float64(C)/float64(totalBases)*100
  aContent:= float64(A)/float64(totalBases)*100
  tContent:= float64(T)/float64(totalBases)*100
  averageReadQuality := float64(qualitySum)/float64(totalBases)
  averageReadLengths := float64(totalBases)/float64(numReads)
  output := `G/C Content: G: `+ strconv.FormatFloat(gContent,'f',2,64) + ` C: ` + strconv.FormatFloat(cContent,'f',2,64)+ `&#13;&#10;`+
`Per Base Sequnce Content: G: ` + strconv.FormatFloat(gContent,'f',2,64) + ` C: ` + strconv.FormatFloat(cContent,'f',2,64) + ` A: ` + strconv.FormatFloat(aContent,'f',2,64) + ` T: ` + strconv.FormatFloat(tContent,'f',2,64) + `&#13;&#10;` +
`Average Read Quality: ` + strconv.FormatFloat(averageReadQuality,'f',2,64) +`&#13;$#10;` +
`Average Read Length: ` + strconv.FormatFloat(averageReadLengths,'f',2,64) + `&#13;$#10;` +
`Analyzed ` + strconv.Itoa(numReads) + `reads`
  MakeHistogram(prefix,"readlengths",readLengths)
  MakeHistogram(prefix,"Quality Scores of Raw Reads", qScores)
  DrawBaseContentGraph(gContent,cContent,tContent,aContent,prefix)
  DrawAdapterContent(adapterMap,prefix)
  return output
}
