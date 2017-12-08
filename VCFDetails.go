package main

import(
  "fmt"
  "gonum.org/v1/plot"
  "gonum.org/v1/plot/plotter"
  "gonum.org/v1/plot/plotutil"
  "gonum.org/v1/plot/vg"
)

/*
VCF details takes a path to a VCF file. If an anlysis has not been run on it already,
it predicts resistances based on it. Then for every gene a bar graph is made
and the number of mutations is drawn.
*/
func VCFDetails(vcfFile string){
  if len(allGenes)==0{
      reference,_ := HandleReference()
      PredictResistance(reference,vcfFile)
  }

  prefix := cwd+"/resources/VCFPlots/"
  //mutations by gene
  names := make([]string,len(allGenes))
  numMutations := make([]float64,len(allGenes))
  for gene := range(allGenes){
    names[gene] = allGenes[gene].name
    numMutations[gene] = float64(len(allGenes[gene].mutations))
  }
  fmt.Println(numMutations)
  plot := DrawBarGraph("Number of Mutations", "Mutations By Gene", numMutations)
  plot.NominalX(names...)
  plot.Save(10*vg.Inch, 10*vg.Inch, prefix+"MutationsByGene.png")
}

func DrawBarGraph(yAxis string, title string, data []float64) *plot.Plot{
  plot,err := plot.New()
  CheckError(err)
  plot.Title.Text = "Base Content"
	plot.Y.Label.Text = yAxis
  var dataPoints plotter.Values
  dataPoints = data
  bar,err := plotter.NewBarChart(dataPoints, vg.Points(400./float64(len(data))-4.))
  CheckError(err)
  bar.Offset= vg.Points(160./float64(len(data))-4.)
  bar.Color= plotutil.Color(1)
  plot.Add(bar)
  return plot
}
