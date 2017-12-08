package main

import (
	"bufio"
	"os"
	"strings"

	"github.com/murlokswarm/app"
	_ "github.com/murlokswarm/mac"
)

//TURN BACK! THIS CODE IS HIDEOUS!!!!! Quite possibly the most naive code I have
//ever written.

type FinalTable struct {
	OverallResistanceDescription string
	Title                        string
	Rif                          bool
	Rif1gene                     []string
	Rif2gene                     []string
	Rif3gene                     []string
	Inh                          bool
	Inh1gene                     []string
	Inh2gene                     []string
	Inh3gene                     []string
	Emb                          bool
	Emb1gene                     []string
	Emb2gene                     []string
	Emb3gene                     []string
	Pza                          bool
	Pza1gene                     []string
	Pza2gene                     []string
	Pza3gene                     []string
	Sm                           bool
	Sm1gene                      []string
	Sm2gene                      []string
	Sm3gene                      []string
	Ami                          bool
	Ami1gene                     []string
	Ami2gene                     []string
	Ami3gene                     []string
	Eth                          bool
	Eth1gene                     []string
	Eth2gene                     []string
	Eth3gene                     []string
	Flq                          bool
	Flq1gene                     []string
	Flq2gene                     []string
	Flq3gene                     []string
	Pas                          bool
	Pas1gene                     []string
	Pas2gene                     []string
	Pas3gene                     []string
}

// Render returns the HTML describing the FileSummary component content.
// It contains a link to show how to navigate to an other component (Sf).
func (p *FinalTable) Render() string {
	return `
	   <body style="background-color:rgb(45,52,145);">
	   <table style="width:100%">
	     <tr>
	       <th colspan="5">Drug Resistance Profile</th>
	     </tr>
	     <tr>
	     	<th colspan="5"> {{.OverallResistanceDescription}}</th>
	     </tr>
	     <tr>
	       <td>Drug</td>
	       <td>Gene</td>
	       <td>Amino Acid</td>
	       <td>Confidence</td>
	       <td>Description</td>
	     </tr>{{if eq (len .Rif1gene) 5}}
	     <tr>
	   		<td>{{index .Rif1gene 0}}</td>
	       <td>{{index .Rif1gene 1}}</td>
	       <td>{{index .Rif1gene 2}}</td>
	       <td>{{index .Rif1gene 3}}</td>
				 <td>{{index .Rif1gene 4}}</td>
	     </tr>{{end}}{{if eq (len .Rif2gene) 5}}
	   	<tr>
	     	<td>{{index .Rif2gene 0}}</td>
	   		<td>{{index .Rif2gene 1}}</td>
	       <td>{{index .Rif2gene 2}}</td>
	       <td>{{index .Rif2gene 3}}</td>
	       <td>{{index .Rif2gene 4}}</td>
	     </tr>{{end}}{{if eq (len .Rif3gene) 5}}
	   	<tr>
	   		<td>{{index .Rif3gene 0}}</td>
	       <td>{{index .Rif3gene 1}}</td>
	       <td>{{index .Rif3gene 2}}</td>
	       <td>{{index .Rif3gene 3}}</td>
				 <td>{{index .Rif3gene 4}}</td>
	     </tr>{{end}}{{if eq (len .Inh1gene) 5}}
	   	<tr>
	   		<td>{{index .Inh1gene 0}}</td>
	       <td>{{index .Inh1gene 1}}</td>
	       <td>{{index .Inh1gene 2}}</td>
	       <td>{{index .Inh1gene 3}}</td>
				 <td>{{index .Inh1gene 4}}</td>
	     </tr>{{end}}{{if eq (len .Inh2gene) 5}}
	   	<tr>
	   		<td>{{index .Inh2gene 0}}</td>
	       <td>{{index .Inh2gene 1}}</td>
	       <td>{{index .Inh2gene 2}}</td>
	       <td>{{index .Inh2gene 3}}</td>
				 <td>{{index .Inh2gene 4}}</td>
	     </tr>{{end}}{{if eq (len .Inh3gene) 5}}
	   	<tr>
	   		<td>{{index .Inh3gene 0}}</td>
	       <td>{{index .Inh3gene 1}}</td>
	       <td>{{index .Inh3gene 2}}</td>
	       <td>{{index .Inh3gene 3}}</td>
				 <td>{{index .Inh3gene 3}}</td>
	     </tr>{{end}}{{if eq (len .Emb1gene) 5}}
	   	<tr>
	   		<td>{{index .Emb1gene 0}}</td>
	       <td>{{index .Emb1gene 1}}</td>
	       <td>{{index .Emb1gene 2}}</td>
	       <td>{{index .Emb1gene 3}}</td>
				 <td>{{index .Emb1gene 4}}</td>
		  </tr>{{end}}{{if eq (len .Emb2gene) 5}}
	   	<tr>
	   		<td>{{index .Emb2gene 0}}</td>
	       <td>{{index .Emb2gene 1}}</td>
	       <td>{{index .Emb2gene 2}}</td>
	       <td>{{index .Emb2gene 3}}</td>
				 <td>{{index .Emb2gene 4}}</td>
		  </tr>{{end}}{{if eq (len .Emb3gene) 5}}
	   	<tr>
	   		<td>{{index .Emb3gene 0}}</td>
	       <td>{{index .Emb3gene 1}}</td>
	       <td>{{index .Emb3gene 2}}</td>
	       <td>{{index .Emb3gene 3}}</td>
				 <td>{{index .Emb3gene 4}}</td>
				</tr>{{end}}{{if eq (len .Pza1gene) 5}}
	   	<tr>
	   		<td>{{index .Pza1gene 0}}</td>
	       <td>{{index .Pza1gene 1}}</td>
	       <td>{{index .Pza1gene 2}}</td>
	       <td>{{index .Pza1gene 3}}</td>
				 <td>{{index .Pza1gene 4}}</td>
				</tr>{{end}}{{if eq (len .Pza2gene) 5}}
	   	<tr>
	   		<td>{{index .Pza2gene 0}}</td>
	       <td>{{index .Pza2gene 1}}</td>
	       <td>{{index .Pza2gene 2}}</td>
	       <td>{{index .Pza2gene 3}}</td>
				 <td>{{index .Pza2gene 4}}</td>
				</tr>{{end}}{{if eq (len .Pza3gene) 5}}
	   	<tr>
	   		<td>{{index .Pza3gene 0}}</td>
	       <td>{{index .Pza3gene 1}}</td>
	       <td>{{index .Pza3gene 2}}</td>
	       <td>{{index .Pza3gene 3}}</td>
				 <td>{{index .Pza4gene 4}}</td>
				</tr>{{end}}{{if eq (len .Sm1gene) 5}}
	   	<tr>
	   		<td>{{index .Sm1gene 0}}</td>
	       <td>{{index .Sm1gene 1}}</td>
	       <td>{{index .Sm1gene 2}}</td>
	       <td>{{index .Sm1gene 3}}</td>
				 <td>{{index .Sm1gene 4}}</td>
				</tr>{{end}}{{if eq (len .Sm2gene) 5}}
	   	<tr>
	   		<td>{{index .Sm2gene 0}}</td>
	       <td>{{index .Sm2gene 1}}</td>
	       <td>{{index .Sm2gene 2}}</td>
	       <td>{{index .Sm2gene 3}}</td>
				 <td>{{index .Sm2gene 4}}</td>
				</tr>{{end}}{{if eq (len .Sm3gene) 5}}
	   	<tr>
	   		<td>{{index .Sm3gene 0}}</td>
	       <td>{{index .Sm3gene 1}}</td>
	   		<td>{{index .Sm3gene 2}}</td>
	       <td>{{index .Sm3gene 3}}</td>
				 <td>{{index .Sm3gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Ami1gene) 5}}
		 <tr>
			 <td>{{index .Ami1gene 0}}</td>
				<td>{{index .Ami1gene 1}}</td>
				<td>{{index .Ami1gene 2}}</td>
				<td>{{index .Ami1gene 3}}</td>
				<td>{{index .Ami1gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Ami2gene) 5}}
		 <tr>
			 <td>{{index .Ami2gene 0}}</td>
				<td>{{index .Ami2gene 1}}</td>
				<td>{{index .Ami2gene 2}}</td>
				<td>{{index .Ami2gene 3}}</td>
				<td>{{index .Ami2gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Ami3gene) 5}}
		 <tr>
			 <td>{{index .Ami3gene 0}}</td>
				<td>{{index .Ami3gene 1}}</td>
				<td>{{index .Ami3gene 2}}</td>
				<td>{{index .Ami3gene 3}}</td>
				<td>{{index .Ami3gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Eth1gene) 5}}
		 <tr>
			 <td>{{index .Eth1gene 0}}</td>
				<td>{{index .Eth1gene 1}}</td>
				<td>{{index .Eth1gene 2}}</td>
				<td>{{index .Eth1gene 3}}</td>
				<td>{{index .Eth1gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Eth2gene) 5}}
		 <tr>
			 <td>{{index .Eth2gene 0}}</td>
				<td>{{index .Eth2gene 1}}</td>
				<td>{{index .Eth2gene 2}}</td>
				<td>{{index .Eth2gene 3}}</td>
				<td>{{index .Eth2gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Eth3gene) 5}}
		 <tr>
			 <td>{{index .Eth3gene 0}}</td>
				<td>{{index .Eth3gene 1}}</td>
				<td>{{index .Eth3gene 2}}</td>
				<td>{{index .Eth3gene 3}}</td>
				<td>{{index .Eth3gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Flq1gene) 5}}
		 <tr>
			 <td>{{index .Flq1gene 0}}</td>
				<td>{{index .Flq1gene 1}}</td>
				<td>{{index .Flq1gene 2}}</td>
				<td>{{index .Flq1gene 3}}</td>
				<td>{{index .Flq1gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Flq2gene) 5}}
		 <tr>
			 <td>{{index .Flq2gene 0}}</td>
				<td>{{index .Flq2gene 1}}</td>
				<td>{{index .Flq2gene 2}}</td>
				<td>{{index .Flq2gene 3}}</td>
				<td>{{index .Flq2gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Flq3gene) 5}}
		 <tr>
			 <td>{{index .Flq3gene 0}}</td>
				<td>{{index .Flq3gene 1}}</td>
				<td>{{index .Flq3gene 2}}</td>
				<td>{{index .Flq3gene 3}}</td>
				<td>{{index .Flq3gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Pas1gene) 5}}
		 <tr>
			 <td>{{index .Pas1gene 0}}</td>
				<td>{{index .Pas1gene 1}}</td>
				<td>{{index .Pas1gene 2}}</td>
				<td>{{index .Pas1gene 3}}</td>
				<td>{{index .Pas1gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Pas2gene) 5}}
		 <tr>
			 <td>{{index .Pas2gene 0}}</td>
				<td>{{index .Pas2gene 1}}</td>
				<td>{{index .Pas2gene 2}}</td>
				<td>{{index .Pas2gene 3}}</td>
				<td>{{index .Pas2gene 4}}</td>
			 </tr>{{end}}{{if eq (len .Pas3gene) 5}}
		 <tr>
			 <td>{{index .Pas3gene 0}}</td>
				<td>{{index .Pas3gene 1}}</td>
				<td>{{index .Pas3gene 2}}</td>
				<td>{{index .Pas3gene 3}}</td>
				<td>{{index .Pas3gene 4}}</td>
			 </tr>{{end}}
	   </table>
	   </body>
`
}

// HasString returns true if the array contains the string
func HasString(x []string, y string) bool {
	for i := range x {
		if x[i] == y {
			return true
		}
	}
	return false
}

func (p *FinalTable) DisplayFinalTable(resultFile string) {
	file, err := os.Open(resultFile)
	CheckError(err)
	scanner := bufio.NewScanner(file)
	scanner.Scan()
	count := 0
	var overall []string
	for scanner.Scan() {
		items := strings.Split(scanner.Text(), "\t")
		if !HasString(overall, items[0]) {
			overall = append(overall, items[0])
		}
		switch items[0] {
		case "Rifampin":
			if !p.Rif { //Rifampin first time
				p.Rif = true
				count = 1
				p.Rif1gene = append(items[:4], items[6])
			} else if count == 1 && p.Rif {
				p.Rif2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Rif {
				p.Rif3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Isoniazid":
			if !p.Inh { //Inhampin first time
				p.Inh = true
				count = 1
				p.Inh1gene = append(items[:4], items[6])
			} else if count == 1 && p.Inh {
				p.Inh2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Inh {
				p.Inh3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "PyrazinAmide":
			if !p.Pza { //Rifampin first time
				p.Pza = true
				count = 1
				p.Pza1gene = append(items[:4], items[6])
			} else if count == 1 && p.Pza {
				p.Pza2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Pza {
				p.Pza3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Streptomycin":
			if !p.Sm { //Rifampin first time
				p.Sm = true
				count = 1
				p.Sm1gene = append(items[:4], items[6])
			} else if count == 1 && p.Sm {
				p.Sm2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Sm {
				p.Sm3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Aminoglycosides":
			if !p.Ami { //Rifampin first time
				p.Ami = true
				count = 1
				p.Ami1gene = append(items[:4], items[6])
			} else if count == 1 && p.Ami {
				p.Ami2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Ami {
				p.Ami3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Ethambutol":
			if !p.Emb { //Rifampin first time
				p.Emb = true
				count = 1
				p.Emb1gene = append(items[:4], items[6])
			} else if count == 1 && p.Emb {
				p.Emb2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Emb {
				p.Emb3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "EthionAmide":
			if !p.Eth { //Rifampin first time
				p.Eth = true
				count = 1
				p.Eth1gene = append(items[:4], items[6])
			} else if count == 1 && p.Eth {
				p.Eth2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Eth {
				p.Eth3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Fluoroquinolones":
			if !p.Flq { //Rifampin first time
				p.Flq = true
				count = 1
				p.Flq1gene = append(items[:4], items[6])
			} else if count == 1 && p.Flq {
				p.Flq2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Flq {
				p.Flq3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Para-Aminosalicylic Acid":
			if !p.Pas { //Rifampin first time
				p.Pas = true
				count = 1
				p.Pas1gene = append(items[:4], items[6])
			} else if count == 1 && p.Pas {
				p.Pas2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.Pas {
				p.Pas3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		}
	}
	description := "Overall Resistance: "
	for i := range overall {
		description += overall[i] + "; "
	}
	p.OverallResistanceDescription = description
}

// /!\ Register the component. Required to use the component into a context.
func init() {
	app.RegisterComponent(&FinalTable{})
}
