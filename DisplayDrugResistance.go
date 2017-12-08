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
	Title    string
	RIF      bool
	rif1gene []string
	rif2gene []string
	rif3gene []string
	INH      bool
	inh1gene []string
	inh2gene []string
	inh3gene []string
	EMB      bool
	emb1gene []string
	emb2gene []string
	emb3gene []string
	PZA      bool
	pza1gene []string
	pza2gene []string
	pza3gene []string
	SM       bool
	sm1gene  []string
	sm2gene  []string
	sm3gene  []string
	AMI      bool
	ami1gene []string
	ami2gene []string
	ami3gene []string
	ETH      bool
	eth1gene []string
	eth2gene []string
	eth3gene []string
	FLQ      bool
	flq1gene []string
	flq2gene []string
	flq3gene []string
	PAS      bool
	pas1gene []string
	pas2gene []string
	pas3gene []string
}

// Render returns the HTML describing the FileSummary component content.
// It contains a link to show how to navigate to an other component (Sf).
func (p *FinalTable) Render() string {
	return `
	<!DOCTYPE html>
	<html>
	<head>
	<style>

	table, th, td {
	    border: 3px solid black;
	    border-collapse: collapse;
	}
	th {
		font-size: 30px;
	    color: rgb(205,184,158);
	    padding: 10px;
	    text-align: center;
	}
	td {
		padding: 10px;
	    text-align: left;
	}
	table tr:nth-child(2) {
	    background-color: rgb(233,0,0);
		color: black;
	}
	table tr:nth-child(odd) {
	   background-color:rgb(136,131,188);
	}
	table tr:nth-child(3) {
	    background-color: rgb(233,122,130);
		color: black;
	}

	</style>
	</head>
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
	  </tr>{{if .rif1gene}}
  <tr>
  	<td>RIF</td>
		<td>{{index .rif1gene 0}}</td>
    <td>{{index .rif1gene 1}}</td>
    <td>{{index .rif1gene 2}}</td>
    <td>{{index .rif1gene 3}}</td>
  </tr>{{end}}{{if .rif2gene}}
	<tr>
  	<td>RIF</td>
		<td>{{index .rif2gene 0}}</td>
    <td>{{index .rif2gene 1}}</td>
    <td>{{index .rif2gene 2}}</td>
    <td>{{index .rif2gene 3}}</td>
  </tr>{{end}}{{if .rif3gene}}
	<tr>
  	<td>RIF</td>
		<td>{{index .rif3gene 0}}</td>
    <td>{{index .rif3gene 1}}</td>
    <td>{{index .rif3gene 2}}</td>
    <td>{{index .rif3gene 3}}</td>
  </tr>{{end}}{{if .inh1gene}}
	<tr>
  	<td>INH</td>
		<td>{{index .inh1gene 0}}</td>
    <td>{{index .inh1gene 1}}</td>
    <td>{{index .inh1gene 2}}</td>
    <td>{{index .inh1gene 3}}</td>
  </tr>{{end}}{{if .inh2gene}}
	<tr>
  	<td>INH</td>
		<td>{{index .inh2gene 0}}</td>
    <td>{{index .inh2gene 1}}</td>
    <td>{{index .inh2gene 2}}</td>
    <td>{{index .inh2gene 3}}</td>
  </tr>{{end}}{{if .inh3gene}}
	<tr>
  	<td>INH</td>
		<td>{{index .inh3gene 0}}</td>
    <td>{{index .inh3gene 1}}</td>
    <td>{{index .inh3gene 2}}</td>
    <td>{{index .inh3gene 3}}</td>
  </tr>{{end}}{{if .emb1gene}}
	<tr>
  	<td>EMB</td>
		<td>{{index .emb1gene 0}}</td>
    <td>{{index .emb1gene 1}}</td>
    <td>{{index .emb1gene 2}}</td>
    <td>{{index .emb1gene 3}}</td>
  </tr>{{end}}{{if .emb2gene}}
	<tr>
  	<td>EMB</td>
		<td>{{index .emb2gene 0}}</td>
    <td>{{index .emb2gene 1}}</td>
    <td>{{index .emb2gene 2}}</td>
    <td>{{index .emb2gene 3}}</td>
  </tr>{{end}}{{if .emb3gene}}
	<tr>
  	<td>EMB</td>
		<td>{{index .emb3gene 0}}</td>
    <td>{{index .emb3gene 1}}</td>
    <td>{{index .emb3gene 2}}</td>
    <td>{{index .emb3gene 3}}</td>
  </tr>{{end}}{{if .pza1gene}}
	<tr>
  	<td>PZA</td>
		<td>{{index .pza1gene 0}}</td>
    <td>{{index .pza1gene 1}}</td>
    <td>{{index .pza1gene 2}}</td>
    <td>{{index .pza1gene 3}}</td>
  </tr>{{end}}{{if .pza2gene}}
	<tr>
  	<td>PZA</td>
		<td>{{index .pza2gene 0}}</td>
    <td>{{index .pza2gene 1}}</td>
    <td>{{index .pza2gene 2}}</td>
    <td>{{index .pza2gene 3}}</td>
  </tr>{{end}}{{if .pza3gene}}
	<tr>
  	<td>PZA</td>
		<td>{{index .pza3gene 0}}</td>
    <td>{{index .pza3gene 1}}</td>
    <td>{{index .pza3gene 2}}</td>
    <td>{{index .pza3gene 3}}</td>
  </tr>{{end}}{{if .sm1gene}}
	<tr>
  	<td>SM</td>
		<td>{{index .sm1gene 0}}</td>
    <td>{{index .sm1gene 1}}</td>
    <td>{{index .sm1gene 2}}</td>
    <td>{{index .sm1gene 3}}</td>
  </tr>{{end}}{{if .sm2gene}}
	<tr>
  	<td>SM</td>
		<td>{{index .sm2gene 0}}</td>
    <td>{{index .sm2gene 1}}</td>
    <td>{{index .sm2gene 2}}</td>
    <td>{{index .sm2gene 3}}</td>
  </tr>{{end}}{{if .sm3gene}}
	<tr>
  	<td>SM</td>
		<td>{{index .sm3gene 0}}</td>
    <td>{{index .sm3gene 1}}</td>
    <td>{{index .sm3gene 2}}</td>
    <td>{{index .sm3gene 3}}</td>
  </tr>{{end}}{{if .ami1gene}}
	<tr>
  	<td>AMI</td>
		<td>{{index .ami1gene 0}}</td>
    <td>{{index .ami1gene 1}}</td>
    <td>{{index .ami1gene 2}}</td>
    <td>{{index .ami1gene 3}}</td>
  </tr>{{end}}{{if .ami2gene}}
	<tr>
  	<td>AMI</td>
		<td>{{index .ami2gene 0}}</td>
    <td>{{index .ami2gene 1}}</td>
    <td>{{index .ami2gene 2}}</td>
    <td>{{index .ami2gene 3}}</td>
  </tr>{{end}}{{if .ami3gene}}
	<tr>
  	<td>AMI</td>
		<td>{{index .ami3gene 0}}</td>
    <td>{{index .ami3gene 1}}</td>
    <td>{{index .ami3gene 2}}</td>
    <td>{{index .ami3gene 3}}</td>
  </tr>{{end}}{{if .eth1gene}}
	<tr>
  	<td>ETH</td>
		<td>{{index .eth1gene 0}}</td>
    <td>{{index .eth1gene 1}}</td>
    <td>{{index .eth1gene 2}}</td>
    <td>{{index .eth1gene 3}}</td>
  </tr>{{end}}{{if .eth2gene}}
	<tr>
  	<td>ETH</td>
		<td>{{index .eth2gene 0}}</td>
    <td>{{index .eth2gene 1}}</td>
    <td>{{index .eth2gene 2}}</td>
    <td>{{index .eth2gene 3}}</td>
  </tr>{{end}}{{if .eth3gene}}
	<tr>
  	<td>ETH</td>
		<td>{{index .eth3gene 0}}</td>
    <td>{{index .eth3gene 1}}</td>
    <td>{{index .eth3gene 2}}</td>
    <td>{{index .eth3gene 3}}</td>
  </tr>{{end}}{{if .flq1gene}}
	<tr>
  	<td>FLQ</td>
		<td>{{index .flq1gene 0}}</td>
    <td>{{index .flq1gene 1}}</td>
    <td>{{index .flq1gene 2}}</td>
    <td>{{index .flq1gene 3}}</td>
  </tr>{{end}}{{if .flq2gene}}
	<tr>
  	<td>FLQ</td>
		<td>{{index .flq2gene 0}}</td>
    <td>{{index .flq2gene 1}}</td>
    <td>{{index .flq2gene 2}}</td>
    <td>{{index .flq2gene 3}</td>
  </tr>{{end}}{{if .flq3gene}}
	<tr>
  	<td>FLQ</td>
		<td>{{index .flq3gene 0}}</td>
    <td>{{index .flq3gene 1}}</td>
    <td>{{index .flq3gene 2}}</td>
    <td>{{index .flq3gene 3}}</td>
  </tr>{{end}}{{if .pas1gene}}
	<tr>
  	<td>PAS</td>
		<td>{{index .pas1gene 0}}</td>
    <td>{{index .pas1gene 1}}</td>
    <td>{{index .pas1gene 2}}</td>
    <td>{{index .pas1gene 3}}</td>
  </tr>{{end}}{{if .pas2gene}}
	<tr>
  	<td>PAS</td>
		<td>{{index .pas2gene 0}}</td>
    <td>{{index .pas2gene 1}}</td>
    <td>{{index .pas2gene 2}}</td>
    <td>{{index .pas2gene 3}}</td>
  </tr>{{end}}{{if .pas3gene}}
	<tr>
  	<td>PAS</td>
		<td>{{index .pas3gene 0}}</td>
    <td>{{index .pas3gene 1}}</td>
    <td>{{index .pas3gene 2}}</td>
    <td>{{index .pas3gene 3}}</td>
  </tr>{{end}}
</table>

</body>
</html>
`
}

func (p *FinalTable) DisplayFinalTable(resultFile string) {
	file, err := os.Open("/Users/jinkeliu/Documents/GitHub/PathogenID/Analysis/Results/DrugResistance.txt")
	CheckError(err)
	scanner := bufio.NewScanner(file)
	scanner.Scan()
	count := 0
	for scanner.Scan() {
		items := strings.Split(scanner.Text(), "\t")
		switch items[0] {
		case "Rifampin":
			if !p.RIF { //rifampin first time
				p.RIF = true
				count = 1
				p.rif1gene = append(items[:4], items[6])
			} else if count == 1 && p.RIF {
				p.rif2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.RIF {
				p.rif3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Isoniazid":
			if !p.INH { //INHampin first time
				p.INH = true
				count = 1
				p.inh1gene = append(items[:4], items[6])
			} else if count == 1 && p.INH {
				p.inh2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.INH {
				p.inh3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Pyrazinamide":
			if !p.PZA { //rifampin first time
				p.PZA = true
				count = 1
				p.pza1gene = append(items[:4], items[6])
			} else if count == 1 && p.PZA {
				p.pza2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.PZA {
				p.pza3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Streptomycin":
			if !p.SM { //rifampin first time
				p.SM = true
				count = 1
				p.sm1gene = append(items[:4], items[6])
			} else if count == 1 && p.SM {
				p.sm2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.SM {
				p.sm3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Aminoglycosides":
			if !p.AMI { //rifampin first time
				p.AMI = true
				count = 1
				p.ami1gene = append(items[:4], items[6])
			} else if count == 1 && p.AMI {
				p.ami2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.AMI {
				p.ami3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Ethambutol":
			if !p.EMB { //rifampin first time
				p.EMB = true
				count = 1
				p.emb1gene = append(items[:4], items[6])
			} else if count == 1 && p.EMB {
				p.emb2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.EMB {
				p.emb3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Ethionamide":
			if !p.ETH { //rifampin first time
				p.ETH = true
				count = 1
				p.eth1gene = append(items[:4], items[6])
			} else if count == 1 && p.ETH {
				p.eth2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.ETH {
				p.eth3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Fluoroquinolones":
			if !p.FLQ { //rifampin first time
				p.FLQ = true
				count = 1
				p.flq1gene = append(items[:4], items[6])
			} else if count == 1 && p.FLQ {
				p.flq2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.FLQ {
				p.flq3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		case "Para-Aminosalicylic Acid":
			if !p.PAS { //rifampin first time
				p.PAS = true
				count = 1
				p.pas1gene = append(items[:4], items[6])
			} else if count == 1 && p.PAS {
				p.pas2gene = append(items[:4], items[6])
				count += 1
			} else if count == 2 && p.PAS {
				p.pas3gene = append(items[:4], items[6])
				count += 1
			} else { //count is too high
				continue
			}
		}
	}
}

// /!\ Register the component. Required to use the component into a context.
func init() {
	app.RegisterComponent(&FinalTable{})
}
