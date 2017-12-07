package main

import "github.com/murlokswarm/app"

// FileSummary is the component displaying the summary of the input file.
type FileSummary struct {
	Output string
}

// Render returns the HTML describing the FileSummary component content.
// It contains a link to show how to navigate to an other component (Sf).
func (f *FileSummary) Render() string {
	return `
<div class="Home">
	<div class="Example">
		<h1>File Summary</h1>
		<textarea placeholder="Copy/Paste your input file. &#13;&#10; Or use Ctrl + O to open a file :)" oncontextmenu="OnContextMenu">
			{{if .Output}}
		  {{ .Output}}
			{{end}}
		</textarea>
    <button type = "button" onclick="OnFastqButtonClick">Show .fastq file</button>
    <button onclick="OnSAMButtonClick">Show .sam file</button>
    <button onclick="OnVariantButtonClick">Show .vcf file</button>
    <button onclick="OnClinicalButtonClick">Drug database</button>
    <button onclick="OnResistenceButtonClick">Drug-Resistence result</button>
    <button_go_back onclick="OnGoBackButtonClick">Go Back to Analysis</button_go_back>
	</div>
</div>
	`
}

// Graph is the component displaying the graph of the summary.
type Graph struct {
	Figure string
}

// Render returns the HTML describing the Graph component content.
// It contains a link to show how to navigate to an other component (Sf).
func (g *Graph) Render() string {
	return `
<div class="Home">
	<div class="Example">
		<h1>Base Content</h1>
		<figure>
  		<img
  			src="Fastq_Plots/BaseContent.png"
  			alt="An awesome picture"/>
		</figure>
		<button onclick = "CloseSubWin">Close</button>
	</div>
</div>
	`
}

// newSubWindow creates a new small window to display the plots
func newSubWindow() app.Contexter {
	// Creates a window context.
	subwin := app.NewWindow(app.Window{
		Title:          "Pathogen ID",
		Width:          640,
		Height:         360,
		TitlebarHidden: true,
		MinimizeHidden: false,
		OnClose: func() bool {
			subWin = nil
			return true
		},
	})
	graph := &Graph{}   // Creates a Hello component.
	subwin.Mount(graph) // Mounts the Hello component into the window context.
	//	graph.Figure = "profile-pictures/logo.png" // parse the source of the image
	//	app.Render(graph)
	return subwin
}

// OnFastqButtonClick is called when "show .fastq file" is clicked
func (f *FileSummary) OnFastqButtonClick() {
	if len(fastq) == 0 { // if there's no .fastq files
		f.Output = "There's no .fastq files. &#13;&#10; Please use Ctrl + O to open up .fastq files"
		app.Render(f) // Tells the app to update the rendering of the component.
	} else {
		f.Output = FastDetails(fastq) // Read in a slice of .fastq filename strings and return its summary as a string
		app.Render(f)                 // Tells the app to update the rendering of the component.
		subWin = newSubWindow()       // Pop up a new window for the graph display
	}
}

// OnSAMButtonClick is called when "show .sam file" is clicked
func (f *FileSummary) OnSAMButtonClick() {
	if len(sam) == 0 { // if there's no .fastq files
		f.Output = "There's no .sam files. &#13;&#10; Please use Ctrl + O to open up .sam files"
		app.Render(f) // Tells the app to update the rendering of the component.
	} else {
		f.Output = SamDetails(sam[0]) // Read in the filename of .sam file and return its summary as a tring
		app.Render(f)                 // Tells the app to update the rendering of the component.
		subWin = newSubWindow()       // Pop up a new window for the graph display
	}
}

// OnVariantButtonClick is called when "show .vcf file" is clicked
func (f *FileSummary) OnVariantButtonClick() {
	if len(vcf) == 0 { // if there's no .fastq files
		f.Output = "There's no .vcf files. &#13;&#10; Please use Ctrl + O to open up .vcf files"
		app.Render(f) // Tells the app to update the rendering of the component.
	} else {
		f.Output = "Fastq summay present here!  &#13;&#10; We can start a new line here too."
		app.Render(f)           // Tells the app to update the rendering of the component.
		subWin = newSubWindow() // Pop up a new window for the graph display
	}
}

// OnClinicalButtonClick is called when "Drug database" is clicked
func (f *FileSummary) OnClinicalButtonClick() {
	if len(db) == 0 { // if there's no .fastq files
		f.Output = "There's no database files. &#13;&#10; Please use Ctrl + O to open up database files"
		app.Render(f) // Tells the app to update the rendering of the component.
	} else {
		f.Output = "Fastq summay present here!  &#13;&#10; We can start a new line here too."
		app.Render(f)           // Tells the app to update the rendering of the component.
		subWin = newSubWindow() // Pop up a new window for the graph display
	}
}

// OnResistenceButtonClick is called when "Drug-Resistence result" is clicked
func (f *FileSummary) OnResistenceButtonClick() {
	if len(resistence) == 0 { // if there's no .fastq files
		f.Output = "There's no result files. &#13;&#10; Please use Ctrl + O to open up .txt result files"
		app.Render(f) // Tells the app to update the rendering of the component.
	} else {
		/*subwin := app.NewWindow(app.Window{
			Title:          "Pathogen ID",
			Width:          640,
			Height:         360,
			TitlebarHidden: true,
			MinimizeHidden: false,
			OnClose: func() bool {
				subWin = nil
				return true
			},
		})
		table := &FinalTable{}   // Creates a Hello component.
		subwin.Mount(table) // Mounts the Hello component into the window context.
		// graph.Figure = "profile-pictures/logo.png" // parse the source of the image

		app.Render(table)*/
		f.Output = "Fastq summay present here!  &#13;&#10; We can start a new line here too."
		app.Render(f)           // Tells the app to update the rendering of the component.
		subWin = newSubWindow() // Pop up a new window for the graph display
	}
}

// On OnGoBackButtonClick is called when "Go Back to analysis" is clicked
// It leads us to the previous interface of analysis
func (f *FileSummary) OnGoBackButtonClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
}

// /!\ Register the component. Required to use the component into a context.
func init() {
	app.RegisterComponent(&FileSummary{})
	app.RegisterComponent(&Graph{})
}
