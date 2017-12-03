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
		<textarea  placeholder="Copy/Paste your input file or use meta + o to open a file :)" oncontextmenu="OnContextMenu">{{if .Output}}{{ .Output}}{{end}}</textarea>
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
type Graph struct{}

// Render returns the HTML describing the Graph component content.
// It contains a link to show how to navigate to an other component (Sf).
func (g *Graph) Render() string {
	return `
<div class="Home">
	<div class="Example">
		<h1>Visualize the file</h1>
		<figure>
  		<img
  			src="profile-pictures/logo.png"
  			alt="An awesome picture"/>
		</figure>
	</div>
</div>
	`
}

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
		OnMinimize:       func() {},
		OnDeminimize:     func() {},
		OnFullScreen:     func() {},
		OnExitFullScreen: func() {},
		OnMove:           func(x float64, y float64) {},
		OnResize:         func(width float64, height float64) {},
		OnFocus:          func() {},
		OnBlur:           func() {},
	})
	graph := &Graph{}   // Creates a Hello component.
	subwin.Mount(graph) // Mounts the Hello component into the window context.
	return subwin
}

func (f *FileSummary) OnFastqButtonClick() {
	f.Output = "Fastq summay present here!"
	app.Render(f)           // Tells the app to update the rendering of the component.
	subWin = newSubWindow() // Pop up a new window for the graph display
}

func (f *FileSummary) OnSAMButtonClick() {

}

func (f *FileSummary) OnVariantButtonClick() {

}
func (f *FileSummary) OnClinicalButtonClick() {

}
func (f *FileSummary) OnResistenceButtonClick() {

}

func (f *FileSummary) OnGoBackButtonClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
}

// /!\ Register the component. Required to use the component into a context.
func init() {
	app.RegisterComponent(&FileSummary{})
	app.RegisterComponent(&Graph{})
}
