package main

import "github.com/murlokswarm/app"

// FileSummary is the component displaying the summary of the input file.
type AnalyzeButton struct {
	Info string
}

// Render returns the HTML describing the FileSummary component content.
// It contains a link to show how to navigate to an other component (Sf).
func (p *AnalyzeButton) Render() string {
	return `
<div class="Home">
	<div class="Example">
		<h1>Analysis menu</h1>
		<textarea  placeholder="Please load input files before analysis :)" >{{if .Info}}{{ .Info}}{{end}}</textarea>
		<button onclick="OnPipelineButtonClick">Run Pipeline</button>
    <button onclick="OnAlignmentButtonClick">Alignment</button>
    <button onclick="OnVariantButtonClick">Identify variant</button>
    <button onclick="OnClinicalButtonClick">Drug resistance</button>
    <button_previous onclick="OnPreviousButtonClick">Previous</button_previous>
    <button_next onclick="OnNextButtonClick">Next</button_next>
	</div>
</div>
	`
}

//PickFile is called when the users try to do analysis while there are no
//corresponding files loaded in app
func PickFile(h *AnalyzeButton, files *[]string, command string) {
	app.NewFilePicker(app.FilePicker{ // read the files using a file picker
		MultipleSelection: true,
		NoDir:             true,
		NoFile:            false,
		OnPick: func(filenames []string) {
			for i := range filenames { // read all the files into the filepaths array
				*files = append(*files, filenames[i])
			}
			// run different analysis according to their tags
			switch command {
			case "Pipeline":
				h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
				app.Render(h)
				RunPipeline(*files)
			case "Alignment":
				h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
				app.Render(h)
				samFile := OnlyAlign(*files) // alignment function; generate the sam file; return the filepath
				sam = append(sam, samFile)   // append the sam file path to the list
			case "Variant":
				h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
				app.Render(h)
				vcfFile := OnlyVCF(*files) // call VCF function; generate .vcf file; return the path
				vcf = append(vcf, vcfFile) // append the vcf file path to the list
			case "Drug":
				h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
				app.Render(h)
				reference, _ := HandleReference()
				txtFile := PredictResistance(reference, (*files)[0]) // function that analyzes the .vcf file
				resistence = append(resistence, txtFile)             // append txt file to the list for result file paths
			}
			fileSummary := &FileSummary{} // Creates a FileSummary component
			win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
		},
	})
}

// OnAlignmentButtonClick is called when when the button is clicked.
func (h *AnalyzeButton) OnPipelineButtonClick() {
	if len(fastq) == 0 { // there's no fastq file loaded
		h.Info = "Please input .fastq files by using Open file (Ctrl+O)"
		app.Render(h)
		PickFile(h, &fastq, "Pipeline")
	} else {
		h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
		app.Render(h)
		RunPipeline(fastq)
		fileSummary := &FileSummary{} // Creates a FileSummary component
		win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	}
}

// OnAlignmentButtonClick is called when when the button is clicked.
func (h *AnalyzeButton) OnAlignmentButtonClick() {
	if len(fastq) == 0 { // there's no fastq file loaded
		h.Info = "Please input .fastq files by using Open file (Ctrl+O)"
		app.Render(h)
		PickFile(h, &fastq, "Alignment")
	} else {
		h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
		app.Render(h)
		samFile := OnlyAlign(fastq)   // alignment function; generate the sam file; return the filepath
		sam = append(sam, samFile)    // append the sam file path to the list
		fileSummary := &FileSummary{} // Creates a FileSummary component
		win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	}
}

// OnAlignmentButtonClick is called when when the button is clicked.
func (h *AnalyzeButton) OnVariantButtonClick() {
	if len(sam) == 0 { // there's no sam file loaded
		h.Info = "Please input .sam files by using Open file (Ctrl+O)"
		app.Render(h)
		PickFile(h, &sam, "Variant")
	} else {
		h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
		app.Render(h)
		vcfFile := OnlyVCF(sam)       // call VCF function; generate .vcf file; return the path
		vcf = append(vcf, vcfFile)    // append the vcf file path to the list
		fileSummary := &FileSummary{} // Creates a FileSummary component
		win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	}
}

// OnAlignmentButtonClick is called when when the buttton is clicked.
func (h *AnalyzeButton) OnClinicalButtonClick() {
	if len(vcf) == 0 { // there's no vcf file loaded
		h.Info = "Please input .vcf files by using Open file (Ctrl+O)"
		app.Render(h)
		PickFile(h, &vcf, "Drug")
	} else {
		h.Info = "Running... Please wait patiently for jumping. This may take a few minutes"
		app.Render(h)
		reference, _ := HandleReference()
		txtFile := PredictResistance(reference, vcf[0]) // function that analyzes the .vcf file
		resistence = append(resistence, txtFile)        // append txt file to the list for result file paths
		fileSummary := &FileSummary{}                   // Creates a FileSummary component
		win.Mount(fileSummary)                          // Mounts the FileSummary component into the window context
	}
}

// OnPreviousButtonClick is called when the "Previous" button is clicked.
// It leads the users to the greeting interface.
func (a *AnalyzeButton) OnPreviousButtonClick() {
	hello := &Hello{}
	win.Mount(hello)
}

// OnNextButtonClick is called when the "Next" button is clicked.
// It leads the users to the display interface.
func (a *AnalyzeButton) OnNextButtonClick() {
	display := &FileSummary{}
	win.Mount(display)
}

// /!\ Register the component. Required to use the component into a context.
func init() {
	app.RegisterComponent(&AnalyzeButton{})
}
