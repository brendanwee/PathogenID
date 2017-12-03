package main

import "github.com/murlokswarm/app"

// FileSummary is the component displaying the summary of the input file.
type AnalyzeButton struct{}

// Render returns the HTML describing the FileSummary component content.
// It contains a link to show how to navigate to an other component (Sf).
func (p *AnalyzeButton) Render() string {
	return `
<div class="Home">
	<div class="Example">
		<h1>Analysis menu</h1>
		<button onclick="OnPipelineButtonClick">Run Pipeline</button>
    <button onclick="OnAlignmentButtonClick">Alignment</button>
    <button onclick="OnVariantButtonClick">Identify variant</button>
    <button onclick="OnClinicalButtonClick">Drug resistence</button>
    <button_previous onclick="OnPreviousButtonClick">Previous</button_previous>
    <button_next onclick="OnNextButtonClick">Next</button_next>
	</div>
</div>
	`
}

func PickFile(files []string) {
	app.NewFilePicker(app.FilePicker{ // read the files using a file picker
		MultipleSelection: true,
		NoDir:             true,
		NoFile:            false,
		OnPick: func(filenames []string) {
			for i := range filenames { // read all the files into the filepaths array
				files = append(files, filenames[i])
			}
			fileSummary := &FileSummary{} // Creates a FileSummary component
			win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
		},
	})
}

// OnAlignmentButtonClick is called when when the buttton is clicked.
func (h *AnalyzeButton) OnPipelineButtonClick() {
	if len(fastq) == 0 { // there's no fastq file loaded
		PickFile(fastq)
	}
	RunPipeline(fastq)
}

// OnAlignmentButtonClick is called when when the buttton is clicked.
func (h *AnalyzeButton) OnAlignmentButtonClick() {
	if len(fastq) == 0 { // there's no fastq file loaded
		PickFile(fastq)
	}
	samFile := OnlyAlign(fastq) // alignment function; generate the sam file; return the filepath
	sam = append(sam, samFile)  // append the sam file path to the list

}

// OnAlignmentButtonClick is called when when the buttton is clicked.
func (h *AnalyzeButton) OnVariantButtonClick() {
	if len(sam) == 0 { // there's no sam file loaded
		PickFile(sam)
	}
	/*
		vcfFile := OnlyVCF(sam)    // call VCF function; generate .vcf file; return the path
		vcf = append(vcf, vcfFile) // append the vcf file path to the list
	*/
}

// OnAlignmentButtonClick is called when when the buttton is clicked.
func (h *AnalyzeButton) OnClinicalButtonClick() {

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
