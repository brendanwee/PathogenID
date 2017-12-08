package main

import (
	"github.com/murlokswarm/app"
	"github.com/murlokswarm/log"
)

// AppMainMenu implements app.Componer interface.
type AppMainMenu struct {
	CustomTitle string
	Disabled    bool
}

// Render returns the HTML markup that describes the appearance of the component.
// In this case, the component will be mounted into a menu context.
// This restrict the markup to a compositon of menu and menuitem.
func (m *AppMainMenu) Render() string {
	return `
<menu>
    <menu label="app">
        <menuitem label="{{if .CustomTitle}}{{.CustomTitle}}{{else}}Custom item{{end}}"
                  shortcut="ctrl+c"
                  onclick="OnCustomMenuClick"
                  icon="star.png"
                  separator="true"
                  disabled="{{.Disabled}}" />
        <menuitem label="Quit" shortcut="meta+q" selector="terminate:" />
    </menu>
    <WindowMenu />
		<FileMenu />
		<AnalyzeMenu />
		<DisplayMenu />
		<menu label="Help"></menu>
</menu>
    `
}

// OnCustomMenuClick is the handler called when an onclick event occurs in a menuitem.
func (m *AppMainMenu) OnCustomMenuClick() {
	log.Info("OnCustomMenuClick")
}

// WindowMenu implements app.Componer interface.
// It's another component which will be nested inside the AppMenu component.
type WindowMenu struct {
}

func (m *WindowMenu) Render() string {
	return `
<menu label="Window">
		<menuitem label="New Window" onclick="OnNewWindowMenuClick" shortcut="meta+n"/>
		<menuitem label="Close" selector="performClose:" shortcut="meta+w" separator="true"/>
		<menuitem label="Minimize" selector="performMiniaturize:" shortcut="meta+m" />
		<menuitem label="Zoom" selector="performZoom:" separator="true" />
		<menuitem label="Bring All to Front" selector="arrangeInFront:" />
</menu>
    `
}

// OnNewWindowMenuClick is called when "New Window" is clicked
func (m *WindowMenu) OnNewWindowMenuClick() {
	win = newMainWindow()
}

// FileMenu implements app.Componer interface.
// It's another component which will be nested inside the AppMenu component.
type FileMenu struct {
}

func (m *FileMenu) Render() string {
	return `
<menu label="File">
    <menuitem label="Open..." shortcut="meta+o" onclick="OnOpenFileMenuClick"/>
		<menuitem label="Close" shortcut="meta+w" onclick="OnCloseFileMenuClick"/>
		<menuitem label="New File"  shortcut="shift+meta+n" onclick="OnNewFileMenuClick" separator="true"/>
		<menuitem label="Save"  shortcut="meta+s" onclick="OnSaveFileMenuClick"/>
		<menuitem label="Save As..." shortcut="shift+meta+s" onclick="OnSaveAsFileMenuClick"/>
</menu>
    `
}

/*
// This struct is for a textarea that our users can copy paste
// their data and save them somewhere
type NewFile struct {
	Content string
}

func (n *NewFile) Render() string {
	return `
<div class="Home">
	<div class="Example">
		<textarea placeholder = "Copy/Paste your new file here">{{if .Content}}{{. Content}}{{end}}</textarea>
	</div>
</div>
	`
}
*/

// OnOpenFileMenuClick is called when "Open..." is clicked
// It recognizes the suffix of the files and load their absolute paths
func (m *FileMenu) OnOpenFileMenuClick() {
	// Our users may visualize the data once they are loaded into the app
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	app.NewFilePicker(app.FilePicker{
		MultipleSelection: true,
		NoDir:             true,
		NoFile:            false,
		OnPick: func(filenames []string) {
			// add files to file path arrays
			for i := range filenames {
				length := len(filenames[i])
				if filenames[i][length-5:length] == "fastq" {
					// The file chosen is a .fastq file
					fastq = append(fastq, filenames[i])
				} else if filenames[i][length-3:length] == "sam" {
					// The file chosen is a .sam file
					sam = append(sam, filenames[i])
				} else if filenames[i][length-3:length] == "vcf" {
					// The file chosen is a .vcf file
					vcf = append(vcf, filenames[i])
				} else if filenames[i][length-3:length] == "txt" {
					// The file chosen is a .txt file
					resistence = append(resistence, filenames[i])
				} else { // The file chosen has invalid suffix
					fileSummary.Output = "The file chosen is invalid. We take only .fastq, .sam, or .vcf files"
					app.Render(fileSummary)
				}
			}
		},
	})
}

// OnCloseFileMenuClick is called when the "Close" is clicked
func (m *FileMenu) OnCloseFileMenuClick() {
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	fastq = fastq[:0]             // clear fastq files
	sam = sam[:0]                 // clear sam files
	vcf = vcf[:0]                 // clear vcf files
	resistence = resistence[:0]   // clear txt files
	fileSummary.Output = "All files are closed!"
	app.Render(fileSummary)
}

// OnNewFileMenuClick is called when "New File" is clicked
func (m *FileMenu) OnNewFileMenuClick() {
}

// OnSaveFileMenuClick is called when "Save" is clicked
func (m *FileMenu) OnSaveFileMenuClick() {
}

// OnSaveAsFileMenuClick is called when "Save As" is clicked
func (m *FileMenu) OnSaveAsFileMenuClick() {
}

// AnalyzeMenu implements app.Componer interface.
// It's another component which will be nested inside the AppMenu component.
type AnalyzeMenu struct {
}

func (m *AnalyzeMenu) Render() string {
	return `
<menu label="Analyze">
		<menuitem label="Run Pipeline" onclick="OnPipelineMenuClick" separator="true"/>
	  <menuitem label="Alignment" onclick="OnAlignmentMenuClick"/>
		<menuitem label="Identify variants"  onclick="OnVariantMenuClick" separator="true"/>
		<menuitem label="Clinical analysis"  onclick="OnClinicalMenuClick"/>
</menu>
    `
}

// OnPipelineMenuClick is clicked when "Run Pipeline" is clicked
// The same as OnPipelineButtonClick is clicked
func (m *AnalyzeMenu) OnPipelineMenuClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
	analyze.OnPipelineButtonClick()
}

// OnAlignmentMenuClick is clicked when "Alignment" is clicked
// The same as OnAlignmentButtonClick is clicked
func (m *AnalyzeMenu) OnAlignmentMenuClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
	analyze.OnAlignmentButtonClick()
}

// OnVariantMenuClick is clicked when "Identify variants" is clicked
// The same as OnVariantButtonClick is clicked
func (m *AnalyzeMenu) OnVariantMenuClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
	analyze.OnVariantButtonClick()
}

// OnClinicalMenuClick is clicked when "Clinical analysis" is clicked
// The same as OnClinicalButtonClick is clicked
func (m *AnalyzeMenu) OnClinicalMenuClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
	analyze.OnClinicalButtonClick()
}

// DisplayMenu implements app.Componer interface.
// It's another component which will be nested inside the AppMenu component.
type DisplayMenu struct {
}

func (m *DisplayMenu) Render() string {
	return `
<menu label="Display">
    <menuitem label="Fastq file"  onclick="OnFastqMenuClick"/>
		<menuitem label="SAM file"  onclick="OnSAMMenuClick"/>
		<menuitem label="VCF file"   onclick="OnVCFMenuClick" separator="true"/>
		<menuitem label="Drug database"   onclick="OnDrugDBMenuClick"/>
		<menuitem label="Resistence"  onclick="OnResistenceMenuClick"/>
</menu>
    `
}

// OnFastqMenuClick is called when the "Fastq" button is clicked.
// The same as OnFastqButtonClick is clicked
func (m *DisplayMenu) OnFastqMenuClick() {
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	fileSummary.OnFastqButtonClick()
}

// OnSAMMenuClick is called when the "SAM" button is clicked.
// The same as OnSAMButtonClick is clicked
func (m *DisplayMenu) OnSAMMenuClick() {
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	fileSummary.OnSAMButtonClick()
}

// OnVCFMenuClick is called when the "VCF" button is clicked.
// The same as OnVariantButtonClick is clicked
func (m *DisplayMenu) OnVCFMenuClick() {
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	fileSummary.OnVariantButtonClick()
}

// OnDrugDBMenuClick is called when the "Drug database" button is clicked.
// The same as OnClinicalButtonClick is clicked
func (m *DisplayMenu) OnDrugDBMenuClick() {
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	fileSummary.OnClinicalButtonClick()
}

// OnResistenceMenuClick is called when the "Resistence" button is clicked.
// The same as OnResistenceButtonClick is clicked
func (m *DisplayMenu) OnResistenceMenuClick() {
	fileSummary := &FileSummary{} // Creates a FileSummary component
	win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
	fileSummary.OnResistenceButtonClick()
}

func init() {
	// Allows the app to create a AppMainMenu and WindowMenu components when it finds its declaration
	// into a HTML markup.
	app.RegisterComponent(&AppMainMenu{})
	app.RegisterComponent(&WindowMenu{})
	app.RegisterComponent(&FileMenu{})
	app.RegisterComponent(&AnalyzeMenu{})
	app.RegisterComponent(&DisplayMenu{})
	//	app.RegisterComponent(&NewFile{})

}
