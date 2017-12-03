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
		<menuitem label="Close" selector="performClose:" shortcut="meta+w" onclick="OnCloseFileMenuClick"/>
		<menuitem label="New File"  shortcut="shift+meta+n" onclick="OnNewFileMenuClick" separator="true"/>
		<menuitem label="Save"  shortcut="meta+s" onclick="OnSaveFileMenuClick"/>
		<menuitem label="Save As..." shortcut="shift+meta+s" onclick="OnSaveAsFileMenuClick"/>
</menu>
    `
}

func (m *FileMenu) OnOpenFileMenuClick() {
	app.NewFilePicker(app.FilePicker{
		MultipleSelection: true,
		NoDir:             true,
		NoFile:            false,
		OnPick: func(filenames []string) {
			// handle files here
			fileSummary := &FileSummary{} // Creates a FileSummary component
			win.Mount(fileSummary)        // Mounts the FileSummary component into the window context
		},
	})
}

func (m *FileMenu) OnCloseFileMenuClick() {

}

func (m *FileMenu) OnNewFileMenuClick() {

}

func (m *FileMenu) OnSaveFileMenuClick() {

}

func (m *FileMenu) OnSaveAsFileMenuClick() {

}

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

func (m *AnalyzeMenu) OnPipelineMenuClick() {

}

func (m *AnalyzeMenu) OnAlignmentMenuClick() {

}

func (m *AnalyzeMenu) OnVariantMenuClick() {

}
func (m *AnalyzeMenu) OnClinicalMenuClick() {

}

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
func (m *DisplayMenu) OnFastqMenuClick() {

}
func (m *DisplayMenu) OnSAMMenuClick() {

}
func (m *DisplayMenu) OnVCFMenuClick() {

}
func (m *DisplayMenu) OnDrugDBMenuClick() {

}
func (m *DisplayMenu) OnResistenceMenuClick() {

}

func init() {
	// Allows the app to create a AppMainMenu and WindowMenu components when it finds its declaration
	// into a HTML markup.
	app.RegisterComponent(&AppMainMenu{})
	app.RegisterComponent(&WindowMenu{})
	app.RegisterComponent(&FileMenu{})
	app.RegisterComponent(&AnalyzeMenu{})
	app.RegisterComponent(&DisplayMenu{})
}
