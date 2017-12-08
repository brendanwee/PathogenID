package main

import (
	"path/filepath"
	"runtime"

	"github.com/murlokswarm/app"
	_ "github.com/murlokswarm/mac"
)

// create global variables
var (
	fbpC, Rv0340, iniB, iniA, iniC, mabA, inhA, Rv1592c, Rv1772, ndh, katG, furA, srmR, fabD, kasA, accD6, oxyR, aphC, efpA, fadE24, nhoA Gene

	pncA, rpsL, rrs, gidB, rpoB, embB, tlyA, embR, Rv3124, Rv3125c, Rv3264c, Rv3266c, embC, embA, ethA, gyrB, gyrA, thyA Gene

	RIF, INH, PZA, SM, AMI, EMB, ETH, FLQ, PAS Drug
	allDrug                                    []Drug
	allMutant                                  MutationFile //not a file but ok

	cwd          string
	outputPath   string
	win          app.Contexter
	subWin       app.Windower
	fastq        []string
	sam          []string
	vcf          []string
	db           []string
	resistence   []string
	pathogen     string
	PathogenList = []string{"M. tuberculosis", "Staphylococcus aureus"} //A list of pathogens
	numProcs     int
	allGenes     []Gene
)

func main() {
	numProcs := runtime.NumCPU()
	if numProcs > 1 { //use all but one Processor just in case
		numProcs -= 1
	}
	cwd = pwd()                         //  the current work directory
	outputPath = MakeFolder("Analysis") // the directory of the result analysis folder
	MakeFolder("Analysis/Results")      // make the folder to store all the result files
	MakeFolder("resources/Fastq_Plots") // make the folder to store plots for fastq display
	MakeFolder("resources/SamPlots")    // make the folder to store plots for sam display
	MakeFolder("resources/VCFPlots")    // make the folder to store plots for vcf display

	PrepareBin() // creates command chmod 755 for each file in PathogenID/bin

	// OnLaunch is a handler which is called when the app is initialized and ready.
	// The main window is created here.
	app.OnLaunch = func() {
		appMenu := &AppMainMenu{}             // Creates the AppMainMenu component.
		if menuBar, ok := app.MenuBar(); ok { // Mounts the AppMainMenu component into the application menu bar.
			menuBar.Mount(appMenu)
		}

		appMenuDock := &DockMenu{}      // Creates the DockMenu component.
		if dock, ok := app.Dock(); ok { // Mounts the DockMenu component into the application dock.
			dock.Mount(appMenuDock)
			icon := filepath.Join(app.Resources(), "logo.png")
			dock.SetIcon(icon)
		}
		win = newMainWindow() // Create the main window.
		return
	}

	// OnReopen is a handler which is called when the app is reopened.
	// For example, when the dock icon is clicked.
	// The main window is created another time here.
	app.OnReopen = func() {
		if win != nil {
			return
		}
		win = newMainWindow() // Create the main window again.
		return
	}

	// OnTerminate is a handler which (if set) is called when the app is
	// requested to terminates. Return false cancels the termination request.
	app.OnTerminate = func() bool {
		if win != nil {
			return true
		} else {
			return false
		}
	}
	// Run the application
	app.Run()
}

// newMainWindow() creates a new application main window
func newMainWindow() app.Contexter {
	// Creates a window context.
	win := app.NewWindow(app.Window{
		Title:          "Pathogen ID",
		Width:          1280,
		Height:         720,
		TitlebarHidden: true,
		MinimizeHidden: false,
		OnClose: func() bool {
			win = nil
			return true
		},
	})
	hello := &Hello{} // Creates a Hello component as the welcome interface
	win.Mount(hello)  // Mounts the Hello component into the window context.
	return win
}
