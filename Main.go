package main

import (
	"path/filepath"
	"runtime"
	"github.com/murlokswarm/app"
	_ "github.com/murlokswarm/mac"
)

// create global variables
var (
	cwd          string
	outputPath   string
	win          app.Contexter
	subWin       app.Contexter
	fastq        []string
	sam          []string
	vcf          []string
	db           []string
	resistence   []string
	pathogen     string
	PathogenList = []string{"M. tuberculosis", "Staphylococcus aureus"} //A list of pathogens
	numProcs     int
)

func main() {
	numProcs := runtime.NumCPU()
	if numProcs > 1 { //use all but one Processor just in case
		numProcs -= 1
	}
	cwd = pwd()                       //  the current work directory
	outputPath = MakeFolder("Analysis") // the directory of the result analysis folder
	PrepareBin()

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

	// OnFinalize is a handler which (if set) is called when the app is about
	// to be terminated.
	// It should be used to perform any final cleanup before the application
	// terminates.
	app.OnFinalize = func() {

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
		OnMinimize:       func() {},
		OnDeminimize:     func() {},
		OnFullScreen:     func() {},
		OnExitFullScreen: func() {},
		OnMove:           func(x float64, y float64) {},
		OnResize:         func(width float64, height float64) {},
		OnFocus:          func() {},
		OnBlur:           func() {},
	})

	hello := &Hello{} // Creates a Hello component.
	win.Mount(hello)  // Mounts the Hello component into the window context.
	return win
}
