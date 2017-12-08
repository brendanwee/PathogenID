package main

import (
	"github.com/murlokswarm/app"
)

// This is the initial interface when opening up the UI
type Hello struct {
	Greeting string
}

// Render returns the HTML markup that describes the appearance of the component.
// It supports standard HTML and extends it slightly to handle other component
// declaration or Golang callbacks.
func (h *Hello) Render() string {
	return `
<div class="WindowLayout">
    <div class="HelloBox">
        <h1>
            About tuberculosis
          <span>
                {{if .Greeting}}
                    {{html .Greeting}}
                {{else}}
                    pathogen
                {{end}}
            </span>
        </h1>
        <input type="text"
               value="{{html .Greeting}}"
               placeholder="Which pathogen are you interested in?"
               autofocus="true"
               onkeydown="Greeting"
               onkeyup="Greeting"
               onchange="OnInputChange" />
    </div>
		<button2 onclick="OnButtonClick">Next</button2>
</div>
    `
}

// OnInputChange is the handler called when an onchange event occurs.
// Once the users enter/select the pathogen of interest, move on to the analysis phase.
func (h *Hello) OnInputChange(arg app.ChangeArg) {
	pathogen = arg.Value // select the pathogen
	if Include(pathogen, PathogenList) {
		h.Greeting = arg.Value
		app.Render(h) // Tells the app to update the rendering of the component.
	} else {
		h.Greeting = "Sorry: The pathogen is not included in our local database."
		h.Greeting = "pathogen"
		app.Render(h)
	}
}

// Include returns true if what the user inputs is in our pathagen list
func Include(pathogen string, PathogenList []string) bool {
	for _, p := range PathogenList {
		if p == pathogen {
			return true
		}
	}
	return false
}

// OnButtonclick is the handler called when the "Next" button is clicked.
// It leads the user to the interface of analysis.
func (h *Hello) OnButtonClick() {
	analyze := &AnalyzeButton{} // creates a AnalyzeButton component.
	win.Mount(analyze)          // Mounts the AnalyzeButton component into the window context.
}

func init() {
	// Registers the Hello component.
	// Allows the app to create a Hello component when it finds its declaration into a HTML markup.
	app.RegisterComponent(&Hello{})
}
