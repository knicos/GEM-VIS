function Menu(element) {
	// Populate the buttons
	this.element = element;

	// Initiate from history...
	this.reloadHistory();

	// Model
	this.button("object", ()=> {
		// Show model and reaction selection dialog
	});

	// Data
	this.button("database", ()=> {
		// Show data upload and setup dialog
	});

	// Options

	// Remove Node

	// Expand Node

	// Unfix Node
}

