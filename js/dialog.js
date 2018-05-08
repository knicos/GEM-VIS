function Dialog(element, cb) {
	let me = this;

	this.element = element; //document.createElement("div");
	this.element.className = "gem-vis dialog";
	this.model = null;
	this.csv = null;
	this.filename = null;
	this.reaction_select = null;
	this.subsystem_select = null;
	this.label_column_select = null;
	this.data_column_select = null;
	this.cb = cb;

	// Google Sheets API

	// Select reaction and metabolite columns

	// Choose model
	this.modelSelector();

	this.dataFile();
	this.attributeSelector();
	this.selectColumns();

	// Choose reactions or reaction sets.
	this.reactionSelector();

	// TODO Visualisation options
	// Don't show metabolites
	// Color for data, or node size, or both

	// Buttons
	this.buttons();
	
	document.body.appendChild(this.element);

	this.label_col = -1;
	this.data_col = -1;
	this.data = null;
	this.attribute = null;
	this.reaction_list = null;
}

Dialog.prototype.generateDataFromCSV = function() {
	if (!this.model || this.label_col < 0 || this.data_col < 0 || !this.csv) return;
	let isgene = this.attribute == "gene";

	let data = {};
	for (var i=1; i<this.csv.length; i++) {
		let l = this.csv[i][this.label_col];
		let v = parseFloat(this.csv[i][this.data_col]);

		//console.log("Gene", l);
		if (l.trim() == "") continue;

		if (isNaN(v)) console.error("NAN",this.csv[i][this.data_col]);

		if (isgene) {
			let r = this.model.genes[l];
			if (!r) continue;
			for (var j=0; j<r.length; j++) {
				data[r[j].id] = v;
			}
		} else {
			data[l] = v;
		}
	}

	console.log("DATA", data);
	this.data = data;
}

let models = {
	iSynCJ816: "https://knicos.github.io/SBML/samples/iSynCJ816-mmc12.xml",
	"iSynCJ816 (corrected)": "https://knicos.github.io/SBML/samples/iSynCJ816.xml"
};

function formEntry(label, element) {
	let f = document.createElement("div");
	let l = document.createElement("label");
	l.textContent = label;
	f.appendChild(l);
	f.appendChild(element);
	f.className = "question";
	return f;
}

Dialog.prototype.modelSelector = function() {
	let select = document.createElement("select");
	select.appendChild(document.createElement("option"));

	for (var x in models) {
		let opt = document.createElement("option");
		opt.textContent = x;
		opt.value = models[x];
		select.appendChild(opt);
	}
	this.element.appendChild(formEntry("Model",select));

	let me = this;
	select.onchange = function(e) {
		SBML.fromURL(select.value, function(model) {
			me.model = model;
			console.log(model);

			// Update reactions list
			me.updateReactions(model);
			// Identify allowed data key attributes
		});
	}
}

Dialog.prototype.attributeSelector = function() {
	let select = document.createElement("select");
	let opt = document.createElement("option");
	opt.textContent = "Reaction";
	opt.value = "reaction";
	select.appendChild(opt);

	opt = document.createElement("option");
	opt.textContent = "Metabolite";
	opt.value = "metabolite";
	select.appendChild(opt);

	opt = document.createElement("option");
	opt.textContent = "Gene";
	opt.value = "gene";
	select.appendChild(opt);

	let me = this;
	select.onchange = function(e) {
		me.attribute = e.target.value;
		me.generateDataFromCSV();
	}

	this.element.appendChild(formEntry("Attribute", select));
}

Dialog.prototype.reactionSelector = function() {
	let me = this;

	let select = document.createElement("select");
	select.setAttribute("multiple",true);
	this.reaction_select = select;
	this.element.appendChild(formEntry("Reactions",select));

	let subselect = document.createElement("select");
	this.subsystem_select = subselect;
	this.element.appendChild(formEntry("Subsystems",subselect));
	subselect.onchange = function(e) {
		if (e.target.value == "all") {
			me.reaction_list = me.model.reactions;
		} else if (e.target.value == "allenz") {
			me.reaction_list = me.model.reactions.filter(a => a.ec != "");
		} else if (e.target.value == "SKIP") {
			// Nothing
		} else if (e.target.value != "alldata") {
			me.reaction_list = me.model.subsystems[e.target.value].reactions;
		} else {
			var rl = [];
			if (me.data) {
				for (var x in me.data) {
					if (x === undefined || x.trim() == "") continue;
					rl.push(x);
				}
			}
			console.log(rl);
			me.reaction_list = rl;
		}
	}
}

Dialog.prototype.updateReactions = function(model) {
	while (this.reaction_select.lastChild) this.reaction_select.removeChild(this.reaction_select.lastChild);

	for (var i=0; i<model.reactions.length; i++) {
		let opt = document.createElement("option");
		opt.textContent = model.reactions[i].name;
		opt.value = model.reactions[i].id;
		this.reaction_select.appendChild(opt);
	}

	// Subsystems
	while (this.subsystem_select.lastChild) this.subsystem_select.removeChild(this.subsystem_select.lastChild);

	this.subsystem_select.appendChild(document.createElement("option"));
	let allopt = document.createElement("option");
	allopt.textContent = "All";
	allopt.value = "all";
	this.subsystem_select.appendChild(allopt);
	allopt = document.createElement("option");
	allopt.textContent = "All in data";
	allopt.value = "alldata";
	this.subsystem_select.appendChild(allopt);
	allopt = document.createElement("option");
	allopt.textContent = "All Enzymes";
	allopt.value = "allenz";
	this.subsystem_select.appendChild(allopt);

	allopt = document.createElement("option");
	allopt.textContent = "----------";
	allopt.value = "SKIP";
	this.subsystem_select.appendChild(allopt);

	for (var x in model.subsystems) {
		let opt = document.createElement("option");
		opt.textContent = x;
		opt.value = x;
		this.subsystem_select.appendChild(opt);
	}
}

Dialog.prototype.dataFile = function() {
	var me = this;

	let file = document.createElement("input");
	file.setAttribute("type", "file");
	this.element.appendChild(formEntry("Data",file));

	file.onchange = function(e) {
		let file = e.target.files[0];
		me.filename = file.name;
		
		let reader = new FileReader();
		reader.onload = function(e) {
			let csv = e.target.result;
			csv = csv.split("\n");
			for (var i=0; i<csv.length; i++) {
				let cols = csv[i].split(",");
				csv[i] = cols;
			}

			me.csv = csv;

			me.updateColumnsCSV(csv);
		}
		reader.readAsText(file);
	}
}

Dialog.prototype.selectColumns = function() {
	let me = this;

	let select = document.createElement("select");
	this.label_column_select = select;
	this.element.appendChild(formEntry("Label Column",select));
	select.onchange = function(e) {
		me.label_col = parseInt(e.target.value);
		me.generateDataFromCSV();
	}


	let dselect = document.createElement("select");
	this.data_column_select = dselect;
	this.element.appendChild(formEntry("Data Column",dselect));
	dselect.onchange = function(e) {
		console.log(e);
		me.data_col = parseInt(e.target.value);
		me.generateDataFromCSV();
	}
}

Dialog.prototype.updateColumnsCSV = function(csv) {
	//console.log("Update Columns", csv);

	while (this.label_column_select.lastChild) this.label_column_select.removeChild(this.label_column_select.lastChild);

	this.label_column_select.appendChild(document.createElement("option"));
	for (var i=0; i<csv[0].length; i++) {
		let opt = document.createElement("option");
		opt.textContent = csv[0][i];
		opt.value = i;
		this.label_column_select.appendChild(opt);
	}

	while (this.data_column_select.lastChild) this.data_column_select.removeChild(this.data_column_select.lastChild);

	this.data_column_select.appendChild(document.createElement("option"));
	for (var i=0; i<csv[0].length; i++) {
		let opt = document.createElement("option");
		if (i == 1) opt.setAttribute("default", true);
		opt.textContent = csv[0][i];
		opt.value = i;
		this.data_column_select.appendChild(opt);
	}
}

Dialog.prototype.buttons = function() {
	let me = this;	
	let go = document.createElement("button");
	go.textContent = "Go";
	go.onclick = function(e) {
		me.cb({
			model: me.model,
			reactionData: (me.attribute != "metabolite") ? me.data : null,
			metaboliteData: (me.attribute == "metabolite") ? me.data : null,
			reactions: me.reaction_list,
			skipMissing: false,
			hideMetabolites: false,
			hideMetaboliteNames: true,
			hideSpecials: true,
			showReactionNames: true,
			removeIsolatedMetabolites: true
		});
	}
	this.element.appendChild(go);
}

module.exports = Dialog;

