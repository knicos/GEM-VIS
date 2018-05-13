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
	this.label_column_anno = null;
	this.data_column_anno = null;
	this.cb = cb;
	this.annotationsCSV = null;
	this.annotations = null;

	// Google Sheets API

	// Select reaction and metabolite columns

	// Choose model
	this.modelSelector();
	this.attributeSelector();
	// Choose reactions or reaction sets.
	this.reactionSelector();

	this.element.appendChild(document.createElement("hr"));

	this.dataFile();
	this.selectColumns();

	this.element.appendChild(document.createElement("hr"));

	this.annotationFile();
	this.annotationColumns();

	this.element.appendChild(document.createElement("hr"));

	this.optelements = [];
	this.addOptions();

	this.element.appendChild(document.createElement("hr"));

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

	this.element.appendChild(formEntry("Index", select));
}

// Return an array of the selected opion values
// select is an HTML select element
function getSelectValues(select) {
  var result = [];
  var options = select && select.options;
  var opt;

  for (var i=0, iLen=options.length; i<iLen; i++) {
    opt = options[i];

    if (opt.selected) {
      result.push(opt.value || opt.text);
    }
  }
  return result;
}

Dialog.prototype.reactionSelector = function() {
	let me = this;

	/*let select = document.createElement("select");
	select.setAttribute("multiple",true);
	this.reaction_select = select;
	this.element.appendChild(formEntry("Reactions",select));*/

	let subselect = document.createElement("select");
	this.subsystem_select = subselect;
	subselect.setAttribute("multiple",true);
	this.element.appendChild(formEntry("Subsystems",subselect));
	subselect.onchange = function(e) {
		let values = getSelectValues(e.target);
		me.reaction_list = [];

		console.log("SELECTED",values);

		for (var i=0; i<values.length; i++) {
			if (values[i] == "all") {
				me.reaction_list.push.apply(me.reaction_list,me.model.reactions);
			} else if (values[i] == "allenz") {
				me.reaction_list.push.apply(me.reaction_list,me.model.reactions.filter(a => a.ec != ""));
			} else if (values[i] == "SKIP") {
				// Nothing
			} else if (values[i] == "allanno") {
				var rl = [];
				if (me.annotations) {
					for (var x in me.annotations) {
						if (x === undefined || x.trim() == "") continue;
						rl.push(x);
					}
				}
				console.log(rl);
				me.reaction_list.push.apply(me.reaction_list,rl);
			} else if (values[i] != "alldata") {
				me.reaction_list.push.apply(me.reaction_list,me.model.subsystems[values[i]].reactions);
			} else {
				var rl = [];
				if (me.data) {
					for (var x in me.data) {
						if (x === undefined || x.trim() == "") continue;
						rl.push(x);
					}
				}
				console.log(rl);
				me.reaction_list.push.apply(me.reaction_list,rl);
			}
		}

		console.log("REACTIONS", me.reaction_list);
	}
}

Dialog.prototype.updateReactions = function(model) {
	/*while (this.reaction_select.lastChild) this.reaction_select.removeChild(this.reaction_select.lastChild);

	for (var i=0; i<model.reactions.length; i++) {
		let opt = document.createElement("option");
		opt.textContent = model.reactions[i].name;
		opt.value = model.reactions[i].id;
		this.reaction_select.appendChild(opt);
	}*/

	// Subsystems
	while (this.subsystem_select.lastChild) this.subsystem_select.removeChild(this.subsystem_select.lastChild);

	//this.subsystem_select.appendChild(document.createElement("option"));
	let allopt = document.createElement("option");
	allopt.textContent = "All";
	allopt.value = "all";
	this.subsystem_select.appendChild(allopt);
	allopt = document.createElement("option");
	allopt.textContent = "All in data";
	allopt.value = "alldata";
	this.subsystem_select.appendChild(allopt);
	allopt = document.createElement("option");
	allopt.textContent = "All in annotations";
	allopt.value = "allanno";
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

Dialog.prototype.annotationFile = function() {
	var me = this;

	let file = document.createElement("input");
	file.setAttribute("type", "file");
	this.element.appendChild(formEntry("Annotations",file));

	file.onchange = function(e) {
		let file = e.target.files[0];
		//me.filename = file.name;
		
		let reader = new FileReader();
		reader.onload = function(e) {
			let csv = e.target.result;
			csv = csv.split("\n");
			for (var i=0; i<csv.length; i++) {
				let cols = csv[i].split(",");
				csv[i] = cols;
			}

			me.annotationsCSV = csv;

			me.updateAnnoColumns(csv);
		}
		reader.readAsText(file);
	}
}

Dialog.prototype.annotationColumns = function() {
	let me = this;

	let select = document.createElement("select");
	let dselect = document.createElement("select");

	this.label_column_anno = select;
	this.element.appendChild(formEntry("Index Column",select));
	select.onchange = function(e) {
		//me.label_col = parseInt(e.target.value);
		//me.generateDataFromCSV();
		me.generateAnnotations(select.value, dselect.value);
	}

	this.data_column_anno = dselect;
	this.element.appendChild(formEntry("Text Column",dselect));
	dselect.onchange = function(e) {
		//console.log(e);
		//me.data_col = parseInt(e.target.value);
		//me.generateDataFromCSV();
		me.generateAnnotations(select.value, dselect.value);
	}
}

Dialog.prototype.generateAnnotations = function(index, text) {
	let data = {};
	let isgene = this.attribute == "gene";

	for (var i=1; i<this.annotationsCSV.length; i++) {
		let l = this.annotationsCSV[i][index].trim();
		let t = this.annotationsCSV[i][text];

		if (isgene) {
			let r = this.model.genes[l];
			if (!r) continue;
			for (var j=0; j<r.length; j++) {
				data[r[j].id] = t;
			}
		} else {
			data[l] = t;
		}
	}

	console.log("Annotations", data);
	this.annotations = data;
}

Dialog.prototype.selectColumns = function() {
	let me = this;

	let select = document.createElement("select");
	this.label_column_select = select;
	this.element.appendChild(formEntry("Index Column",select));
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

Dialog.prototype.updateAnnoColumns = function(csv) {
	//console.log("Update Columns", csv);

	while (this.label_column_anno.lastChild) this.label_column_anno.removeChild(this.label_column_anno.lastChild);

	this.label_column_anno.appendChild(document.createElement("option"));
	for (var i=0; i<csv[0].length; i++) {
		let opt = document.createElement("option");
		opt.textContent = csv[0][i];
		opt.value = i;
		this.label_column_anno.appendChild(opt);
	}

	while (this.data_column_anno.lastChild) this.data_column_anno.removeChild(this.data_column_anno.lastChild);

	this.data_column_anno.appendChild(document.createElement("option"));
	for (var i=0; i<csv[0].length; i++) {
		let opt = document.createElement("option");
		if (i == 1) opt.setAttribute("default", true);
		opt.textContent = csv[0][i];
		opt.value = i;
		this.data_column_anno.appendChild(opt);
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

Dialog.prototype.addOptions = function() {
	let outer = document.createElement("div");
	outer.className = "options-list";

	let opts = {
		"skipMissing": false,
		"hideMetaboliteNames": false,
		"hideSpecials": true,
		"showReactionNames": false,
		//"removeIsolated": true,
		//"subsystemCluster": false,
		"hideZero": false,
		"showAnnotationLabels": true,
		"clearCache": false,
		"noGravity": false,
		"multipleLinks": true,
		"linkThreshold": 1.0,
		"hideLoose": false
	};

	for (var x in opts) {
		let o = document.createElement("div");
		o.className = "question2";
		let lab = document.createElement("label");
		lab.setAttribute("for", "opt-"+x);
		lab.textContent = x;
		o.appendChild(lab);
		let opt = document.createElement("input");

		if (typeof opts[x] == "boolean") {
			opt.setAttribute("type", "checkbox");
			opt.checked = opts[x];
		} else {
			opt.setAttribute("type", "number");
			opt.value = opts[x];
		}
		opt.name = x;
		opt.id = "opt-"+x;
		o.appendChild(opt);
		outer.appendChild(o);
		this.optelements.push(opt);
	}

	//this.optelement = outer;
	this.element.appendChild(outer);
}

Dialog.prototype.saveToSVG = function() {
	let svge = document.querySelector("svg");
	var me = this;
	var svgString = new XMLSerializer().serializeToString(svge);

	//var canvas = document.createElement("canvas");
	//canvas.width = 3000;
	//canvas.height = 3000;
	//var ctx = canvas.getContext("2d");
	var DOMURL = self.URL || self.webkitURL || self;
	//var img = new Image();
	var svg = new Blob([svgString], {type: "image/svg+xml;charset=utf-8"});
	var url = DOMURL.createObjectURL(svg);

	/*img.onload = function() {
		ctx.drawImage(img, 0, 0);
		var png = canvas.toDataURL("image/png");
		//document.querySelector('#png-container').innerHTML = '<img src="'+png+'"/>';
		let dlink = (me.element.lastChild.nodeName == "A") ? me.element.lastChild : document.createElement("a");
		dlink.href = png;
		dlink.download = true;
		dlink.textContent = "Download";
		me.element.appendChild(dlink);
		DOMURL.revokeObjectURL(png);
	};
	img.src = url;*/

	let dlink = (me.element.lastChild.nodeName == "A") ? me.element.lastChild : document.createElement("a");
	dlink.href = url;
	dlink.download = "visualisation.svg";
	dlink.textContent = "Download SVG";
	me.element.appendChild(dlink);
}

Dialog.prototype.buttons = function() {
	let me = this;	
	let go = document.createElement("button");
	go.textContent = "Go";
	go.onclick = function(e) {
		let opts = {
			model: me.model,
			reactionData: (me.attribute != "metabolite") ? me.data : null,
			metaboliteData: (me.attribute == "metabolite") ? me.data : null,
			reactions: me.reaction_list,
			annotations: me.annotations,
			style: "metabolite"
		}

		for (var i=0; i<me.optelements.length; i++) {
			opts[me.optelements[i].name] = (me.optelements[i].getAttribute("type") == "checkbox") ? me.optelements[i].checked : me.optelements[i].value;
		}

		me.cb(opts);
	}
	this.element.appendChild(go);

	let save = document.createElement("button");
	save.textContent = "Save";
	save.onclick = function(e) {
		me.saveToSVG();
	}
	this.element.appendChild(save);
}

module.exports = Dialog;

