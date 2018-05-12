//const ForceGraph3D = require('3d-force-graph');
//const metanal = require('../metabolite-analysis.js');
//const saveSvgAsPng = require("save-svg-as-png");
const biocyc = require('./compounds.js');
const wrap = require('word-wrap');

function processMetabolite(m) {
	for (var x in biocyc) {
		let b = biocyc[x];
		if (b && b.name && m.name.toLowerCase() == b.name.toLowerCase()) {
			return b;
		}

		if (b.synonyms) {
			for (var i=0; i<b.synonyms.length; i++) {
				if (b.synonyms[i].toLowerCase() == m.name.toLowerCase()) {
					//console.log("SYNO MATCH",b);
					return b;
				}
			}
		}
	}

	//console.log("NO MATCH",m);
}

function MetabolicGraph(parent, model, options) {
	this.model = model;
	this.options = (options) ? options : {};
	this.force = null;

	// TODO Process metabolites with BioCyc Data
	for (var i=0; i<model.metabolites.length; i++) {
		model.metabolites[i].biocyc = processMetabolite(model.metabolites[i]);
	}

	// Remove all in parent
	while (parent.lastChild) parent.removeChild(parent.lastChild);

	this.element = document.createElement("div");
	this.element.className = "synechocystis-graph";

	var width = 3000,
    height = 3000;
	this.svg = d3.select(this.element).append("svg")
		.attr("width", width)
		.attr("height", height);

	if (parent) parent.appendChild(this.element);

	parent.scrollLeft = 1000;
	parent.scrollTop = 1000;

	let button = document.createElement("button");
	var me = this;
	button.textContent = "Save";
	button.onclick = function() {
		saveSvgAsPng(me.svg, "model.png");
	}
	this.element.appendChild(button);

	if (!this.options.clearCache && window.localStorage && window.localStorage.gemvis_cached) {
		this.cached = JSON.parse(window.localStorage.gemvis_cached);
	} else {
		this.cached = {};
	}
}

// TODO Generate intelligently
const specialMetabolites = {
	"M_atp_c": true,
	"M_nad_c": true,
	"M_h_c": true,
	"M_pi_c": true,
	"M_adp_c": true,
	"M_h2o_c": true,
	"M_o2_c": true,
	"M_nadph_c": true,
	"M_nadp_c": true,
	"M_co2_c": true,
	"M_nadh_c": true,
	"M_hco3_c": true,

	"M_atp_x": true,
	"M_nad_x": true,
	"M_h_x": true,
	"M_pi_x": true,
	"M_adp_x": true,
	"M_h2o_x": true,
	"M_o2_x": true,
	"M_nadph_x": true,
	"M_co2_x": true,
	"M_nadh_x": true,
	"M_hco3_x": true,

	"M_atp_p": true,
	"M_nad_p": true,
	"M_h_p": true,
	"M_pi_p": true,
	"M_adp_p": true,
	"M_h2o_p": true,
	"M_o2_p": true,
	"M_nadph_p": true,
	"M_co2_p": true,
	"M_nadh_p": true,
	"M_hco3_p": true,

	"M_h2o_l": true,
	"M_h2o_u": true
};

let excluded = {
	"General": true,
	"Biomass": true,
	"Others": true,
	"Exchange reactions": true,
	"Modeling": true,
	"Transport": true
};

function calculateEdge(data1, data2, m, mdata) {
	let t1 = (Math.pow(Math.abs(data1)+1,1)-1) / m.producers.length;
	let t2 = (Math.pow(Math.abs(data2)+1,1)-1) / m.consumers.length;
	let sig = (t1) >= 0.01*mdata; //*mdata;
	let value = (data1 / m.consumers.length) * (data2*0.5 + ((data2 < 0) ? -1 : 1)); // / m.producers.length

	return [value, sig];
}

function rotate(cx, cy, x, y, angle) {
    var radians = (Math.PI / 180) * angle,
        cos = Math.cos(radians),
        sin = Math.sin(radians),
        nx = (cos * (x - cx)) + (sin * (y - cy)) + cx,
        ny = (cos * (y - cy)) - (sin * (x - cx)) + cy;
    return [nx, ny];
}

MetabolicGraph.prototype.filterReactions = function(list,data) {
	let reactions = [];

	for (var i=0; i<list.length; i++) {
		let r = null;

		// Convert to reaction object if necessary
		if (typeof list[i] == "string") {
			r = this.model.getReactionById(list[i]);
		} else {
			r = list[i];
		}

		// Is this reaction in the data?
		if (!r) continue;
		if (data && this.options.skipMissing && !data.hasOwnProperty(r.id)) continue;
		reactions.push(r);
	}

	return reactions;
}

MetabolicGraph.prototype.choosePrimary = function(set) {
	let rank = {};
	let maxr = 10000;
	let maxm = null;

	for (var x in set) {
		let m = this.model.index_metabolites[x];
		let r = 0;
		if (specialMetabolites[m.id]) continue;
		/*else if (!m.biocyc) r = 1000;
		else if (m.biocyc.inchi) r = m.biocyc.inchi.length;
		else {
			r = 900;
			//console.log("MISSING INCHI",m);
		}*/

		r = (m.producers.length+1) * (m.consumers.length+1);

		if (r < maxr) {
			maxr = r;
			maxm = m;
		}
	}

	if (maxm == null) console.error("MISSING PRIMARY");
	return maxm;
}

function wrapName(n) {
	if (n.length <= 20) return n;
	return wrap(n, {width: 20, cut: true});
}

MetabolicGraph.prototype.createMetabolite = function(m, primary) {
	let n = {
		type: "metabolite",
		id: m.id,
		primary: primary,
		colour: "white",
		value: 1,
		data: 1,
		name: wrapName(m.name),
		blocked: false,
		extra: !primary
	};

	if (m.name == "Pyruvate") {
		console.log("FOUND PYRUVATE", n);
		n.fixed = true;
		n.x = 1500;
		n.y = 1500;
		n.fx = 1500;
		n.fy = 1500;
		n.colour = "blue";
	} else if (this.cached[m.id]) {
		let c = this.cached[m.id];
		n.x = c[0];
		n.y = c[1];
		n.fx = c[0];
		n.fy = c[1];
		n.fixed = true;
		n.colour = "#eee";
	}
	return n;
}

MetabolicGraph.prototype.setReactions = function(list) {
	let metabData = this.options.metaboliteData;
	let reactData = this.options.reactionData;
	let reactions = this.filterReactions(list, reactData);

	//var mdatascale = 5.0 / maxmdata;

	// Generate nodes from metabolites
	/*var nodes = {};
	for (var x in metabs) {
		nodes[x] = {
			origin: metabs[x],
			"id": x,
			//"name": x,
			"val": 3,
			count: 0,
			type: "metabolite",
			name: metabs[x].name,
			data: (metabData && metabData.hasOwnProperty(x)) ? metabData[x] * mdatascale : 0,
			dead: (metabData && (!metabData.hasOwnProperty(x) || metabData[x] == 0)) || !metabData
		};
	}*/

	// Find max data value in reactions list
	var maxdata = 0.001;
	if (reactData) {
		for (var i=0; i<reactions.length; i++) {
			let data = reactData[reactions[i].id];
			if (Math.abs(data) > maxdata) maxdata = Math.abs(data);
		}
	}

	var datascale = 15.0 / maxdata;
	let maxcount = 0;

	let nodes = {};
	let links = [];

	for (var i=0; i<reactions.length; i++) {
		let data = (reactData && reactData.hasOwnProperty(reactions[i].id)) ? reactData[reactions[i].id] : 0
		let scaled = data*datascale;
		let blocked = Math.abs(data) < 0.00000000000001;
		let missing = reactData && !reactData.hasOwnProperty(reactions[i].id);

		if (this.options.hideZero && blocked) continue;

		if (missing && excluded[reactions[i].subsystem]) continue;

		// Choose primary in and out metabolites
		let prim_in = this.choosePrimary(reactions[i].inputs);
		let prim_out = this.choosePrimary(reactions[i].outputs);

		if (prim_in == null || prim_out == null) continue;

		if (!nodes[prim_in.id]) nodes[prim_in.id] = this.createMetabolite(prim_in, true);
		if (!nodes[prim_out.id]) nodes[prim_out.id] = this.createMetabolite(prim_out, true);
		//links[reactions[i].id] = createReaction(reactions[i], prim_in, prim_out);

		//if (prim_in === prim_out) continue;

		if (this.options.annotations && this.options.annotations[reactions[i].id]) {
			nodes[prim_out.id].annotated = true;
		}

		links.push({
			id: reactions[i].id,
			origin: reactions[i],
			source: nodes[prim_in.id],
			target: nodes[prim_out.id],
			input: false,
			val: scaled + ((scaled < 0) ? -1 : 1),
			blocked: blocked,
			missing: missing
		});

		// TODO Deal with non-primary nodes

		// Create links to each metabolite input
		/*for (var x in reactions[i].inputs) {
			if (specialMetabolites.hasOwnProperty(x)) {
				if (this.options.hideSpecials) continue;
				nodes[x+reactions[i].id] = {"id": x+reactions[i].id, "name": x, "val": 3, type: "metabolite", blocked: blocked, data: 0.0};
				links.push({origin: reactions[i], source: nodes[x+reactions[i].id], target: nodes[reactions[i].id], input: true, val: value, special: true, blocked: blocked, missing: missing});
			} else {
				nodes[x].count++;
				if (nodes[x].count > maxcount) maxcount = nodes[x].count;
				links.push({
					origin: reactions[i],
					source: nodes[x],
					target: nodes[reactions[i].id],
					input: true,
					val: value,
					blocked: blocked,
					missing: missing
				});
			}
		}*/

		// Create links to each metabolite output
		/*for (var x in reactions[i].outputs) {
			if (specialMetabolites.hasOwnProperty(x)) {
				if (this.options.hideSpecials) continue;
				nodes[x+reactions[i].id] = {"id": x+reactions[i].id, "name": x, "val": 3, type: "metabolite", blocked: blocked, data: 0.0};
				links.push({origin: reactions[i], source: nodes[reactions[i].id], target: nodes[x+reactions[i].id], val: value, special: true, blocked: blocked, missing: missing});
			} else {
				nodes[x].count++;
				if (nodes[x].count > maxcount) maxcount = nodes[x].count;
				links.push({origin: reactions[i], source: nodes[reactions[i].id], target: nodes[x], val: value, blocked: blocked, missing: missing});
			}
		}*/
	}

	// Filter isolated nodes and links (with only an in or out)
	/*if (this.options.removeIsolated) {
		for (var i=0; i<links.length; i++) {
			let n = (links[i].target.type == "metabolite") ? links[i].target : links[i].source;
			if (n.count <= 1) {
				links.splice(i,1);
				i--;
				delete nodes[n.id];
			}
		}
	}*/

	// TODO? Group metabolites by subsystem if possible

	// Add invisible links between reactions in same subsystem
	/*if (this.options.subsystemCluster) {
		let excluded = {
			"General": true,
			"Biomass": true,
			"Others": true,
			"Exchange reactions": true,
			"Modeling": true,
			"Transport": true
		};

		for (var x in this.model.subsystems) {
			if (excluded[x]) continue;

			let s = this.model.subsystems[x].reactions;
			for (var i=0; i<s.length; i++) {
				if (!nodes[s[i].id]) continue;
				for (var j=i+1; j<s.length; j++) {
					if (!nodes[s[j].id]) continue;
					links.push({source: nodes[s[i].id], target: nodes[s[j].id], val: 1, subsys: true, input: true});
				}
				//console.log("ADDING SUBSYS LINKS");
			}
		}
	}*/

	console.log("DATA",nodes,links);

	// Update the d3 graph
	this.graphData({nodes: nodes, links: links}, 40, -100);
}

function selectColor(colorNum, colors){
    if (colors < 1) colors = 1; // defaults to one color - avoid divide by zero
    return "hsl(" + Math.floor(colorNum * (360 / colors) % 360) + ",100%,50%)";
}

function pathStyle(d) {
	if (d.invisible) return "stroke: none";
	//if (d.subsys && !d.sig) return "stroke: rgba(0,0,255,0.4); stroke-width: 1px;";
	if (d.missing) return "stroke: rgba(100,100,100,0.5); stroke-width: 1px; stroke-dasharray: 2, 2;";
	if (d.blocked) return "stroke: rgb(0,0,255); stroke-width: 1px";
	let val = (d.special) ? 0.08 : (Math.abs(d.val)/15 + 0.33);
	if (d.val < 0) return "stroke-width: " + Math.ceil(Math.abs(d.val)) + "px; stroke: rgb(0,255,0)";
	return "stroke-width: " + Math.ceil(Math.abs(d.val)) + "px; stroke: rgb(255,0,0)";
}

function nodeStyle(d) {
	return (d.type == "reaction") ? "fill: "+d.colour+";" + ((d.annotated) ? " stroke-width: 5px;" : "") : "";
}

MetabolicGraph.prototype.displaySubsystems = function(force, svg, colours) {
	let nodes = force.nodes();
	let subsyspos = {};

	for (var i=0; i<nodes.length; i++) {
		let s = nodes[i].origin.subsystem;
		if (excluded[s]) continue;
		//if (nodes[i].extra) continue;

		if (!subsyspos[s]) {
			subsyspos[s] = [nodes[i].x,nodes[i].y,1,0];
		} else {
			subsyspos[s][0] += nodes[i].x;
			subsyspos[s][1] += nodes[i].y;
			subsyspos[s][2]++;
		}
	}

	for (var x in subsyspos) {
		subsyspos[x][0] /= subsyspos[x][2];
		subsyspos[x][1] /= subsyspos[x][2];
	}

	// Now find max distance from centre
	for (var i=0; i<nodes.length; i++) {
		let s = nodes[i].origin.subsystem;
		if (excluded[s]) continue;
		//if (nodes[i].extra) continue;

		var dx = nodes[i].x - subsyspos[s][0];
		var dy = nodes[i].y - subsyspos[s][1];
		var d = Math.sqrt(dx*dx + dy*dy);

		if (d+nodes[i].radius > subsyspos[s][3]) subsyspos[s][3] = d+nodes[i].radius;
	}

	let paddingLR = 10;
	let paddingTB = 10;

	for (var x in subsyspos) {
		let g = svg.insert("g",".node");
		g.append("circle")
		.attr("style", "fill: rgba(200,200,200,0.5)")
		.attr("stroke", "none")
		.attr("r", Math.floor(subsyspos[x][3])+20)
		.attr("title", x)
		.attr("cx", Math.floor(subsyspos[x][0]))
		.attr("cy", Math.floor(subsyspos[x][1]))
		.append("title").text(x);

		let r = g.append("rect");

		let t = g.append("text")
		.attr("y", Math.floor(subsyspos[x][1]) + Math.floor(subsyspos[x][3])+35)
		.attr("x", Math.floor(subsyspos[x][0]))
		.attr("text-anchor","middle")
		.text(x);
		let bb = null; 
		t.each(function() { bb = this.getBBox(); }); 

		r.attr("x", subsyspos[x][0] - bb.width/2 - paddingLR)
		.attr("y", subsyspos[x][1] + Math.floor(subsyspos[x][3])+35 - bb.height/2 - paddingTB)
		.attr("rx", 5)
		.attr("ry", 5)
		.attr("stroke", "#222")
		.attr("stroke-width", 1)
		.attr("width", bb.width + paddingLR*2)
		.attr("height", bb.height + paddingTB)
		.attr("fill", "white");
	}
}

MetabolicGraph.prototype.graphData = function(data, dist, charge, cols) {
	var width = 3000,
    height = 3000;
	let paddingLR = 10;
	let paddingTB = 5;


	var force = d3.layout.force()
		.nodes(d3.values(data.nodes))
		.links(data.links)
		.size([width, height])
		.gravity((this.options.noGravity) ? 0 : 0.01)
		.linkDistance(dist) //60
		.charge(charge) // -300
		.linkStrength(link => (link.missing) ? 0.3 : (link.blocked) ? 0.5 : 1)
		.on("tick", tick)
		.start();
	this.force = force;

	let svg = this.svg;

	svg.selectAll("*").remove();

	svg.append("rect")
		.attr("fill", "none")
		.attr("pointer-events", "all")
		.attr("width", 3000)
		.attr("height", 3000)
		.call(d3.behavior.zoom()
		    .scaleExtent([0.2, 8])
		    .on("zoom", zoom));

	function zoom() {
	  g.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
	}

	// build the arrow.
	svg.append("svg:defs").selectAll("marker")
		.data(["small","medium","large"])      // Different link/path types can be defined here
	  .enter().append("svg:marker")    // This section adds in the arrows
		.attr("id", String)
		.attr("viewBox", "0 -5 10 10")
		.attr("refX", d => (d == "small") ? 4 : (d == "medium") ? 6 : 10) //15
		.attr("refY", 0) //-1.5
		.attr("markerWidth", d => (d == "small") ? 8 : (d == "medium") ? 12 : 20)
		.attr("markerHeight", d => (d == "small") ? 8 : (d == "medium") ? 12 : 20)
		.attr("markerUnits","userSpaceOnUse")
		//.attr("style","stroke-width: 3px")
		//.attr("style", "stroke-width: 1.5px")
		.attr("orient", "auto")
	  .append("svg:path")
		.attr("d", "M0,-5L10,0L0,5");

	let g = svg.append("svg:g");

	let me = this;
	force.on("end", function() {
		console.log("END", force);
		//me.displaySubsystems(force,g,cols);
	});

	// add the links and the arrows
	var path = g.append("svg:g").selectAll("path")
		.data(force.links())
	  .enter().append("svg:path")
	//    .attr("class", function(d) { return "link " + d.type; })
		.attr("class", function(d) { return (((d.special) ? "link special" : "link") + ((d.blocked) ? " blocked" : "") + ((d.val < 0) ? " negative" : "") + ((d.missing) ? " missing" : "")); })
		.attr("style", pathStyle)
		.attr("fill","none")
		.attr("marker-mid", function(d) { return (d.input) ? "" : (d.val <= 5) ? "url(#small)" : (d.val <= 9) ? "url(#medium)" : "url(#large)"; });

	var drag = force.drag().on("dragstart", function(d) {
		d.fixed = true;
		d.fx = d.x;
		d.fy = d.y;
		d.colour = "#eee";
		this.setAttribute("fill", "#eee");
	}).on("dragend", function(d) {
		me.cached[d.id] = [d.x,d.y];
		window.localStorage.gemvis_cached = JSON.stringify(me.cached);
	});

	// define the nodes
	var node = g.selectAll(".node")
		.data(force.nodes())
	  .enter().append("g")
		.attr("class", function(d) { return ((d.special) ? "node special" : "node") + " " + d.type + ((d.blocked) ? " blocked" : "") + ((d.dead) ? " dead" : "") + ((d.val < 0) ? " negative" : ""); })
		.call(drag);

	// add the nodes
	/*node.append("circle")
		.attr("style", d => nodeStyle(d))
		.attr("fill", d => (d.colour) ? d.colour : "#eee")
		.attr("stroke", d=> (d.extra || d.blocked) ? "none" : (d.annotated) ? "yellow" : "#222")
		.attr("r", function(d) { return (d.type == "metabolite") ? 3 + Math.floor(Math.abs(d.data)*4) : 4 + Math.floor(Math.abs(d.data)*8); })
		.attr("title", function(d) { return d.name; })
		.append("title").text(function(d) { return (d.type == "metabolite" && d.fullname) ? d.fullname : d.name; });*/

	let r = node.append("rect");

	//if (!this.options.hideMetaboliteNames) {
		// add the text 
		let t = node.append("text")
		//.attr("y", Math.floor(subsyspos[x][1]) + Math.floor(subsyspos[x][3])+35)
		//.attr("x", Math.floor(subsyspos[x][0]))
		.attr("text-anchor","middle");
		//.text(d => d.name);

		t.selectAll("tspan")
		.data(function(d) { return d.name.split("\n"); })
		.enter()
		.append("tspan")
		.attr("x",0)
		.attr("dy","1.2em")
		.text(d => d);

		//let bb = null; 
		//t.each(function() { bb = this.getBBox(); }); 

		//r.attr("x", - bb.width/2 - paddingLR)
		//.attr("y", - bb.height/2 - paddingTB)
		r.attr("rx", 5)
		.attr("ry", 5)
		.attr("stroke", d => (d.annotated) ? "orange" : "#222")
		.attr("stroke-width", d => (d.annotated) ? 3 : 2)
		//.attr("width", bb.width + paddingLR*2)
		//.attr("height", bb.height + paddingTB)
		.attr("fill", d => d.colour);
	//}

	if (this.options.showReactionNames) {
		// add the text 
		node.append("text")
			.attr("x", 12)
			.attr("dy", ".35em")
			.text(function(d) { return (d.type == "metabolite") ? "" : ((d.name.length > 10) ? d.name.substring(0,10)+"..." : d.name); });
	}

	// add the curvy lines
	function tick() {
		path.attr("d", function(d) {
		    var dx = d.target.x - d.source.x,
		        dy = d.target.y - d.source.y,
		        dr = Math.sqrt(dx * dx + dy * dy);
		    return "M" + 
		        d.source.x + "," + 
		        d.source.y + "A" + 
		        dr + "," + dr + " 0 0," + ((d.input) ? "0" : "1") + " " + 
		        d.target.x + "," + 
		        d.target.y;
		});

		t.each(function() {
			let bb = this.getBBox();
			let r = this.previousSibling;

			this.setAttribute("y", -bb.height/2 - paddingTB);

			r.setAttribute("x", - bb.width/2 - paddingLR);
			r.setAttribute("y", -bb.height/2 - paddingTB);
			r.setAttribute("width", bb.width + paddingLR*2);
			r.setAttribute("height", bb.height + 2*paddingTB);
		});

		node
		    .attr("transform", function(d) { 
	  	    return "translate(" + d.x + "," + d.y + ")"; });
	}
}

module.exports = MetabolicGraph;


