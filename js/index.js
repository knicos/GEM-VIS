//const ForceGraph3D = require('3d-force-graph');
//const metanal = require('../metabolite-analysis.js');
//const saveSvgAsPng = require("save-svg-as-png");

function MetabolicGraph(parent, model, options) {
	this.model = model;
	this.options = (options) ? options : {};
	this.force = null;

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

	"M_atp_p": true,
	"M_nad_p": true,
	"M_h_p": true,
	"M_pi_p": true,
	"M_adp_p": true,
	"M_h2o_p": true,
	"M_o2_p": true,
	"M_nadph_p": true,
	"M_co2_p": true,
	"M_nadh_p": true
};

MetabolicGraph.prototype.setReactions = function(list) {
	if (this.options.hideMetabolites) {
		this.setReactions_hidemetabs(list);
	} else {
		this.setReactions_showmetabs(list);
	}
}

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
	let value = (data1 / m.consumers.length) * (data2*0.5 + 1); // / m.producers.length

	return [value, sig];
}

MetabolicGraph.prototype.setReactions_hidemetabs = function(list) {
	let reactions = [];
	let reactData = this.options.reactionData;

	for (var i=0; i<list.length; i++) {
		let r = null;

		if (typeof list[i] == "string") {
			r = this.model.getReactionById(list[i]);
		} else {
			r = list[i];
		}

		if (!r) continue;
		if (reactData && this.options.skipMissing && !reactData.hasOwnProperty(r.id)) continue;
		reactions.push(r);

		/*for (var x in r.metabolites) {
			if (!specialMetabolites.hasOwnProperty(x)) {
				if (r.metabolites[x]) {
					metabs[x] = r.metabolites[x];
				} else {
					metabs[x] = {
						name: "MISSING"
					};
				}
			}
		}*/
	}

	var nodes = {};
	var links = [];
	var maxdata = 0.001;
	if (reactData) {
		for (var i=0; i<reactions.length; i++) {
			let data = reactData[reactions[i].id];
			if (Math.abs(data) > maxdata) maxdata = Math.abs(data);
		}
	}

	// Tally represented subsystems.
	let subsystems = {};
	for (var i=0; i<reactions.length; i++) {
		subsystems[reactions[i].subsystem] = true;
	}

	// Generate colours and positions for subsystems
	let subsys_colours = {};
	let subsys_positions = {};
	let sscount = 0;
	let sscount2 = 0;
	for (var x in subsystems) {
		sscount++;
	}
	for (var x in subsystems) {
		subsys_colours[x] = selectColor(sscount2,sscount);
		subsys_positions[x] = rotate(1500,1500, 1900, 1500, ((2*Math.PI) / sscount) * sscount2);
		sscount2++;
	}

	var datascale = 10.0 / maxdata;
	let maxcount = 0;

	for (var i=0; i<reactions.length; i++) {
		let data = (reactData && reactData.hasOwnProperty(reactions[i].id)) ? reactData[reactions[i].id] : 0
		let value = data*datascale + 1;
		let blocked = Math.abs(data) < 0.00000000000001;
		//let col = (excluded[reactions[i].subsystem]) ? "#eee" : subsys_colours[reactions[i].subsystem];
		let col = (blocked) ? "blue" : (data < 0) ? "green" : "red";
		let pos = subsys_positions[reactions[i].subsystem];

		if (this.options.hideZero && blocked) continue;

		// Create a reaction node
		nodes[reactions[i].id] = {
			origin: reactions[i],
			id: reactions[i].id,
			colour: col,
			name: reactions[i].name,
			data: data,
			val: value,
			radius: 4 + Math.floor(Math.abs(data)*8),
			count: 0,
			type: "reaction",
			blocked: blocked,
			annotated: this.options.annotations && this.options.annotations[reactions[i].id],
			x: pos[0], y: pos[1]
		};
	}

	for (var i=0; i<reactions.length; i++) {
		let data = (reactData && reactData.hasOwnProperty(reactions[i].id)) ? reactData[reactions[i].id] : 0;
		//let value = data*datascale + 1;
		let blocked = Math.abs(data) < 0.00000000000001;
		let missing = reactData && !reactData.hasOwnProperty(reactions[i].id);

		if (this.options.hideZero && blocked) continue;

		for (var x in reactions[i].outputs) {
			if (!reactions[i].metabolites[x]) continue;

			/*let value = data; //(data / Math.sqrt(reactions[i].metabolites[x].producers.length));
			//if (Math.abs(value) < 0.01*maxdata) continue;
			let insig = (Math.pow(Math.abs(value)+1,2)-1) / reactions[i].metabolites[x].producers.length < 0.01*maxdata*maxdata;
			value /= reactions[i].metabolites[x].consumers.length;
			if (insig) continue;*/

			let visited = {};

			//if (specialMetabolites.hasOwnProperty(x)) {
				//nodes[x+reactions[i].id] = {"id": x+reactions[i].id, "name": x, "val": 3, type: "metabolite", blocked: blocked, data: 0.0};
				//links.push({source: nodes[x+reactions[i].id], target: nodes[reactions[i].id], input: true, val: value, special: true, blocked: blocked, missing: missing});
			//} else {
				for (var j=0; j<reactions[i].metabolites[x].consumers.length; j++) {
					let con = reactions[i].metabolites[x].consumers[j];

					if (visited[con.id]) continue;
					visited[con.id] = true;

					let cdata = (reactData && reactData.hasOwnProperty(con.id)) ? reactData[con.id] : 0;
					let [value,edgeSig] = calculateEdge(data, cdata, reactions[i].metabolites[x], maxdata)

					if (!edgeSig) continue;

					// Insert missing target reactions
					if (!nodes[con.id]) {

						// Exclude some (virtual) reactions always
						if (con.subsystem == "Biomass" || con.subsystem == "General") continue;

						// Create a reaction node
						nodes[con.id] = {
							extra: true,
							origin: con,
							id: con.id,
							colour: "#222", //subsys_colours[con.subsystem],
							name: con.name,
							radius: 4 + Math.floor(0.1*8),
							data: 0.1,
							val: 0,
							count: 0,
							type: "reaction"
						};
					}

					nodes[con.id].count++;
					nodes[reactions[i].id].count++;
					//links.push({source: nodes[reactions[i].metabolites[x].reactions[j].id], target: nodes[reactions[i].id], input: true, val: 1, blocked: blocked, missing: true});
					if (nodes[con.id].count > maxcount) maxcount = nodes[con.id].count;
									
					links.push({
						origin: reactions[i],
						target: nodes[con.id],
						source: nodes[reactions[i].id],
						input: false,
						val: value*3*datascale, // + ((value < 0) ? -1 : 1),
						blocked: Math.abs(value) < 0.001,
						missing: missing,
						special: nodes[con.id].extra == true
					});
				}
			//}
		}
	}

	// Filter isolated nodes and links (with only an in or out)
	if (this.options.removeIsolated) {
		for (var x in nodes) {
			let n = nodes[x];
			if (n.count <= 0) {
				//nodes.splice(i,1);
				//i--;
				delete nodes[x];
			}
		}
	}

	// TODO? Group metabolites by subsystem if possible

	let subsyscount = {};

	// Add invisible links between reactions in same subsystem
	// Also add subsystem nodes and link reactions to subsystem nodes
	if (this.options.subsystemCluster) {
		for (var x in this.model.subsystems) {
			if (excluded[x]) continue;

			subsyscount[x] = 0;
			let s = this.model.subsystems[x].reactions;
			for (var i=0; i<s.length; i++) {
				if (!nodes[s[i].id]) continue;
				for (var j=i+1; j<s.length; j++) {
					if (!nodes[s[j].id]) continue;
					subsyscount[x]++;
					links.push({
						source: nodes[s[i].id],
						name: x,
						target: nodes[s[j].id],
						val: 1,
						subsys: true,
						input: true});
				}
				//console.log("ADDING SUBSYS LINKS");
			}
		}
	}

	// Update the d3 graph
	let me = this;
	this.graphData({nodes: nodes, links: links}, function(link) {
		let count = ((link.target.count > link.source.count) ? link.target.count-1 : link.source.count-1) / (maxcount-1);
		return (link.subsys) ? Math.floor(10 + 5*Math.sqrt(subsyscount[link.name])) : Math.floor(40 + count*20);
	}, -300, subsys_colours);
}

function rotate(cx, cy, x, y, angle) {
    var radians = (Math.PI / 180) * angle,
        cos = Math.cos(radians),
        sin = Math.sin(radians),
        nx = (cos * (x - cx)) + (sin * (y - cy)) + cx,
        ny = (cos * (y - cy)) - (sin * (x - cx)) + cy;
    return [nx, ny];
}

//----------------------------------------------------------------
// Include METABOLITES
//----------------------------------------------------------------

MetabolicGraph.prototype.setReactions_showmetabs = function(list) {
	let reactions = [];
	let metabs = {};
	let maxmdata = 0.0;
	let metabData = this.options.metaboliteData;
	let reactData = this.options.reactionData;

	for (var i=0; i<list.length; i++) {
		let r = null;

		if (typeof list[i] == "string") {
			r = this.model.getReactionById(list[i]);
		} else {
			r = list[i];
		}

		if (!r) continue;
		if (reactData && this.options.skipMissing && !reactData.hasOwnProperty(r.id)) continue;
		reactions.push(r);

		for (var x in r.metabolites) {
			if (!specialMetabolites.hasOwnProperty(x)) {
				if (r.metabolites[x]) {
					metabs[x] = r.metabolites[x];
					//metanal.productionFlux(metabs[x]);
					if (metabData) {
						let data = metabData[x];
						if (data > maxmdata) maxmdata = data;
					}	
				} else {
					metabs[x] = {
						name: "MISSING"
					};
				}
			}
		}
	}
	var mdatascale = 5.0 / maxmdata;

	// Generate nodes from metabolites
	var nodes = {};
	for (var x in metabs) {
		nodes[x] = {
			origin: metabs[x],
			"id": x,
			"name": x,
			"val": 3,
			count: 0,
			type: "metabolite",
			fullname: metabs[x].name,
			data: (metabData && metabData.hasOwnProperty(x)) ? metabData[x] * mdatascale : 0,
			dead: (metabData && (!metabData.hasOwnProperty(x) || metabData[x] == 0)) || !metabData
		};
	}

	// Find max data value in reactions list
	var links = [];
	var maxdata = 0.001;
	if (reactData) {
		for (var i=0; i<reactions.length; i++) {
			let data = reactData[reactions[i].id];
			if (Math.abs(data) > maxdata) maxdata = Math.abs(data);
		}
	}

	// Generate colours and positions for subsystems
	let subsys_colours = {};
	let subsys_positions = {};
	let sscount = 0;
	let sscount2 = 0;
	for (var x in this.model.subsystems) {
		sscount++;
	}
	for (var x in this.model.subsystems) {
		subsys_colours[x] = selectColor(sscount2,sscount);
		subsys_positions[x] = rotate(1500,1500, 1900, 1500, ((2*Math.PI) / sscount) * sscount2);
		sscount2++;
	}

	var datascale = 10.0 / maxdata;
	let maxcount = 0;

	for (var i=0; i<reactions.length; i++) {
		let data = (reactData && reactData.hasOwnProperty(reactions[i].id)) ? reactData[reactions[i].id] : 0
		let value = data*datascale + 1;
		let blocked = Math.abs(data) < 0.00000000000001;
		let missing = reactData && !reactData.hasOwnProperty(reactions[i].id);
		let col = subsys_colours[reactions[i].subsystem];
		let pos = subsys_positions[reactions[i].subsystem];

		if (this.options.hideZero && blocked) continue;

		// Create a reaction node
		nodes[reactions[i].id] = {
			origin: reactions[i],
			id: reactions[i].id,
			colour: col,
			name: reactions[i].name,
			data: data,
			val: value,
			type: "reaction"
			//x: pos[0], y: pos[1]
		};

		// Create links to each metabolite input
		for (var x in reactions[i].inputs) {
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
		}

		// Create links to each metabolite output
		for (var x in reactions[i].outputs) {
			if (specialMetabolites.hasOwnProperty(x)) {
				if (this.options.hideSpecials) continue;
				nodes[x+reactions[i].id] = {"id": x+reactions[i].id, "name": x, "val": 3, type: "metabolite", blocked: blocked, data: 0.0};
				links.push({origin: reactions[i], source: nodes[reactions[i].id], target: nodes[x+reactions[i].id], val: value, special: true, blocked: blocked, missing: missing});
			} else {
				nodes[x].count++;
				if (nodes[x].count > maxcount) maxcount = nodes[x].count;
				links.push({origin: reactions[i], source: nodes[reactions[i].id], target: nodes[x], val: value, blocked: blocked, missing: missing});
			}
		}
	}

	// Filter isolated nodes and links (with only an in or out)
	if (this.options.removeIsolated) {
		for (var i=0; i<links.length; i++) {
			let n = (links[i].target.type == "metabolite") ? links[i].target : links[i].source;
			if (n.count <= 1) {
				links.splice(i,1);
				i--;
				delete nodes[n.id];
			}
		}
	}

	// TODO? Group metabolites by subsystem if possible

	// Add invisible links between reactions in same subsystem
	if (this.options.subsystemCluster) {
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
	}

	// Update the d3 graph
	this.graphData({nodes: nodes, links: links}, function(link) {
		/*if (link.target.type == "metabolite") {
			if (link.target.origin.subsystems && link.target.origin.subsystems.hasOwnProperty(link.origin.subsystem)) {
				return 100 - 70 * link.target.origin.subsystems[link.origin.subsystem];
			} else {
				return 100;
			}
		} else {
			if (link.source.origin.subsystems && link.source.origin.subsystems.hasOwnProperty(link.origin.subsystem)) {
				return 100 - 70 * link.source.origin.subsystems[link.origin.subsystem];
			} else {
				return 100;
			}
		}*/
		//return 60; //130 - 10 * Math.abs(link.val);
		//return 100; //(link.subsys) ? 100 : 40;

		let count = ((link.target.type == "metabolite") ? link.target.count-1 : link.source.count-1) / (maxcount-1);
		return (link.subsys) ? 100 : Math.floor(40 + count*100);
	}, -200);
}

function selectColor(colorNum, colors){
    if (colors < 1) colors = 1; // defaults to one color - avoid divide by zero
    return "hsl(" + Math.floor(colorNum * (360 / colors) % 360) + ",100%,50%)";
}

function pathStyle(d) {
	if (d.subsys) return "";
	if (d.missing) return "";
	if (d.blocked) return "stroke: rgba(0,0,255,0.1); stroke-width: 1px";
	let val = (d.special) ? 0.08 : (Math.abs(d.val)/15 + 0.33);
	if (d.val < 0) return "stroke-width: " + Math.ceil(Math.abs(d.val)) + "px; stroke: rgba(0,255,0,"+val+")";
	return "stroke-width: " + Math.ceil(Math.abs(d.val)) + "px; stroke: rgba(255,0,0,"+val+")";
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

	var force = d3.layout.force()
		.nodes(d3.values(data.nodes))
		.links(data.links)
		.size([width, height])
		//.gravity(0.05)
		.linkDistance(dist) //60
		.charge(charge) // -300
		.linkStrength(link => (link.subsys) ? 2 : (link.blocked) ? 0 : 0.05)
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
		.data(["end"])      // Different link/path types can be defined here
	  .enter().append("svg:marker")    // This section adds in the arrows
		.attr("id", String)
		.attr("viewBox", "0 -5 10 10")
		.attr("refX", 15) //15
		.attr("refY", -1.5)
		.attr("markerWidth", 8)
		.attr("markerHeight", 8)
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
		me.displaySubsystems(force,g,cols);
	});

	// add the links and the arrows
	var path = g.append("svg:g").selectAll("path")
		.data(force.links())
	  .enter().append("svg:path")
	//    .attr("class", function(d) { return "link " + d.type; })
		.attr("class", function(d) { return (d.subsys) ? "link invisible" : (((d.special) ? "link special" : "link") + ((d.blocked) ? " blocked" : "") + ((d.val < 0) ? " negative" : "") + ((d.missing) ? " missing" : "")); })
		.attr("style", pathStyle)
		.attr("fill","none");
		//.attr("marker-mid", function(d) { return (d.input) ? "" : "url(#end)"; });

	// define the nodes
	var node = g.selectAll(".node")
		.data(force.nodes())
	  .enter().append("g")
		.attr("class", function(d) { return ((d.special) ? "node special" : "node") + " " + d.type + ((d.blocked) ? " blocked" : "") + ((d.dead) ? " dead" : "") + ((d.val < 0) ? " negative" : ""); })
		.call(force.drag);

	// add the nodes
	node.append("circle")
		.attr("style", d => nodeStyle(d))
		.attr("fill", d => (d.colour) ? d.colour : "#eee")
		.attr("stroke", d=> (d.extra || d.blocked) ? "none" : (d.annotated) ? "yellow" : "#222")
		.attr("r", function(d) { return (d.type == "metabolite") ? 3 + Math.floor(Math.abs(d.data)*4) : 4 + Math.floor(Math.abs(d.data)*8); })
		.attr("title", function(d) { return d.name; })
		.append("title").text(function(d) { return (d.type == "metabolite" && d.fullname) ? d.fullname : d.name; });

	if (!this.options.hideMetaboliteNames) {
		// add the text 
		node.append("text")
			.attr("x", 12)
			.attr("dy", ".35em")
			.text(function(d) { return (d.type == "metabolite" && !d.noname) ? d.name.substring(2) : ""; });
	}

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

		node
		    .attr("transform", function(d) { 
	  	    return "translate(" + d.x + "," + d.y + ")"; });
	}
}

module.exports = MetabolicGraph;


