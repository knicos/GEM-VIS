//const ForceGraph3D = require('3d-force-graph');
//const metanal = require('../metabolite-analysis.js');
//const saveSvgAsPng = require("save-svg-as-png");
const biocyc = require('./compounds.js');
const wrap = require('word-wrap');

function processMetabolite(m) {
	// First add id2
	m.id2 = m.id.substring(0,m.id.length-2);

	for (var x in biocyc) {
		let b = biocyc[x];
		let pname = m.name.toLowerCase().replace(/\-/g," ");

		if (b && b.name && pname == b.name.toLowerCase().replace(/\-/g," ").replace("&beta;","beta").replace("&alpha;","alpha")) {
			return b;
		}

		if (b.synonyms) {
			for (var i=0; i<b.synonyms.length; i++) {
				if (b.synonyms[i].toLowerCase().replace(/\-/g," ").replace("&beta;","beta").replace("&alpha;","alpha") == pname) {
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

	this.tert_nodes = [];
	this.metaProducers = {};
	this.metaConsumers = {};
	this.metaRank = {};
	this.metaRatio = {};

	// TODO Process metabolites with BioCyc Data
	for (var i=0; i<model.metabolites.length; i++) {
		let m = model.metabolites[i];

		m.biocyc = processMetabolite(m);
		if (!this.metaProducers[m.id2]) this.metaProducers[m.id2] = 0;
		if (!this.metaConsumers[m.id2]) this.metaConsumers[m.id2] = 0;
		this.metaProducers[m.id2] += m.producers.length;
		this.metaConsumers[m.id2] += m.consumers.length;
		
		//if (!this.metaRank[m.id2]) this.metaRank[m.id2] = 0;
		//this.metaRank[m.id2] += (m.producers.length+1) * (m.consumers.length+1);
	}

	// Calculate metabolite metrics
	//let maxratio = 0;
	for (var x in this.metaProducers) {
		this.metaRank[x] = (this.metaProducers[x]+1) * (this.metaConsumers[x]+1);
		this.metaRatio[x] = 1.0 / this.metaRank[x];
		//if (this.metaRatio[x] > maxratio) maxratio = this.metaRatio[x];
	}
	// Normalise (and log) ratio
	//for (var x in this.metaRatio) this.metaRatio[x] = this.metaRatio[x] / maxratio;

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
	let list = [];

	for (var x in set) {
		let m = this.model.index_metabolites[x];
		let r = 0;
		//if (specialMetabolites[m.id]) continue;

		// TODO Consider the set number, divide by it?
		// TODO Must consider compartment independent numbers.
		//r = ((m.producers.length+1) * (m.consumers.length+1)); // * set[x]; // / set[x];
		r = this.metaRank[m.id2];

		rank[x] = r;
		list.push(m);

		if (r < maxr) {
			maxr = r;
			maxm = m;
		}
	}

	list = list.sort((a,b) => rank[a.id] - rank[b.id]);
	if (this.options.linkThreshold) {
		let thresh = this.options.linkThreshold;
		let listA = list.filter((a) => rank[a.id] / thresh <= maxr);
		let listB = list.filter((a) => rank[a.id] / thresh > maxr);
		return [listA,listB];
	}

	if (maxm == null) console.error("MISSING PRIMARY");
	return [list,[]];
}

function wrapName(n) {
	if (n.length <= 20) return n;

	let s = n.split(/(?=[\s\-]+)/);
	let r = [""];
	for (var i=0; i<s.length; i++) {
		if (r[r.length-1].length + s[i].length > 20) r.push("");
		r[r.length-1] += s[i];
	}
	return r.join("\n");


	//n = n.replace("alpha-","&alpha;-").replace("beta-","&beta;-").replace("gamma-","&gamma;-");
	//if (n.length <= 20) return n;
	//return wrap(n, {width: 20, cut: true});
}

MetabolicGraph.prototype.scanMetabolites = function(m, visited) {
	let upcount = 0;
	let downcount = 0;

	if (visited[m.id]) return 0;
	visited[m.id] = true;

	// Calculate metabolite colour from data
	for (var i=0; i<m.producers.length; i++) {
		if (this.options.reactionData[m.producers[i].id] > 0) upcount++;
		else if (this.options.reactionData[m.producers[i].id] < 0) downcount++;
		else {
			for (var x in m.producers[i].inputs) {
				let a = this.scanMetabolites(m.producers[i].metabolites[x], visited);
				if (a > 0) upcount++;
				else if (a < 0) downcount++;
			}
		}
	}

	if (upcount > 0 && downcount == 0) return 1;
	else if (downcount > 0 && upcount == 0) return -1;
	else return 0;
}

MetabolicGraph.prototype.colourFromData = function(m) {
	let colour = "white";
	if (this.options.reactionData) {
		let upcount = 0;
		let downcount = 0;

		// Calculate metabolite colour from data
		for (var i=0; i<m.producers.length; i++) {
			if (this.options.reactionData[m.producers[i].id] > 0) upcount++;
			else if (this.options.reactionData[m.producers[i].id] < 0) downcount++;
		}

		if (upcount > 0 && downcount == 0) colour = "rgb(255,100,100)";
		else if (downcount > 0 && upcount == 0) colour = "rgb(100,255,100)";
		else {
			let a = this.scanMetabolites(m,{});
			if (a > 0) colour = "rgb(255,180,180)";
			if (a < 0) colour = "rgb(180,255,180)";
		}
	}
	return colour;
}

MetabolicGraph.prototype.hasData = function(m) {
	if (!this.options.reactionData) return false;
	for (var i=0; i<m.reactions.length; i++) {
		if (this.options.reactionData[m.reactions[i].id]) return true;
	}
	return false;
}

MetabolicGraph.prototype.createMetabolite = function(m, primary, pid) {
	let colour = "white";
	let id = (pid) ? pid : m.id2;

	let slist = [];
	for (var x in m.subsystems) {
		slist.push(x);
	}
	// TODO Check sort order.
	slist = slist.sort((a,b) => m.subsystems[b] - m.subsystems[a]);

	switch (this.options.colourRule) {
	case "Subsystem"			: if (slist.length > 0) colour = this.subsys_colours[slist[0]]; break;
	case "Data Boolean"			: colour = (this.hasData(m)) ? "#ddd" : "white"; break;
	case "Data Extrapolated"	: colour = this.colourFromData(m); break;
	case "Significance Boolean"	: colour = (this.metaRatio[m.id2] < 0.01) ? "#aaa" : "white"; break;
	default:
	}

	let name = m.name;
	if (m.biocyc && m.biocyc.abbreviation) {
		name = m.biocyc.abbreviation.toUpperCase();
	} else if (m.biocyc && m.biocyc.synonyms) {
		// Find shortest synonym
		let s = m.biocyc.synonyms.sort((a,b) => a.length - b.length);
		name = (s[0].length < name.length) ? s[0] : name;
	}

	let n = {
		type: "metabolite",
		id: id,
		origin: m,
		primary: primary,
		colour: colour,
		count: 0,
		value: 1,
		data: 1,
		fullname: m.name,
		name: wrapName(name),
		blocked: false,
		extra: !primary,
		subsystem: (slist.length > 0) ? slist[0] : "None"
	};

	/*if (m.name == "Pyruvate") {
		console.log("FOUND PYRUVATE", n);
		n.fixed = true;
		n.x = 1500;
		n.y = 1500;
		n.fx = 1500;
		n.fy = 1500;
		//n.colour = "blue";
	} else */

	if (this.cached[id]) {
		let c = this.cached[id];
		n.x = c[0];
		n.y = c[1];
		n.fx = c[0];
		n.fy = c[1];
		n.fixed = true;
		//n.colour = "#eee";
	}/* else if (this.cached[m.id]) {
		let c = this.cached[m.id];
		n.x = c[0];
		n.y = c[1];
		n.fx = c[0];
		n.fy = c[1];
		n.fixed = true;
		//n.colour = "#eee";
	}*/
	return n;
}

MetabolicGraph.prototype.createMergedMetabolite = function(m, primary, pid) {
	let colour = "white";
	let id = pid;

	let name = "";

	for (var i=0; i<m.length; i++) {
		let tname = m[i].name;
		if (m[i].biocyc && m[i].biocyc.abbreviation) {
			tname = m[i].biocyc.abbreviation.toUpperCase();
		} else if (m[i].biocyc && m[i].biocyc.synonyms) {
			// Find shortest synonym
			let s = m[i].biocyc.synonyms.sort((a,b) => a.length - b.length);
			tname = (s[0].length < tname.length) ? s[0] : tname;
		}
		name += wrapName(tname) + "\n";
	}

	let n = {
		type: "metabolite",
		id: id,
		origin: m[0],
		primary: primary,
		colour: colour,
		count: 0,
		value: 1,
		data: 1,
		fullname: name,
		name: name.trim(),
		blocked: false,
		extra: !primary,
		subsystem: "Mixed"
	};
	return n;
}

MetabolicGraph.prototype.makePrimary = function(reaction, m1, m2, nodes, links, blocked, missing, scaled, data) {
	/*let reactData = this.options.reactionData;
	let data = (reactData && reactData.hasOwnProperty(reaction.id)) ? reactData[reaction.id] : 0
	let scaled = data*datascale;
	let blocked = Math.abs(data) < 0.00000000000001;
	let missing = reactData && !reactData.hasOwnProperty(reactions[i].id);*/

	if (!nodes[m1.id2]) nodes[m1.id2] = this.createMetabolite(m1, true);
	if (!nodes[m2.id2]) nodes[m2.id2] = this.createMetabolite(m2, true);
	let n1 = nodes[m1.id2];
	let n2 = nodes[m2.id2];

	//links[reactions[i].id] = createReaction(reactions[i], prim_in, prim_out);

	//if (prim_in === prim_out) continue;

	if (this.options.hideLoose && !n1.fixed && !n2.fixed) return;

	if (this.options.annotations && this.options.annotations[reaction.id]) {
		n2.annotated = true;
	}

	// If the nodes are fixed and too far apart, create a ghost node
	if (n1.fixed && n2.fixed) {
		var dx = n1.x - n2.x,
	        dy = n1.y - n2.y,
	        dr = Math.sqrt(dx * dx + dy * dy);
		if (dr > 500) {
			// Create ghost node for source
			n1 = this.createMetabolite(n1.origin, true);
			n1.colour = "#bbb";
			n1.ghost = true;
			n1.id += "_ghost_"+reaction.id+"_"+n2.id;

			// Does this ghost node have a layout cache?
			if (this.cached[n1.id]) {
				n1.fixed = true;
				n1.fx = this.cached[n1.id][0];
				n1.fy = this.cached[n1.id][1];
				n1.x = n1.fx;
				n1.y = n1.fy;
			} else {
				n1.fixed = false;
				delete n1.fx;
				delete n2.fy;
			}
			if (nodes[n1.id]) console.error("EXISTING GHOST", nodes[n1.id]);
			nodes[n1.id] = n1;
		}
	}

	n1.count++;
	n2.count++;

	let l = {
		id: reaction.id,
		origin: reaction,
		source: n1,
		target: n2,
		input: false,
		val: scaled + ((scaled < 0) ? -1 : 1),
		blocked: blocked,
		missing: missing,
		primary: true,
		data: data
	};
	links.push(l);

	if (reaction.reversible) {
		links.push({
			id: reaction.id+"_rev",
			origin: reaction,
			target: n1,
			source: n2,
			input: false,
			val: scaled + ((scaled < 0) ? -1 : 1),
			blocked: blocked,
			missing: missing,
			primary: true,
			data: data
		});
	}

	return l;
}

let secindex = 0;

MetabolicGraph.prototype.makeSecondary = function(reaction, n_tert, m, nodes, links, master, input) {
	//let id = m.id2 + "_SEC_" + reaction.id + "_" + p1.id2;

	if (!nodes[m.id2]) nodes[m.id2] = this.createMetabolite(m, false);
	//if (!nodes[m2.id2]) nodes[m2.id2] = this.createMetabolite(m2, true);
	let n = nodes[m.id2];
	//n.id = n.id + "_"+reaction.id+"_"+p1.id2;
	//let n2 = nodes[m2.id2];

	//let np1 = nodes[p1.id2];
	//let np2 = nodes[p2.id2];

	//if (!np1 || !np2) {
	//	console.error("MISSING PRIMARY");
	//	return;
	//}

	//links[reactions[i].id] = createReaction(reactions[i], prim_in, prim_out);

	//if (prim_in === prim_out) continue;

	//if (this.options.hideLoose && !n1.fixed && !n2.fixed) return;

	//if (this.options.annotations && this.options.annotations[reaction.id]) {
	//	n2.annotated = true;
	//}

	// If the nodes are fixed and too far apart, create a ghost node
	if (n_tert.source.fixed && n_tert.target.fixed) {
		var dx = n_tert.target.x - n.x,
	        dy = n_tert.target.y - n.y,
	        dr1 = Math.sqrt(dx * dx + dy * dy);
			dx = n_tert.source.x - n.x;
	        dy = n_tert.source.y - n.y;
	        dr2 = Math.sqrt(dx * dx + dy * dy);
		if (dr1 > 500 || dr2 > 500) {
			// Create ghost node for source
			n = this.createMetabolite(n.origin, true);
			n.colour = "#bbb";
			n.ghost = true;
			n.id += "_ghost_"+n_tert.id;

			// Does this ghost node have a layout cache?
			if (this.cached[n.id]) {
				n.fixed = true;
				n.fx = this.cached[n.id][0];
				n.fy = this.cached[n.id][1];
				n.x = n.fx;
				n.y = n.fy;
			} else {
				n.fixed = false;
				delete n.fx;
				delete n.fy;
			}
			if (nodes[n.id]) console.error("EXISTING GHOST", nodes[n.id]);
			nodes[n.id] = n;
		}
	}

	n.count++;
	n_tert.count++;
	//np2.count++;

	links.push({
		id: "SEC_"+secindex,
		origin: reaction,
		source: n,
		target: n_tert,
		master: master,
		input: !input,
		val: master.val, // + ((scaled < 0) ? -1 : 1),
		blocked: master.blocked,
		missing: master.missing,
		primary: false,
		secondary: true
	});
	secindex++;
}

MetabolicGraph.prototype.makeTertiary = function(reaction, p, m, nodes, links, master, input, tert) {
	let id = "TERT_" + master.id + "_" + ((input) ? "IN" : "OUT");

	if (!nodes[id]) nodes[id] = this.createMergedMetabolite(m, false, id);
	let n = nodes[id];
	n.id = id;
	n.tertiary = true;
	n.ghost = true;
	n.colour = "white";

	// Hack to get into good initial position;
	if (!n.fixed) {
		n.x = p.x;
		n.y = p.y;
		n.px = p.x;
		n.py = p.y;
	}

	n.count++;
	p.count++;

	links.push({
		id: "SEC_"+secindex,
		origin: reaction,
		source: n,
		target: p,
		master: master,
		input: !input,
		val: master.val, // + ((scaled < 0) ? -1 : 1),
		blocked: master.blocked,
		missing: master.missing,
		primary: false,
		secondary: true,
		tertiary: true
	});
	secindex++;

	/*links.push({
		id: "SEC_"+secindex,
		origin: reaction,
		target: n,
		source: np2,
		master: master,
		input: false,
		val: master.val, // + ((scaled < 0) ? -1 : 1),
		blocked: master.blocked,
		missing: master.missing,
		primary: false,
		secondary: false
	});
	secindex++;*/
	return n;
}

MetabolicGraph.prototype.setReactions = function(list) {
	let metabData = this.options.metaboliteData;
	let reactData = this.options.reactionData;
	let reactions = this.filterReactions(list, reactData);

	// Generate colours
	// Tally represented subsystems.
	let subsystems = {};
	for (var i=0; i<reactions.length; i++) {
		subsystems[reactions[i].subsystem] = true;
	}

	// Generate colours and positions for subsystems
	this.subsys_colours = {};
	let subsys_positions = {};
	let subsyscount = {};
	let sscount = 0;
	let sscount2 = 0;
	for (var x in subsystems) {
		sscount++;
	}
	for (var x in subsystems) {
		subsyscount[x] = 0;
		this.subsys_colours[x] = selectColor(sscount2,sscount);
		//subsys_positions[x] = rotate(1500,1500, 1900, 1500, ((2*Math.PI) / sscount) * sscount2);
		sscount2++;
	}

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

	var datascale = 20.0 / maxdata;
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
		let [prim_in,tert_in] = this.choosePrimary(reactions[i].inputs);
		let [prim_out,tert_out] = this.choosePrimary(reactions[i].outputs);

		if (prim_in.length == 0 || prim_out.length == 0) continue;

		let lenIn = (this.options.multipleLinks) ? prim_in.length : 1;
		let lenOut = (this.options.multipleLinks) ? prim_out.length : 1;
		let tertLenIn = (tert_in.length > 0) ? 1 : 0;
		let tertLenOut = (tert_out.length > 0) ? 1 : 0;

		// For each secondary, create links to each of the original primaries
		// Only one of the links is visualised, but it must know the other links target
		// Some secondaries and primary-level - use threshold. Still links like a secondary

		// 1) Create primary link
		// 2) For each secondary, create 2 links to each primary
		//       - include ref to master link
		//       - mark as level-2
		// 3) For each tertiary, create 2 links to each primary
		//       - mark as level-3
		// Ignore tertiary for now!!!

		let master = this.makePrimary(reactions[i], prim_in[0], prim_out[0], nodes, links, blocked, missing, scaled, data);
		if (!master) continue;

		// Step 1: create virtual tertiary node
		let n_tert = {
			type: "tertiary",
			id: "TERT_"+reactions[i].id,
			origin: reactions[i],
			invisible: true,
			name: reactions[i].name,
			fullname: reactions[i].name,
			fixed: true,
			target: master.target,
			source: master.source,
			count: 1
		};
		nodes[n_tert.id] = n_tert;
		this.tert_nodes.push(n_tert);

		// Step 2: Position node in approximate location
		let n1 = master.target;
		let n2 = master.source;
		if (n1.fixed && n2.fixed) {
			n_tert.x = n1.x + (n2.x-n1.x)/2;
			n_tert.y = n1.y + (n2.y-n1.y)/2;
			n_tert.fx = n_tert.x;
			n_tert.fy = n_tert.y;
		} else {
			n_tert.x = 1500;
			n_tert.y = 1500;
		}

		for (var j=1; j<lenIn; j++) {
			this.makeSecondary(reactions[i], n_tert, prim_in[j], nodes, links, master, true);
		}
		for (var j=1; j<lenOut; j++) {
			this.makeSecondary(reactions[i], n_tert, prim_out[j], nodes, links, master, false);
		}

		if (this.options.tertiary) {
			// Step 3: Create in and out tertiary nodes + links to virtual tert node
			let t1, t2;
			if (tertLenIn > 0) t1 = this.makeTertiary(reactions[i], n_tert, tert_in, nodes, links, master, true);
			if (tertLenOut > 0) t2 = this.makeTertiary(reactions[i], n_tert, tert_out, nodes, links, master, false);

			// Step 4: Link tertiary nodes primaries.
			if (t1) {
				links.push({
					id: "SEC_"+secindex++,
					origin: reactions[i],
					target: n2,
					source: t1,
					master: master,
					input: false,
					val: 0, // + ((scaled < 0) ? -1 : 1),
					blocked: true,
					missing: false,
					primary: false,
					secondary: true,
					tertiary: true,
					invisible: true
				});
			}

			if (t2) {
				links.push({
					id: "SEC_"+secindex++,
					origin: reactions[i],
					target: n1,
					source: t2,
					master: master,
					input: false,
					val: 0, // + ((scaled < 0) ? -1 : 1),
					blocked: true,
					missing: false,
					primary: false,
					secondary: true,
					tertiary: true,
					invisible: true
				});
			}
		}


		//for (var k=0; k<lenOut; k++) {
		/*let primary = j == 0 && k == 0;

		if (!nodes[prim_in[j].id2]) nodes[prim_in[j].id2] = this.createMetabolite(prim_in[j], j == 0);
		if (!nodes[prim_out[k].id2]) nodes[prim_out[k].id2] = this.createMetabolite(prim_out[k], k == 0);
		//links[reactions[i].id] = createReaction(reactions[i], prim_in, prim_out);

		//if (prim_in === prim_out) continue;

		if (this.options.hideLoose && !nodes[prim_in[j].id2].fixed && !nodes[prim_out[k].id2].fixed) continue;

		if (this.options.annotations && this.options.annotations[reactions[i].id]) {
			nodes[prim_out[k].id2].annotated = true;
		}

		let n1 = nodes[prim_in[j].id2];
		let n2 = nodes[prim_out[k].id2];*/

		// If the nodes are fixed and too far apart, create a ghost node
		/*if (n1.fixed && n2.fixed) {
			var dx = n1.x - n2.x,
		        dy = n1.y - n2.y,
		        dr = Math.sqrt(dx * dx + dy * dy);
			if (dr > 500) {
				// Create ghost node for source
				n1 = this.createMetabolite(prim_in[j], j == 0);
				n1.colour = "#bbb";
				n1.ghost = true;
				n1.id += "_ghost_"+reactions[i].id+"_"+n2.id;

				// Does this ghost node have a layout cache?
				if (this.cached[n1.id]) {
					n1.fixed = true;
					n1.fx = this.cached[n1.id][0];
					n1.fy = this.cached[n1.id][1];
					n1.x = n1.fx;
					n1.y = n1.fy;
				} else {
					n1.fixed = false;
					delete n1.fx;
					delete n2.fy;
				}
				if (nodes[n1.id]) console.error("EXISTING GHOST", nodes[n1.id]);
				nodes[n1.id] = n1;
			}
		}

		n1.count++;
		n2.count++;

		links.push({
			id: reactions[i].id,
			origin: reactions[i],
			source: n1,
			target: n2,
			input: false,
			val: scaled + ((scaled < 0) ? -1 : 1),
			blocked: blocked,
			missing: missing,
			primary: j == 0 && k == 0
		});

		if (reactions[i].reversible) {
			links.push({
				id: reactions[i].id+"_rev",
				origin: reactions[i],
				target: n1,
				source: n2,
				input: false,
				val: scaled + ((scaled < 0) ? -1 : 1),
				blocked: blocked,
				missing: missing,
				primary: j == 0 && k == 0
			});
		}*/

		//}
		//}
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
	if (this.options.hideLoose) {
		for (var x in nodes) {
			if (nodes[x].count == 0) delete nodes[x];
		}
	}

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
	this.graphData({nodes: nodes, links: links}, link => (link.tertiary) ? ((link.invisible) ? 20 : 15) : (link.secondary) ? 10 : 40, -100);
}

function selectColor(colorNum, colors){
	let cols = Math.ceil(colors / 2);
	let lit = (colorNum % 2 == 1) ? "50%" : "70%";
    if (colors < 1) colors = 1; // defaults to one color - avoid divide by zero
    return "hsl(" + Math.floor((colorNum / 2) * (360 / cols) % 360) + ",100%,"+lit+")";
}

function pathStyle(d) {
	let alpha = (d.primary) ? 1.0 : 0.5;
	if (d.invisible) return "stroke: none";

	if (d.secondary || d.tertiary) {
		if (d.missing) return "stroke: rgba(100,100,100,0.8); stroke-width: 1px;"; // stroke-dasharray: 2, 2;";
		if (d.blocked) return "stroke: rgba(0,0,255,"+alpha+"); stroke-width: 1px";
		if (d.val < 0) return "stroke-width: 1px; stroke: rgba(0,255,0,"+alpha+")";
		return "stroke-width: 1px; stroke: rgba(255,0,0,"+alpha+")";
	}

	//if (d.subsys && !d.sig) return "stroke: rgba(0,0,255,0.4); stroke-width: 1px;";
	if (d.missing) return "stroke: rgba(100,100,100,0.8); stroke-width: 1px;"; // stroke-dasharray: 2, 2;";
	if (d.blocked) return "stroke: rgba(0,0,255,"+alpha+"); stroke-width: 2px";

	//let val = (d.special) ? 0.08 : (Math.abs(d.val)/15 + 0.33);
	if (d.val < 0) return "stroke-width: " + Math.ceil(Math.abs(d.val)) + "px; stroke: rgba(0,255,0,"+alpha+")";
	return "stroke-width: " + Math.ceil(Math.abs(d.val)) + "px; stroke: rgba(255,0,0,"+alpha+")";
}

function nodeStyle(d) {
	return (d.type == "reaction") ? "fill: "+d.colour+";" + ((d.annotated) ? " stroke-width: 5px;" : "") : "";
}

/*MetabolicGraph.prototype.displaySubsystems = function(force, svg, colours) {
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
}*/

function angle(cx,cy, px,py) {
    var tx = cx;
	var ty = cy - Math.sqrt(Math.abs(px - cx) * Math.abs(px - cx)
            + Math.abs(py - cy) * Math.abs(py - cy));
    return (2 * Math.atan2(py - ty, px - tx));
}

function getPointAt(cx,cy, radius, angle) {
    return [cx + Math.sin(Math.PI - angle) * radius,
        cy + Math.cos(Math.PI - angle) * radius];
}

function angleDiff(a,b) {
	let c = a - b;
	if (Math.abs(c) > Math.PI) c = 2*Math.PI + c;
	return c;
}

// a = 5
// b = 355
// c -> 10
// c = -350

// a = 355
// b = 5
// c -> 10
// c = 350

MetabolicGraph.prototype.unfixNode = function(node) {
	d.fixed = false;
	delete d.fx;
	delete d.fy;
	delete me.cached[d.id];

	// TODO Now unfix all tertiary nodes.

	window.localStorage.gemvis_cached = JSON.stringify(me.cached, undefined, 2);
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
		.charge(n => (n.tertiary) ? charge : charge) // -300
		.linkStrength(link => (link.tertiary) ? ((link.invisible) ? 1 : 2) : (link.secondary) ? 1 : (link.missing) ? 0.3 : (link.blocked) ? 0.5 : 1)
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
		.data(["tiny","small","medium","large"])      // Different link/path types can be defined here
	  .enter().append("svg:marker")    // This section adds in the arrows
		.attr("id", String)
		.attr("viewBox", "0 -5 10 10")
		.attr("refX", d => (d == "tiny") ? 2 : (d == "small") ? 4 : (d == "medium") ? 6 : 10) //15
		.attr("refY", 0) //-1.5
		.attr("markerWidth", d => (d == "tiny") ? 4 : (d == "small") ? 8 : (d == "medium") ? 12 : 20)
		.attr("markerHeight", d => (d == "tiny") ? 4 : (d == "small") ? 8 : (d == "medium") ? 12 : 20)
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
		.attr("marker-mid", function(d) { return (d.secondary) ? "url(#tiny)" : (d.input) ? "" : (d.val <= 5) ? "url(#small)" : (d.val <= 9) ? "url(#medium)" : "url(#large)"; });

	let title = path.append("title");
	title.append("tspan").text(d => "Name: "+d.origin.name);
	title.append("tspan").text(d => "\nSystem: "+d.origin.subsystem);
	title.append("tspan").text(d => "\nData: "+d.data);
	title.append("tspan").text(d => "\nEC: "+ ((d.origin && d.origin.ec != "") ? d.origin.ec : "N/A"));

	var drag = force.drag().on("dragstart", function(d) {
		d.fixed = true;
		d.fx = d.x;
		d.fy = d.y;
		//d.colour = "#eee";
		//this.setAttribute("fill", "#eee");
	}).on("dragend", function(d) {
		me.cached[d.id] = [d.x,d.y];
		window.localStorage.gemvis_cached = JSON.stringify(me.cached, undefined, 2);
	});

	// define the nodes
	var node = g.selectAll(".node")
		.data(force.nodes())
	  .enter().append("g")
		.attr("class", function(d) { return ((d.special) ? "node special" : "node") + " " + d.type + ((d.blocked) ? " blocked" : "") + ((d.dead) ? " dead" : "") + ((d.val < 0) ? " negative" : ""); })
		.on("dblclick", function(d) {
			me.unfixNode(d);
		})		
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
		//.attr("font-size", d => (d.ghost) ? "6pt" : "8pt")
		.attr("class", d => (d.ghost) ? "ghost" : "")
		.attr("text-anchor","middle");
		//.text(d => d.name);

		t.selectAll("tspan")
		.data(function(d) { return (d.type == "metabolite") ? d.name.split("\n") : ""; })
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
		.attr("stroke", d => (d.annotated) ? "orange" : (d.ghost) ? "none" : "#222")
		.attr("stroke-width", d => (d.annotated) ? 3 : 2)
		//.attr("width", bb.width + paddingLR*2)
		//.attr("height", bb.height + paddingTB)
		.attr("fill", d => d.colour);

	title = r.append("title");
	title.append("tspan").text(d => "Name: "+d.fullname);
	title.append("tspan").text(d => "\nSystem: "+d.subsystem);
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
		// TODO Update virtual tertiary node.
		for (var i=0; i<me.tert_nodes.length; i++) {
			let d = me.tert_nodes[i];
			let x1 = d.target.x,
				x2 = d.source.x,
				y1 = d.target.y,
				y2 = d.source.y;
			let dx = d.target.x - d.source.x,
			    dy = d.target.y - d.source.y,
				dr = Math.sqrt(dx * dx + dy * dy);

			let x3 = (x1+x2) / 2;
			let y3 = (y1+y2) / 2;

			let cx = (false) ? x3 + Math.sqrt(dr*dr-Math.pow(dr/2,2))*(y1-y2)/dr : x3 - Math.sqrt(dr*dr-Math.pow(dr/2,2))*(y1-y2)/dr;
			let cy = (false) ? y3 + Math.sqrt(dr*dr-Math.pow(dr/2,2))*(x2-x1)/dr : y3 - Math.sqrt(dr*dr-Math.pow(dr/2,2))*(x2-x1)/dr;

			let a1 = angle(cx,cy,x1,y1);
			let a2 = angle(cx,cy,x2,y2);
			let a = (angleDiff(a2, a1) / 2) + a1;

			let [sx,sy] = getPointAt(cx,cy,dr,a);

			d.x = sx;
			d.y = sy;
			d.fx = sx;
			d.fy = sy;
			d.centre_x = cx;
			d.centre_y = cy;
			d.radius = dr;
		}

		path.attr("d", function(d) {
			if (d.secondary) {
				let dx2 = d.target.x - d.source.x,
					dy2 = d.target.y - d.source.y,
					dr2 = Math.sqrt(dx2 * dx2 + dy2 * dy2)*0.8;

				let dcx = d.source.x - d.target.centre_x,
					dcy = d.source.y - d.target.centre_y,
					dcr = Math.sqrt(dcx*dcx+dcy*dcy);
				let inside = dcr <= d.target.radius;

				if (isNaN(dr2)) return "";

				if (!d.input) {
					return "M" + 
						d.source.x + "," + 
						d.source.y + "A" + 
						dr2 + "," + dr2 + " 0 0," + ((!inside) ? "0" : "1") + " " + 
						d.target.x + "," + 
						d.target.y;
				} else {
					return "M" + 
						d.target.x + "," + 
						d.target.y + "A" + 
						dr2 + "," + dr2 + " 0 0," + ((inside) ? "1" : "0") + " " + 
						d.source.x + "," + 
						d.source.y;
				}
			} else if (d.primary) {
				var dx = d.target.x - d.source.x,
				    dy = d.target.y - d.source.y,
				    dr = Math.sqrt(dx * dx + dy * dy);

				return "M" + 
					d.source.x + "," + 
					d.source.y + "A" + 
					dr + "," + dr + " 0 0," + ((d.input) ? "0" : "1") + " " + 
					d.target.x + "," + 
					d.target.y;
			}
		});

		t.each(function(d) {
			if (d.invisible) return;
			let pLR = (d.ghost) ? paddingLR*0.5 : (d.tertiary) ? 0 : paddingLR;
			let pTB = (d.ghost) ? paddingTB*0.5 : (d.tertiary) ? 0 : paddingTB;

			let bb = this.getBBox();
			let r = this.previousSibling;

			this.setAttribute("y", -bb.height/2 - pTB);

			r.setAttribute("x", - bb.width/2 - pLR);
			r.setAttribute("y", -bb.height/2 - pTB);
			r.setAttribute("width", bb.width + pLR*2);
			r.setAttribute("height", bb.height + 2*pTB);
		});

		node
		    .attr("transform", function(d) { 
	  	    return (isNaN(d.x)) ? "" : "translate(" + d.x + "," + d.y + ")"; });
	}
}

module.exports = MetabolicGraph;


