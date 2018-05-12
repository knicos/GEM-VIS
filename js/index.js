exports.MetaboliteGraph = require('./metabolite-graph.js');
exports.ReactionGraph = require('./reaction-graph.js');
exports.Dialog = require('./dialog.js');

function graph(element, model, options) {
	if (!options) {
		console.error("Graph options required");
		return null;
	}
	if (!model) {
		console.error("Metabolic model required");
		return null;
	}

	if (options.style == "metabolite") return new exports.MetaboliteGraph(element, model, options);
	if (options.style == "reaction") return new exports.ReactionGraph(element, model, options);
	return null;
}

exports.graph = graph;

