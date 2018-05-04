/**
 * This script runs a simple webserver which is capable of passing through
 * requests to various ruby scripts we use. It's intended as an alternative to
 * running a vagrant box with apache.
 *
 * While running this a jseden instance should be visible at localhost:8000
 */
var finalhandler = require('finalhandler')
var http = require('http')
var serveStatic = require('serve-static');
var httpProxy = require('http-proxy');

var proxy = httpProxy.createProxyServer({});
var serve = serveStatic('./build');

var port = 8000;

var server = http.createServer(function onRequest(req, res) {
	if (req.url.startsWith("/KEGG/")) {
		proxy.web(req, res, { ignorePath: true, target: 'http://rest.kegg.jp/'+req.url.substring(6) });
	} else {
		serve(req, res, finalhandler(req,res));
	}
});

server.listen(port);

