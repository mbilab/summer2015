#! /usr/bin/node
var fs = require('fs');
var sh = require('sh');

// read source files
fs.readdir('games/', function (err, files) {
	if (err) throw err;
	//console.log(files);
	for(i in files){
		console.log(files[i]);
		var files1 = fs.readdirSync('games/'+files[i]+'/'); 
		for(j in files1){
		
		// Do parse work for each html
			console.log(files1[j]);








		}
	}
});
