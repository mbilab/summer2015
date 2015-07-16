require! <[fs]>

(err, data) <-! fs.read-file \test encoding: \utf-8
throw err if err
layout = data.match /<div id=layout.*div>/
console.log layout

# vi:et:sw=2:ts=2:fdm=indent:ft=ls:nowrap
