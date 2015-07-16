require! <[fs]>

(err, data) <-! fs.read-file \index.html encoding: \utf-8
throw err if err
for mat in data.match /<h3.*?h3>/g => console.log mat - /(<.*?>)/g if mat isnt /<h3.*?normal.*>/

# vi:et:fdm=indent:sw=2:ts=2:ft=ls:nowrap
