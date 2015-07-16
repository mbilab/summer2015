require! <[request]>
opt =
    url: \http://zoro.ee.ncku.edu.tw/wp2013/
(err, res, body) <- request.post opt.url
throw err if err
re = new Reg-exp "<h3.*>(.*)<\/h3>\n.*<p", "g"
while (arr = re.exec body) isnt null
    console.log arr[1]
