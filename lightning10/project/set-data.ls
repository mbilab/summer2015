require! <[fs tabletojson]>

(err, data) <-! fs.read-file \test encoding: \utf-8 # read file
throw err if err
table-as-json = tabletojson.convert data # convert table to json

output =
  home: {}
  away: {}
  game-info: {}

set-json!

!function set-json
  get-inning table-as-json[13]
  get-hitter \home table-as-json[19]
  get-hitter \away table-as-json[18]
  get-pitcher \home table-as-json[22]
  get-pitcher \away table-as-json[21]
  get-game-info table-as-json[23]
  build-json output

!function get-inning arr
  home = output.home.inning? = {}
  away = output.away.inning? = {}
  for k, v of arr[2]
    d = if k isnt '0' => v else \name
    home[d] = arr[4][k]
    away[d] = arr[3][k]

!function get-hitter team, arr
  _arr = output[team].hitter? = []
  for d, i in arr
    continue if i is 0
    return if d.0 is \Totals
    _arr.push d

!function get-pitcher team, arr
  _arr = output[team].pitcher? = []
  for d, i in arr
    continue if i is 0
    return if key-count(d) is 1
    _arr.push d

!function get-game-info arr
  for d, i in arr
    continue if i is 0
    return if i is 4
    _arr = d.0.split \-
    output.game-info[_arr[0]] = _arr[1]

!function build-json obj
  j = JSON.stringify obj, null, 2
  fs.write-file-sync \output.json j

function key-count obj
  (Object.keys obj).length

# vi:et:sw=2:ts=2:fdm=indent:ft=ls:nowrap
