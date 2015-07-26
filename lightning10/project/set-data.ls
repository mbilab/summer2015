require! <[fs sh]>
html = require \fast-html-parser

path =
  src: \./src
  res: \./res

sh('ls ' + path.src).result (r) !-> # get all days
  days = r.split '\n'
  days.pop!
  Do days

!function Do days
  return if days.length is 0
  day = days.shift!
  fs.mkdir-sync(path.res + "/#day") # make dir for each day
  (r) <-! sh('ls ' + path.src + "/#day").result
  games = r.split '\n'
  games.pop! # get all games of this day
  <-! build-json day, games
  Do days

!function build-json day, games, cb
  return cb! if games.length is 0
  game = games.shift!
  json = home: {} away: {} info: {}
  (root) <-! get-data day, game
  <-! set root, json
  output day, game, json
  build-json day, games, cb

!function get-data day, game, cb
  (err, data) <-! fs.read-file (path.src + "/#day/#game"), encoding: \utf-8
  throw err if err
  cb html.parse data

!function set data, obj, cb
  <-! set-score data, obj
  return cb! if it?
  <-! set-data data, obj
  cb!

!function set-score data, obj, cb
  period = data.query-selector '.lineScore .periodLabels'
  return cb {+err} if! period? # no game data
  home = data.query-selector '.lineScore .homeTeam'
  away = data.query-selector '.lineScore .awayTeam'
  inning = []
  for i in period.child-nodes => inning.push i.text
  score-map home, \home
  score-map away, \away
  cb!

  !function score-map item, type
    name = item.query-selector '.teamLocation a'
    _obj = obj[type].inning? = name: name.text
    for i, j in item.child-nodes => _obj[inning[j]] = i.text if j isnt 0

!function set-data data, obj, cb
  labels = data.query-selector-all '.data .label'
  return cb! if labels.length is 0 # no game data
  p = data.query-selector-all '.data'
  players = p.filter -> it.class-names.length is 1
  hitter = get-label labels[0]
  pitcher = get-label labels[2]
  set-hitter players[0], \away
  set-hitter players[1], \home
  set-pitcher players[2], \away
  set-pitcher players[3], \home
  set-game-info players[4]
  cb!

  function get-label o, label = []
    for i, j in o.child-nodes => label.push i.text
    label

  !function set-hitter o, type
    hitters = o.query-selector-all \#batter
    _arr = obj[type].hitters? = []
    for h in hitters
      _obj = name: h.child-nodes[0].text - /^\s+/
      for i, j in h.child-nodes => _obj[hitter[j]] = i.text if j isnt 0
      _arr.push _obj

  !function set-pitcher o, type
    pitchers = o.child-nodes
    _arr = obj[type].pitchers? = []
    for p in pitchers
      continue if! p.attributes.align?
      break if p.attributes.id is \footnote
      _obj = name: p.child-nodes[0].text - /^\s+/
      for i, j in p.child-nodes => _obj[pitcher[j]] = i.text if j isnt 0
      _arr.push _obj

  !function set-game-info o
    info = o.query-selector \#footnote
    return if! info?
    arr = info.structured-text.split '\n'
    for i til 3
      _arr = arr[i].split \-
      obj.info[_arr[0]] = _arr[1]
    set-game-umpires arr[3]

  !function set-game-umpires s
    arr = s.split \,
    for i, j in arr
      _arr = i.split \-
      _obj = obj.info[_arr.shift!]? = {} if j is 0
      _obj[_arr[0]] = _arr[1]

!function output day, game, obj
  j = JSON.stringify obj, null, 2
  p = path.res + "/#day/#game.json"
  fs.write-file-sync p, j

# vi:et:sw=2:ts=2:fdm=indent:ft=ls:nowrap
