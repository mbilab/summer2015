require! <[fs]>
opt =
    games: \./games/
cat-dir = fs.readdir-sync \./games
for date-dir in cat-dir
    dir = opt.games + date-dir
    games = fs.readdir-sync "#{dir}/"
    for game in games
        game-path = dir + "/" + game
        tmp = fs.read-file-sync game-path, encoding: \utf-8
