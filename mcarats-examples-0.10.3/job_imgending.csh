#!/bin/csh -f

# Arguments
if ($#argv < 2) then
    echo "usage: $0 lastFilePrefix suffix"
    exit 2
endif
set pfx=$argv[1]
set sfx=$argv[2]


# Black out
foreach p(1 2 3 4 5 6 7 8 9)
    cp $pfx.$sfx $pfx"e0"$p.$sfx
end

foreach p(0 1 2 3 4 5 6 7 8 9)
    convert $pfx.$sfx -fill black -colorize $p"0%" $pfx"e1"$p.$sfx
end
convert $pfx.$sfx -fill black -colorize 99% $pfx"e20".$sfx

# Credit
convert -font helvetica -fill white -pointsize 24 -draw 'text 60,100 "Virtual Flight in Artificial Cloudy Atmosphere"' -pointsize 14 -draw 'text 90,200 "Rendered: June, 2009, on MacBook with Intel Core 2 Duo 2.4 GHz"' -pointsize 14 -draw 'text 90,240 "Using MCARaTS version 0.10 alpha"' -pointsize 14 -draw 'text 90,280 "Copyright (C) 2009 Hironobu Iwabuchi"' -pointsize 14 -draw 'text 90,320 "MCARaTS, a free software distributed under GNU GPL v3"' $pfx"e18.$sfx" $pfx"e30.$sfx"

# Black to credit
foreach p(21 22 23 24 25 26 27 28 29)
    @ pct = 30 - $p
    convert $pfx"e30.$sfx" -fill black -colorize $pct"0%" $pfx"e"$p.$sfx
end
