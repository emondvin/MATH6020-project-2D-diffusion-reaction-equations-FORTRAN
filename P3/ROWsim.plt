set pm3d
set size square  
unset surf
set cbr [0:1]
set zr [0:1]
set terminal gif animate delay 1
set output 'ROWani.gif'
do for [i=1:76] {splot sprintf('proj.%d',i); pause 0.1}
set output
set term wxt
