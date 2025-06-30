reset
set terminal postscript eps color enhanced "Times-Roman, 22"
set output 'cp_avg.eps' 
set size 1,1.2
set ylabel "C_P" offset 3,0,0
set xlabel "{/Symbol q} (deg)"
plot[0:180][] "cp_re100_nc32.dat" u 1:2 with l lw 3 lc -1 title 'Present (strong conserv)',\
     "lit_hwang_20072.dat" u 1:2 w p pt 6 ps 2 title 'Hwang and Yang, 2007', \
     "sharman_2005.dat" u 1:2 w p pt 8 ps 2 lc 7 title 'Sharman, et. al., 2005'

reset
set terminal postscript eps color enhanced "Times-Roman, 22"
set output 'cp_smooth_avg.eps' 
set size 1,1.2
#set key box opaque
set ylabel "C_P" offset 3,0,0
set xlabel "{/Symbol q} (deg)"
plot[0:180][] "cp_re100_nc32.dat" u 1:3 with l lw 3 lc -1 title 'Smooth (strong conserv)',\
     "lit_hwang_20072.dat" u 1:2 w p pt 6 ps 2 title 'Hwang and Yang, 2007', \
     "sharman_2005.dat" u 1:2 w p pt 8 ps 2 lc 7 title 'Sharman, et. al., 2005'






