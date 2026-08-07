set logscale
set xrange [50:70000]
set key top left

linear(x) = 2*x 
lin_log(x) = 9*x*log(x)

pl 'morton_AoS.dat' u 1:2 w lp lc 1 lw 2 dt 2 t 'AoS build',\
   'morton_AoS.dat' u 1:3 w lp lc 1 lw 2 dt 3 t 'AoS computation',\
   'morton_AoS.dat' u 1:4 w lp lc 1 lw 2 t 'AoS total',\
   'morton_SoA.dat' u 1:2 w lp lc 2 lw 2 dt 2 t 'SoA build',\
   'morton_SoA.dat' u 1:3 w lp lc 2 lw 2 dt 3 t 'SoA computation',\
   'morton_SoA.dat' u 1:4 w lp lc 2 lw 2 t 'SoA total',\
    [60:60000] linear(x) lw 2 t 'O(n)',\
    [60:60000] lin_log(x) lw 2 t 'O(nlog(n))'
    