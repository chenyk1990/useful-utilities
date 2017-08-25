./iter_decon -F1/3/-5 -N100 -C-2/-10/80 -T0.1 example/KUL.z example/KUL.[r,t]
#SAC> r ./example/KUL.* ./example/pp.*
#SAC> w over

./k_stack -R20/60/1.5/2.0 -I0.5/0.01 -Gexample/hk.grd example/pp.*.ri
./grdmin -D example/hk.grd
./grdmin -D example/hk.grd | ./hk_plt.pl > junk.ps

open junk.ps

