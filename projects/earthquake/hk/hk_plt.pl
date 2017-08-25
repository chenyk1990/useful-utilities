#!/usr/bin/env perl
#
# plot HK stacking result
#   input hk-stacking results from stdin
#

while(<STDIN>) {
     ($grd, $h, $kapa, $cxx, $cyy, $cxy, $min, $dev) = split;
     # plot H-k stack
     open(AAA,">junk.cpt");
	   printf AAA "%5.3f 100 100 100 0 255 255 255\n",$min;
     close(AAA);
     system("grdimage -JX4i/2i $grd -Cjunk.cpt -Ba5f1:\"H (km)\":/a0.1f0.02:k:WSne -K");
     $hk = `grdinfo $grd | grep Command | cut -d' ' -f 4`; chop($hk);
     ($aa,$minH,$maxH,$mink,$maxk) = split(/[R\/]/,$hk);
     $ratioH=4/($maxH-$minH);$ratiok=2/($maxk-$mink);
     $cxx = $cxx/($ratioH*$ratioH);
     $cyy = $cyy/($ratiok*$ratiok);
     $cxy = $cxy/($ratioH*$ratiok);
     $s1 = 0.5*($cxx-$cyy); $s2 = sqrt($s1*$s1+$cxy*$cxy);
     $theta = 0.5*atan2($cxy,$s1)*180./3.14159;
     $s1 = 0.5*($cxx+$cyy) + $s2;
     $s2 = 0.5*($cxx+$cyy) - $s2;
     if ($s2<0.) {print STDERR "Warning: $s2<0\n";$s2=$s1*1.e-6;}
     open(PLT,"|psxy -JX $hk -O -K -Se -W2/255/255/255");
     printf PLT "%6.2f %8.4f %6.1f %s %s\n",$h,$kapa,$theta,sqrt($dev/$s1),sqrt($dev/$s2);
     close(PLT);
     open(PLT,"|pstext -JX -R0/1/0/1 -O");
     printf PLT "0.65 0.85 14 0 0 1 %4.1f km/%4.2f\n",$h,$kapa;
     close(PLT);
}

exit(0);
