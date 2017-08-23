#! /bin/sh

# estmiate source wavelet by shot and receiver response 
# and deconvolve input data

set -x 

infile=multiple_suppressed_muted_data.su
outfile=shot_receiver_sigdecon_$infile

echo "Signature decon with homomorphic wavelet estimation "

## Assumptions:
## 1) wavelet is constant within a shot gather
## 2) reflectivity and multiple series are random
## 3) wavelet is minimum phase  

## correcting for souce wavelet
# remove output files
rm shot_$outfile
rm $outfile
rm all_shot_wavelets.su
rm all_receiver_responses.su

# split the original data into shot gathers
# split shot data
susort  sx offset < $infile > shot_gathers.su
susplit < shot_gathers.su key=sx


# loop over shot gather files
for i in `ls split_sx* `
do

	# stack real part of the complex log transform (log|A|)
	suclogfft < $i  | suamp mode=real | sustack key=dt > real.su

	# stack imaginary part of the complex log transform (phase)
	suclogfft < $i  | suamp mode=imag | sustack key=dt > imag.su

	# combine real and imag and inverse clog transform
        # and output a minimum phase version of
	suop2 real.su imag.su op=zipper | suiclogfft  |
	suminphase | suwind itmax=199 > wavelet_est.su

	cat wavelet_est.su >> all_wavelets.su

	# deconvolve with sucddecon
	sucddecon sufile=wavelet_est.su < $i >> shot_$outfile


done

rm split_sx*

exit 0
