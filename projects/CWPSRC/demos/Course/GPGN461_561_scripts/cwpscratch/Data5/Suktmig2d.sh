#! /bin/sh

#
vfile=vel_nmo.bin      # rms velocity as a function of cdp and time
dx=25                  # spacing between receivers
indata=radon_mute_filtered_repaired_co.su
outdata=ktmig.su

rm ktmig.su

# split the data into a bunch of common offset gaters
susplit < $indata key=offset




# loop over shot gather files
for i in `ls split_* `
do
	suktmig2d vfile=$vfile dx=$dx  < $i >> $outdata
done

## clean up
# remove shot split files
rm split*


exit 0
