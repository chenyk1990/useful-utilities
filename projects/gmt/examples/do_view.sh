#!/bin/bash
#
#	$Id: do_view.sh 9857 2012-03-13 10:55:26Z fwobbe $
#
#	Simple driver to view all examples using ghostview
#
viewer=${1:-gv}
for f in example_*.ps
do
	$viewer $f
done
