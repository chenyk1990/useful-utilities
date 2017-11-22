function c = redblack(m)
%COPPER Linear copper-tone color map
%   COPPER(M) returns an M-by-3 matrix containing a "copper" colormap.
%   COPPER, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(copper)
%
%   See also HSV, GRAY, HOT, COOL, BONE, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   C. Moler, 8-17-88, 8-19-92.
%   Copyright 1984-2004 The MathWorks, Inc.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
c = min(1,gray(m)*diag([1.500 .5 .8]));