function c = blueblack(m)
%BLUEBLACK Linear blue-tone color map
%   BLUEBLACK(M) returns an M-by-3 matrix containing a "blueish" colormap.
%   BLUEBLACK, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(BLUEBLACK)
%
%   Modelled after Matlab's COPPER colormap.
%
% G.F. Margrave, Devon Energy, 2017


if nargin < 1, m = size(get(gcf,'colormap'),1); end
c = min(1,gray(m)*diag([0.4975 0.7812 1.2500]));
