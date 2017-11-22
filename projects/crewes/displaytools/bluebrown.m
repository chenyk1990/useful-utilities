function c = bluebrown(m)
%BLUEBROWN Linear color map grading from light blue through dark brown
%   BLUEBROWN(M) returns an M-by-3 matrix containing a the colormap.
%   BLUEBROWN, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(BLUEBROWN)
%
%   This uses COPPER for the first half of the colormap and BLUEBLACK for the second half
%
% G.F. Margrave, Devon Energy, 2017


if nargin < 1, m = size(get(gcf,'colormap'),1); end

c1=copper(m);
c2=blueblack(m);
m2=floor(m/2);

c = [c1(1:m2,:); c2(m2+1:m,:)];