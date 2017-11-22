function figtit(str,figno)
% Put a title in the title bar of a Figure
%
% figtit(str,figno)
% 
% str ... string giving figure title
% figno ... figure number
% ******** default = gcf *********
%


if(nargin<2)
figno=gcf;
end

set(figno,'name',str)