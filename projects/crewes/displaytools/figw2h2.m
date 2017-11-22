function figw2h2(figno)
%
% fig2w32h3(figno)
%
% Makes the figure figno 2/3 width of the screen and 2/3 the height while
% allowing for a 50 pixel 'border' all around. If no input, then a new
% blank figure is created with the mentioned size.
%
if(nargin<1)
    figure;
    figno=gcf;
end
ss=get(0,'screensize');
border=50;
set(figno,'position',[border, border, round(ss(3)/2)-2*border, round(ss(4)/2)-2*border]);