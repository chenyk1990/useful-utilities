function figfwh2(figno)
%
% figfwh2(figno)
%
% Makes the figure figno full width of the screen and half the height while
% allowing for a 75 pixel 'border' all around. If no input, then a new
% blank figure is created with the mentioned size.
%
if(nargin<1)
    figure;
    figno=gcf;
end
ss=get(0,'screensize');
border=75;
set(figno,'position',[border, border, ss(3)-2*border round(ss(4)/2)-2*border]);