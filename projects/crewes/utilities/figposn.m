function figposn(x,y,figno)
% FIGPOSN: move a figure to a specific place on the screen based on screensize
%
% figposn(x,y,figno)
%
% x ... x position of the figure center in normalized screen coordinates. 0.5 means the center of
%       the screen while 0 is the left-hand side of the screen and 1 is the right-hand side.
% y ... y position of the figure center in normalized screen coordinates. 0.5 means the center of
%       the screen while 0 is the bottom of the screen and 1 is the top side.
% figno ... handle of the figure to be so moved
% *********** default is gcf (the current figure)***********
%
if(nargin<3)
    figno=gcf;
end

ssize=get(0,'screensize');

pos=get(figno,'position');
w=pos(3);
h=pos(4);
x0=pos(1);
y0=pos(2);
x1=x*ssize(3);
y1=y*ssize(4);
dx=x1-x0-w/2;
dy=y1-y0-h/2;
set(figno,'position',[x0+dx y0+dy w h]);