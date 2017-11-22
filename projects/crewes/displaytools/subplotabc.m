function haxe=subplotabc(arg)
%
% haxe=subplotabc(pane)
%
% breaks figure window into three panes. Pane 'top' is the largest and
% takes up aboput 60% of the figure. It extends the full width and 60% of
% the height. Pane 'bota' is the second largest and is on the bottom left.
% Pane 'botb' is the smallest and is bottom right. All panes resize with
% the figure automatically. Run subplotabc with no inputs for a demo.
%

if(nargin<1)
    arg='demo';
end

if(strcmp(arg,'top'))
    haxe=subplot('position',[.1 .55 .8 .35]);
elseif(strcmp(arg,'bota'))
    haxe=subplot('position',[.1 .1 .5 .35]);
elseif(strcmp(arg,'botb'))
    haxe=subplot('position',[.65 .1 .25 .35]);
elseif(strcmp(arg,'demo'))
    figure
    subplotabc('top');title('top pane')
    subplotabc('bota');title('bota pane')
    subplotabc('botb');title('botb pane')
else
    haxe=[];
end