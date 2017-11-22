function ht=titlein(str,pos,fudge)
% TITLEIN: titles an axis with the title inside the axis boundaries
%
% titlein(str,pos,fudge)
%
% str ... title string or cell array
% pos ... either 't' which puts the title at the top or 'b' which puts it at the bottom, or 'c'
%       which puts it in the middle
% ********** default 't' ********
% fudge ... if the y axis extends from y1 to y2 (y1<y2), let dy =y2-y1, then 't' will cause the
%           title to be at y2-fudge*dy and 'b' will put the title at y1+fudge*dy. Fudge has no
%           effect for 'c'
% ********** default dy=.25 *********
%

if(nargin<3)
    fudge=.25;
end
if(nargin<2)
    pos='t';
end

if(pos~='t' && pos~= 'b' && pos~='c')
    error('pos must be either ''t'' or ''b'' or ''c'' ');
end

yl=get(gca,'ylim');
xl=get(gca,'xlim');
xm=mean(xl);
dy=diff(yl);

switch pos
    case 't'
        ht=title(str,'position',[xm yl(2)-fudge*dy]);
    case 'b'
        ht=title(str,'position',[xm yl(1)+fudge*dy]);
    case 'c'
        ht=title(str,'position',[xm mean(yl) 1000*eps]);
end
        