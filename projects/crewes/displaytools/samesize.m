function samesize(hfig1,hfig2)
% SAMESIZE: make figure 1 the same size as figure 2
%
% samesize(hfig1,hfig2)
%

%pos1=get(hfig1,'position');
pos2=get(hfig2,'position');
set(hfig1,'position',[1 1 pos2(3:4)])