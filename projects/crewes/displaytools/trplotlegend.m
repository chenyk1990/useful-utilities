function trplotlegend(htraces,names)
% TRPLOTLEGEND ... provides a simple interface to legend for use with trplot
%
% trplotlegend(htraces,names)
%
% htraces ... this is the first return from trplot. It is a cell array with one entry per trace
%       plotted. Each cell can have multiple handles as determined by trplot
% names ... a cell array of names, the same length as htraces
%

if(~iscell(htraces))
    error('htraces must be a cell array')
end

if(~iscell(names))
    error('names must be a cell array')
end

if(length(names)~=length(htraces))
    error('htraces and names must be cell arrays of the same length');
end

if(~isgraphics(htraces{1}(1)))
    error('htraces must contain graphics handles')
end

if(~ischar(names{1}))
    error('names must contain strings')
end

nt=length(htraces);
hh=zeros(1,nt);
for k=1:nt
    hh(k)=htraces{k}(1);
end

legend(hh,names);