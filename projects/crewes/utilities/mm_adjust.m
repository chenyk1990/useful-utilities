function mm_adjust(fig, parent)
% function mm_adjust(handle, parent)
%
% Where:
%   fig    = handle to figure to adjust
%   parent = handle to parent figure
%
% 1) Figure out which monitor parent figure is on
% 2) Move fig to the center of that monitor
%
% Examples:
%   mm_adjust(plot(sin(1:10)),parent);
%  
%  Authors: Kevin Hall, 2017
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

narginchk(1,2);

if nargin <2
    parent = [];
end

set(0,'Units','pixels');
mps = get(0,'MonitorPositions');
nmonitors = size(mps,1);

if nmonitors <2
    %never mind
    return
end

funits = get(fig,'Units');
set(fig,'Units','pixels');

if isempty(parent)
    punits = get(0,'Units'); %get groot units
    pps = get(0,'ScreenSize'); %get default monitor screensize
else    
    punits = get(parent,'Units');
    set(parent,'Units','pixels');
    pps = get(parent,'Position'); %get parent position
end

mcenter = mps(:,1:2) +mps(:,3:4)./2; %calculate x,y for center of each monitor
pcenter = pps(1:2) +pps(3:4)./2; %calculate x,y for center of parent figure

dxdy = mcenter-pcenter; %find dx and dy
pmdist = sqrt(dxdy(:,1).^2 +dxdy(:,2).^2); %calculate distance between pcenter and mcenter

mcenter = mcenter(pmdist == min(pmdist),:); %x,y of monitor center parent figure center is closest to

fps = get(fig,'Position'); %get figure position
fps(1) = mcenter(1) - fps(3)/2; %update fig x
fps(2) = mcenter(2) - fps(4)/2; %update fig x
set(fig,'Position',fps); %set updated figure position

set(fig,'Units',funits);
if isempty(parent)
    set(0,'Units',punits);
else
    set(parent,'Units',punits);
end

end %end function
