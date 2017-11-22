function [figpos,monpos] = mm_figpos(h)
%function [figpos,monpos] = mm_figpos(h)
%
% Where:
%    h = handle to a figure; if empty, figure_position returns the position
%        of the Matlab Desktop
%    figpos = [x,y,width,height] of the Figure OR Matlab desktop
%    monpos = [x,y,width,height] of the Monitor that the Matlab figure or 
%        desktop is currently displayed on. This is based on the location 
%        of the center of the figure or desktop
%
%    x, y, width and height are in pixels
%
% Examples:
%   [f,m] = figure_position              %Matlab Desktop
%   [f,m] = figure_position([])          %Matlab Desktop
%   h=figure; [f,m] = figure_position(h) %Matlab Figure
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

narginchk(0,1)

if nargin==0 || isempty(h)
    figflag = 0;
    %get location of MATLAB desktop window; Requires java
    %NOTE: Origin is top left
    desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    
%     desktop.getMainFrame
    x = desktop.getMainFrame.getX;
    y = desktop.getMainFrame.getY;
    width  = desktop.getMainFrame.getWidth;
    height = desktop.getMainFrame.getHeight;
    figpos = [x,y,width,height] %In pixels
    
    %default windows screen resolution
    winres = 96;
    %actual screen resolution (high-res displays are higher than 96
    screenres = java.awt.Toolkit.getDefaultToolkit.getScreenResolution;
    %ratio
    scalefact = screenres/winres
    %adjust figpos to device independent pixels as used by Matlab by
    %default
    figpos = figpos./scalefact    
else
    figflag = 1;    
    %get location of Matlab figure
    %NOTE: Origin is bottom left
    tmp = get(h,'Units');
    set(h,'Units','pixels');
    figpos = get(h,'Position');
    x = figpos(1);
    y = figpos(2);
    width = figpos(3);
    height = figpos(4);
    set(h,'Units',tmp);
end
 
%Get monitor information
MonitorPosition = get(0,'MonitorPositions');
numMonitors = size(MonitorPosition,1);

if numMonitors<2
    monpos = MonitorPosition;
    
    return
end
%disp (['Number of monitors = ' num2str(numMonitors)]);
 
for monitor=1:numMonitors
    MonitorWidth   = MonitorPosition(monitor,3);
    MonitorHeight  = MonitorPosition(monitor,4);
    MonitorLeftX   = MonitorPosition(monitor,1);
    MonitorRightX  = MonitorPosition(monitor,1)+MonitorWidth;
    MonitorTopY    = MonitorPosition(monitor,2);
    MonitorBottomY = MonitorPosition(monitor,2)+MonitorHeight;    

    %Find integer x and y in pixels closest to center of desktop
    FigCenterX = round(x +width/2.);
    FigCenterY = round(y +height/2.); 
    
    if figflag %Matlab Figure; convert to Java origin
        FigCenterY=MonitorHeight-FigCenterY;
    end
    
    if FigCenterX >MonitorLeftX && FigCenterX <MonitorRightX && ...
            FigCenterY >MonitorTopY && FigCenterY <MonitorBottomY
            
                monpos = MonitorPosition(monitor,:);
                break
    end    
end

%return Matlab desktop position in matlab coordinate system 
%(origin is bottom left; up until now it was top left
if ~figflag %Matlab Desktop
    figpos(2) = MonitorHeight-y-height;
end

end %end function