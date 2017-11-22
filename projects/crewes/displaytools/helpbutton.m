function helpbutton(figno,helpmsg,titlestring,fontsize,windowsize,initialstate)
% HELPBUTTON: add a help button to a window (often a plotimage window)
%
% helpbutton(figno,helpmsg,titlestring,fontsize,windowsize,initialstate)
%
% figno ... number of the figure that will have the button installed
% helpmsg ... string containing the help message. Can be anylength and will
%       automatically wrap.
% titlestring ... title string for the help window
% fontsize ... fontsize to use in displaying the help message
% windowsize ... two element ve3ctor giving the size of the help window
%       (width,height) in terms of screen dimensions. [1,1' means it will
%       be the same size as the screen. [.1,.2] meanse its height will be
%       one tenth of the screen height and its width will be two tenths.
% ***************** default [.1 .1] *************
% initialstate ... 1 means the help window is created when the button is
%       created. 0 means the window dows not appear until the button is
%       pushed.
% ***************** default =0 ***************
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

% Add a help button to a plotimage figure
if(isgraphics(figno))
    action='init';
else
    action=figno;
    init=0;
end

if(strcmp(action,'init'))
    %check for existing buttons
    hkids=allchild(figno);
    startpos=.02;
    for k=1:length(hkids)
        tag=get(hkids(k),'tag');
        if(strcmp(tag,'raymigbutton')||strcmp(tag,'clearraysbutton')||...
                strcmp(tag,'clearpicksbutton')||strcmp(tag,'raymodbutton'))
            pos=get(hkids(k),'position');
            rightedge=pos(1)+pos(3);
            if(rightedge>startpos)
                startpos=rightedge+.01;
            end
        end
    end
    wpixels=60;%desired width in pixels
    %determine window size
    pos=get(figno,'position');%should be in pixels
    %
    width=wpixels/pos(3);%width in normaliized units
    uicontrol(figno,'style','pushbutton','string','Help!','units','normalized',...
        'position',[startpos .95 width .05],'callback','helpbutton(''help'')',...
        'tag','helpbutton','userdata',{helpmsg,titlestring,fontsize,windowsize},'tag','helpbutton');
    init=1;
    if(initialstate==1)
        action='help';
    end
end
if(strcmp(action,'help'))
    if(~init)
        hbutt=findobj(gcf,'tag','helpbutton');
        udat=get(hbutt,'userdata');
        helpmsg=udat{1};
        titlestring=udat{2};
        fontsize=udat{3};
        windowsize=udat{4};
        figno=get(hbutt,'parent');
    end
    pos=get(figno,'position');
    screensize=get(0,'screensize');
    width=windowsize(1)*screensize(3);
    height=windowsize(2)*screensize(4);
    x0=pos(1)+pos(3)/2-width/2;
    y0=pos(2)+pos(4)/2-height/2;
    hhelp=figure('position',[x0,y0,width,height],'name',titlestring,...
        'toolbar','none','menubar','none');
    uicontrol(hhelp,'style','text','units','normalized','position',[.1 .1 .8 .8],...
        'fontsize',fontsize,'horizontalalignment','left','string',helpmsg);
end