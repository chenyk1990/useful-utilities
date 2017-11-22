function plotimage3D(seis,t,xline,iline,dname,cmap,gridx,gridy)
% PLOTIMAGE3D: provides interactive ability to browse and study a 3D seismic volume
%
% plotimage3D(seis,t,xline,iline,dname,cmap,gridx,gridy)
%
% Seismic data stored in a 3D matrix are presented in a way that facilites easy viewing. Data can be
% viewed as either inline (iline) or crossline (xline) vertical sections or as time slices. The 3D
% matrix must be organized with time as the first dimension, xline as the second dimension, and
% iline as the third dimension. The views presented are always 2D panels from the 3D volume and are
% presented as images (i.e. not wiggle traces). Controls are provided to adjust cliping and to step
% through the volume sequentially. Additionally, when more than one plotimage3D window is active,
% then they can be grouped (linked) and browsed simultaneously. Changing the view in any window of a
% group causes all members of the group to show the same view.  Other functionality includes: easy
% copying of any view to the Windows clipboard (for pasting into PowerPoint), adjusting the display
% clip level, changing the colormap, saving views of interest for easy return, easy cursor location
% in data coordinates, flipping inline and crossline axes direction. Extended analysis functions are
% available through right-clicking on a displayed image. For cross sections, the available options
% are f-k spectrum, time-variant f sprectrum, and f-x analysis sections.
%
% seis ... 3D seismic matrix. Should be a regular grid with dimension 1
%       (row) as time, dimension 2 (column) as xline, and dimension 3 as
%       iline.
% t ... time coordinate vector for seis. Must has the same length as
%       size(seis,1)
% xline ... crossline coordinate or x coordinate. Must have length equal to
%       size(seis,2)
% iline ... inline coordinate or y coordinate. Must have length equal to
%       size(seis,3)
% dataname ... string giving a name for the dataset
% ************ default = '' ***********
% cmap ... initial colormap. This is a text string and should be one of
%     colormaps={'seisclrs','parula','jet','hsv','copper','autumn','bone'...
%         'gray','cool','winter','spring','alpine','summer','hot'};
% ************ default = 'seisclrs' ***********
% gridx ... physical distance between crosslines
% ************ default abs(xline(2)-xline(1)) **********
% gridy ... physical distance between inlines
% ************ default abs(iline(2)-iline(1)) **********
% NOTE: gridx and gridy are only important if you intend to examine 2D spectra of inline, crossline,
% or timeslice views or do wavenumber filtering (all actions accessed by right-clicking in image
% views). The default values will not give physically correct wavenumbers for these processes.
% 
% NOTE2: Missing traces should be filled with zeros (not nan's). The presence of nan's in the data
% volume will cause the program to fail. 
%
% NOTE3: This function is designed to work with the specific size window that it creates by default.
% If you resize this window significantly, then the axes labels (coordinates) may become incorrect.
% There is currently no workaround for this.
%
% NOTE4: When reading into memory using readsegy, a 3D survey will be stored as a 2D matrix in trace
% sequential mode (i.e. one trace after another). You must move these traces into a 3D matrix in
% order to use plotimage3D. Unless the 3D survey has a perfectly square spatial aperture, this will
% generally involve padding with zero traces. You can use 'make3dvol' for this purpose. Here is an
% axample:
%
% Let path and fname be strings whose concatenation points to the input dataset
% [seis,dt,segfmt,texthdrfmt,byteorder,texthdr,binhdr,exthdr,traceheaders] =readsegy([path fname]);
% 
% t=dt*(0:size(seis,1)-1)'; %time coordinate vector
% ilineall=double(traceheaders.InlineNum); %get inline numbers from trace headers
% xlineall=double(traceheaders.XlineNum); %get xline numbers from trace headers
% cdpx=double(traceheaders.CdpX); %get cdpx from trace headers
% cdpy=double(traceheaders.CdpY); %get cdpy from trace headers
%
% [seis3D,xline,iline,xcdp,ycdp,kxline]=make3Dvol(seis,xlineall,ilineall,cdpx,cdpy);
%
% plotimage3D(seis3D,t,xline,iline,'My 3D dataset');
%
% This example expects the inline and crossline numbers and the inline and crossline cdp values to
% be stored in the SEGY trace headers in the standard places. It will fail if this is not the case
% with your data.
%
% G.F. Margrave, Devon Energy, 2016-2017
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

%USER DATA information
%
% control ... tagstring ... userdata
% hbasemap ... 'basemap' ... {seis, t, xline, iline, dname, [smean sdev smin smax], gridx, gridy}
% hinline ... 'inlinebox' ... ykol (yellow color user to indicate which mode is active)
% htslice ... 'tslicebox' ... dt (time sample rate)
% hb1 ... 'inlinemode' ... [hb2 hb3 hinline hxline htslice] (handles of various controls)
% hampapply ... 'ampapply' ... [hmax hmin hmaxlbl hminlbl] (handles of manual clipping controls)
% hlocate ... 'locate' ... a flag used by locate mechanism
% hinlinebutton ... 'inlinemode' ... 5 handles
% htslicebutton ... 'tslicemode' ... not used
% hxlinebutton ... 'xlinemode' ... not used
% hpreviousbutton ... 'previous' ... info for context menu and other uses cell array
%               {seismic_section, xcoord, ycoord, mode, 3rdcoord, dname, dx, dy} See updateview
% hnextbutton ... 'next' ... not used
% hincrementbutton ... 'increment' ... notused
% htmin ... 'tmin' ... time increment for tmin and tmax 
% htmax ... 'tmax' ... not used
% 
% To find any of these handles, just use h=findobj(gcf,'tag','tagstring')
% 
% 
global PLOTIMAGE3DFIGS PLOTIMAGE3DDATASIZE PLOTIMAGE3DDIFFDIAL PLOTIMAGE3DINFODIAL PLOTIMAGE3DMASTER PLOTIMAGE3DTHISFIG;
% Description of globals
% PLOTIMAGE3DFIGS ... Ordinary array of Figure handles of the grouped figures. The first figure is
%           considered the master figure.
% PLOTIMAGE3DDATASIZE ... ordinary array containing size numbers for each grouped dataset. The first
%           3 are the dataset sizes in t,x, and y. The second 3 are shifts for each of the 3
%           dimensions needed to align with the master figure.
% PLOTIMAGE3DDIFFDIAL ... handle to the difference dialog.
% PLOTIMAGE3DINFODIAL ... handle to the information dialog
% PLOTIMAGE3DMASTER ... handle of the master figure in the group
% PLOTIMAGE3DTHISFIG ... set by SANE when it controls the plotimage3D figure and we are doing a grouping operation
%

%installing an analysis tool
%these are accessed via right clicking on a displayed image plot
% 1) Go to the internal function updateview (in this file)
% 2) There are 3 possible views (modes): 'inline', 'xline', 'tslice'. For each view you can set the
%       right-click action via a context menu entry. There is a switch statment with a case for
%       each type of view.
% 3) Decide which of the views (can be all) that you wish your tool to be active on. Then, inside
%       the case statement for that view, locate the line hcm=uicontextmenu; Below that line are a
%       series of uimenu statements, one for each tool that appears on a right click. Deuplicate one
%       of the existing lines and modify it for your tool. The last entry in the uimenu call is the
%       callback that executes when the menu is selected.
% 4) Choose a unique name for your function callback and then write the correspondiong internal
%       function. For example, the 2D spectra callback is @show2dspectrum and corresponding to this
%       there is the internal function function show2dspectrum(~,~). The data that you want to
%       operate on are found in the userdata for hprevious. See the first few lines of show2dspectrum
%       to see how to get this data.
if(ischar(seis))
    action=seis;
else
    action='init';
end

if(strcmp(action,'init'))
    
    PLOTIMAGE3DDIFFDIAL=[];
    
    if(nargin<5)
        dname='';
    end
    if(nargin<6)
        cmap='seisclrs';
    end
    if(nargin<7)
        gridx=abs(xline(2)-xline(1));
    end
    if(nargin<8)
        gridy=abs(iline(2)-iline(1));
    end
    if(~isa(seis,'single'))
        seiss=single(seis);
        clear seis;
    else
        seiss=seis;
        clear seis;
    end
    [nt,nx,ny]=size(seiss);
    
    if(length(t)~=nt)
        error('t is the wrong size')
    end
    if(length(xline)~=nx)
        error('xline is wrong size');
    end
    if(length(iline)~=ny)
        error('iline is wrong size');
    end
    xline=xline(:)';%row vector
    iline=iline(:);%column vector
    xx=xline(ones(size(iline)),:);
    yy=iline(:,ones(size(xline)));
    
    figure
    bigfig
    pos=get(gcf,'position');
    %get plotting statistics
    %first determine where the zero traces are
    set(gcf,'name',['plotimage3D ... ' dname]);
    map=squeeze(sum(abs(seiss(end-100,:,:)),1))';
    %ideadtr=find(map==0);
    ilivetr=find(map~=0);
    %map(ixlive,iylive)=1;
    xlive=xx(ilivetr);
    ylive=yy(ilivetr);
    %determine max min mean and std by examining 10 inlines and 10 xlines
    n1=min([10 length(xline)]);
    n2=min([10 length(iline)]);
    idel=round(length(xline)/n1);
    ix1=10:idel:length(xline);
    idel=round(length(iline)/n2);
    iy1=10:idel:length(iline);
    amax1=zeros(1,length(ix1));
    amin1=amax1;
    std1=amax1;
    amean1=amax1;
    ns1=amax1;
    for k=1:length(ix1)
       tmp=squeeze(seiss(:,ix1(k),:));
       ilive=find(tmp~=0);
       ns1(k)=length(ilive);
       amax1(k)=max(tmp(ilive));
       amin1(k)=min(tmp(ilive));
       std1(k)=std(tmp(ilive));
       amean1(k)=mean(tmp(ilive));
    end
    amax2=zeros(1,length(iy1));
    amin2=amax2;
    std2=amax2;
    amean2=amax2;
    ns2=amax2;
    for k=1:length(iy1)
       tmp=squeeze(seiss(:,:,iy1(k)));
       ilive=find(tmp~=0);
       ns2(k)=length(ilive);
       amax2(k)=max(tmp(ilive));
       amin2(k)=min(tmp(ilive));
       std2(k)=std(tmp(ilive));
       amean2(k)=mean(tmp(ilive));
    end
    smean=sum([amean1.*ns1 amean2.*ns2])/sum([ns1 ns2]);
    sdev=sqrt(sum([std1.^2.*ns1 std2.^2.*ns2])/sum([ns1 ns2]));
    smax=max([amax1 amax2]);
    smin=min([amin1 amin2]);
%     ilive=find(seiss);%live samples (really just the nonzero samples)
%     smean=mean(seiss(ilive));
%     sdev=std(seiss(ilive));
%     smax=max(seiss(ilive));
%     smin=min(seiss(ilive));
    
    %make the basemap axes
    xnot=.05;ynot=.1;
    width=.15;ht=.2;
    ynow=1-ynot-ht;
    xnow=xnot;
    hbmap=axes('position',[xnow,ynow,width,ht],'tag','basemap');
    set(gcf,'nextplot','add');
    hh=plot(xlive(:),ylive(:),'r.','markersize',.1);flipy;
    set(gcf,'nextplot','new');
    set(hbmap,'tag','basemap');
    set(hh,'color',[.5 .5 .5]);
    set(hbmap,'yaxislocation','right');
    xlabel('crossline');ylabel('inline');
    title('basemap')
    xmin=min(xline);
    ymin=min(iline);
    xmax=max(xline);
    ymax=max(iline);
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    set(hbmap,'userdata',{seiss, t, xline, iline, dname, [smean sdev smin smax], gridx, gridy})
    
    %mode buttons
    sep=.05;
    ht=.03;
    width=.05;
    ynow=ynow-sep-ht;
    ykol=[1 1 .5];
    fs=10;
    if(pos(3)<1500)
        fs=7;
    end
    if(pos(3)<900)
        fs=6;
    end

    hb1=uicontrol(gcf,'style','pushbutton','string','Inline','tag','inlinemode',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''inlinemode'')','backgroundcolor',ykol,'fontsize',fs);
    sep=.01;
    xnow=xnow+width+sep;
    hb2=uicontrol(gcf,'style','pushbutton','string','Xline','tag','xlinemode',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''xlinemode'')','fontsize',fs);
    xnow=xnow+width+sep;
    hb3=uicontrol(gcf,'style','pushbutton','string','Tslice','tag','tslicemode',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''tslicemode'')','fontsize',fs);
    
    
    xnow=xnot;
    ynow=ynow-sep-ht;
    inot=round(length(iline)/2);%first inline to display
    hinline=uicontrol(gcf,'style','edit','string',int2str(iline(inot)),'tag','inlinebox',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''inline'')','backgroundcolor',ykol,'fontsize',fs,...
        'userdata',ykol,'tooltipstring',...
        ['inline number to view (min=' int2str(min(iline)) ', max=' int2str(max(iline)) ')']);
    sep=.01;
    xnow=xnow+width+sep;
    hxline=uicontrol(gcf,'style','edit','string','all','tag','xlinebox',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''xline'')','fontsize',fs,'tooltipstring',...
        ['xline number to view (min=' int2str(min(xline)) ', max=' int2str(max(xline)) ')']);
    xnow=xnow+width+sep;
    htslice=uicontrol(gcf,'style','edit','string','all','tag','tslicebox',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''tslice'')','fontsize',fs,'userdata',t(2)-t(1),...
        'tooltipstring',...
        ['timeslice to view (min=' num2str(min(t)) ', max=' num2str(max(t)) ')']);
    
    set(hb1,'userdata',[hb2 hb3 hinline hxline htslice]);
    
    %prev, next and increment
    xnow=xnot;
    ynow=ynow-sep-ht;
    uicontrol(gcf,'style','pushbutton','string','previous','tag','previous',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''previous'')','fontsize',fs,'tooltipstring',...
        'increment the view to the previous inline/xline/tslice');
    sep=.01;
    xnow=xnow+width+sep;
    uicontrol(gcf,'style','pushbutton','string','next','tag','next',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''next'')','fontsize',fs,'tooltipstring',...
        'increment the view to the next inline/xline/tslice');
    xnow=xnow+width+sep;
    w2=width/2-.5*sep;
    nudge=.25*ht;
    uicontrol(gcf,'style','text','string','incr:','units','normalized',...
        'position',[xnow,ynow-nudge,w2,ht],'fontsize',fs,...
        'tooltipstring','increment for prev and next');
    xnow=xnow+w2;
    uicontrol(gcf,'style','edit','string','1','tag','increment',...
        'units','normalized','position',[xnow ynow w2+sep ht],'callback',...
        'plotimage3D(''increment'')','fontsize',fs,...
        'tooltipstring','specify increment in samples');
    
    %tmin and tmax
    xnow=xnot;
    ynow=ynow-sep-ht;
    tmin=min(t);
    tmax=max(t);
    tinc=.1;
    tminvalues=tmin:tinc:tmax-tinc;
    tminlabels=num2strcell(tminvalues,-1);
    width2=.5*width;
    uicontrol(gcf,'style','text','string','Tmin:','tag','tminlabel',...
        'units','normalized','position',[xnow ynow-nudge width2 ht],'fontsize',fs,...
        'tooltipstring','Choose minimum display time.');
    sep=.01;
    xnow=xnow+width2;
    uicontrol(gcf,'style','popupmenu','string',tminlabels,'tag','tmin',...
        'units','normalized','position',[xnow ynow width2 ht],'callback',...
        'plotimage3D(''tmin_tmax'')','fontsize',fs,'value',1,'userdata',tinc);
    xnow=xnow+sep+width2;
    tmin=min(t);
    tmax=max(t);
    tmaxvalues=tmin+tinc:tinc:tmax;
    tmaxlabels=num2strcell(tmaxvalues,-1);
    uicontrol(gcf,'style','text','string','Tmax:','tag','tmaxlabel',...
        'units','normalized','position',[xnow ynow-nudge width2 ht],'fontsize',fs,...
        'tooltipstring','Choose maximum display time.');
    sep=.01;
    xnow=xnow+width2;
    uicontrol(gcf,'style','popupmenu','string',tmaxlabels,'tag','tmax',...
        'units','normalized','position',[xnow ynow width2 ht],'callback',...
        'plotimage3D(''tmin_tmax'')','fontsize',fs,'value',length(tmaxlabels));
    
    %clip control
    xnow=xnot;
    ynow=ynow-3*sep-ht;
    uicontrol(gcf,'style','text','string','Clip level:','tag','cliplabel',...
        'units','normalized','position',[xnow ynow-nudge width ht],'fontsize',fs,...
        'tooltipstring','Choose clipping level, smaller means more clipping.');
    sep=.01;
    xnow=xnow+width;
    cliplevels={'manual','30','20','15','10','8','7','6','5','4','3','2','1','.5','.25','.1','.05'};
    iclip=10;%starting clip level
    uicontrol(gcf,'style','popupmenu','string',cliplevels,'tag','cliplevel',...
        'units','normalized','position',[xnow ynow .75*width ht],'callback',...
        'plotimage3D(''clip'')','fontsize',fs,'value',iclip,...
        'tooltipstring','Value is the number of standard deviations from the mean at which clipping occurs.');
    
    %manual amplitude controls
    xnow=xnow+.5*width+2*sep;
    ynow=ynow+.75*ht;
    vis='off';
    hmaxlbl=uicontrol(gcf,'style','text','string','max:','tag','maxamplbl',...
        'units','normalized','position',[xnow,ynow-nudge,.5*width',ht],...
        'fontsize',fs,'visible',vis,'tooltipstring',...
        'Enter the maximum amplitude to be displayed without clipping.');
    xnow=xnow+.5*width;
    hmax=uicontrol(gcf,'style','edit','string',num2str(smean+str2double(cliplevels{iclip})*sdev),'tag','maxamp',...
        'units','normalized','position',[xnow ynow width ht],'fontsize',...
        fs,'visible',vis,...
        'tooltipstring','Walue shown is the current clipping maximum.');
    xnow=xnow-.5*width;
    ynow=ynow-ht;
    hminlbl=uicontrol(gcf,'style','text','string','min:','tag','minamplbl',...
        'units','normalized','position',[xnow,ynow-nudge,.5*width',ht],...
        'fontsize',fs,'visible',vis,'tooltipstring',...
        'Enter the minimum amplitude to be displayed without clipping.');
    xnow=xnow+.5*width;
    hmin=uicontrol(gcf,'style','edit','string',num2str(smean-str2double(cliplevels{iclip})*sdev),'tag','minamp',...
        'units','normalized','position',[xnow ynow width ht],'fontsize',...
        fs,'visible',vis,...
        'tooltipstring','Value shown is the current clipping minimum.');
    ynow=ynow-ht;
    uicontrol(gcf,'style','pushbutton','string','apply','tag','ampapply',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''manualclipping'')','fontsize',fs,'visible',vis,...
        'tooltipstring','Push to apply manual clipping.',...
        'userdata',[hmax hmin hmaxlbl hminlbl]);
    
    %colormap control
    xnow=xnot;
    ynow=ynow-1.25*ht;
    uicontrol(gcf,'style','text','string','Colormap:','tag','colomaplabel',...
        'units','normalized','position',[xnow ynow-.5*nudge width ht],'fontsize',fs);
    sep=.01;
    xnow=xnow+width;
    if(exist('parula','file')==2)
        colormaps={'seisclrs','redblue','redblue2','redblue3','blueblack','bluebrown','greenblack',...
            'greenblue','jet','parula','copper','bone','gray','winter'};
        icolor=2;%starting colormap
    else
        colormaps={'seisclrs','redblue','redblue2','redblue3','blueblack','bluebrown','greenblack',...
            'greenblue','jet','copper','bone','gray','winter'};
        icolor=2;%starting colormap
    end
    for k=1:length(colormaps)
        if(strcmp(colormaps{k},cmap))
            icolor=k;
        end
    end
    uicontrol(gcf,'style','popupmenu','string',colormaps,'tag','colormap',...
        'units','normalized','position',[xnow ynow width ht],'callback',...
        'plotimage3D(''colormap'')','fontsize',fs,'value',icolor);
    
    %grid lines
    xnow=xnow+width+.5*sep;
    uicontrol(gcf,'style','text','string','Grid lines boldness:','tag','grid',...
        'units','normalized','position',[xnow,ynow-.5*nudge,1.5*width,ht],'fontsize',fs,...
        'tooltipstring','turn coordinate grid on or off');
    %ynow=ynow-ht;
    xnow2=xnow+1.5*width;
    gridoptions={'off','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
    uicontrol(gcf,'style','popupmenu','string',gridoptions,'tooltipstring',...
        'larger number means darker grid lines','units','normalized',...
        'position',[xnow2,ynow,.5*width,ht],'callback','plotimage3D(''grid'')',...
        'value',1,'fontsize',fs,'tag','gridoptions');
    %grid color
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Grid color:','tag','grid',...
        'units','normalized','position',[xnow,ynow-.5*nudge,.75*width,ht],'fontsize',fs,...
        'tooltipstring','Choose grid color');
    gridcolors={'black','white','red','dark red','blue','green','dark green','cyan','magenta','yellow','dark yellow'};
    xnow2=xnow2-.75*width;
    uicontrol(gcf,'style','popupmenu','string',gridcolors,'tooltipstring',...
        'Choose the grid color','units','normalized',...
        'position',[xnow2,ynow,width,ht],'callback','plotimage3D(''grid'')',...
        'value',1,'fontsize',fs,'tag','gridcolors','userdata',{'k','w','r',[.5 0 0],'b','g',[0 .5 0],'c','m','y',[.8 .8 0]});
    
    %brightness
    xnow=xnot;
    %ynow=ynow-ht;
    uicontrol(gcf,'style','text','string','Brightness:','tag','brightnesslabel',...
        'units','normalized','position',[xnow ynow-.5*nudge width ht],'fontsize',fs);
    sep=.01;
    xnow=xnow+width;
    brightnesses={'0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0','-0.1','-0.2','-0.3','-0.4','-0.5','-0.6','-0.7','-0.8'};
    uicontrol(gcf,'style','popupmenu','string',brightnesses,'tag','brighten',...
        'units','normalized','position',[xnow ynow .6*width ht],'callback',...
        'plotimage3D(''colormap'')','fontsize',fs,'value',9);
    %flip colormap
    xnow=xnot;
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Flip Colormap:','tag','flipcolormaplabel',...
        'units','normalized','position',[xnow ynow 1.5*width ht],'fontsize',fs);
    xnow=xnow+1.5*width;
    uicontrol(gcf,'style','radiobutton','tag','flipcolormap','value',0,...
        'units','normalized','position',[xnow ynow+.9*nudge .3*width ht],'callback','plotimage3D(''colormap'')');
    %copy to clipboard
    xnow=xnot;
    ynow=ynow-ht;
    if(isunix)
        msg='To TIFF file without controls';
        msg2='Save current view to a TIFF file';
    else
        msg='To Clipboard without controls';
        msg2='Copy current view to the WINDOWS clipboard without control panel';
    end
    uicontrol(gcf,'style','pushbutton','string',msg,'tag','clipboardalt',...
        'units','normalized','position',[xnow,ynow,2*width ht],'callback',...
        'plotimage3D(''clipboardalt'')','fontsize',fs,'tooltipstring',msg2);
    if(isunix)
        msg='To TIFF file with controls';
        msg2='Save current view to a TIFF file';
    else
        msg='To Clipboard with controls';
        msg2='Copy current view to the WINDOWS clipboard with control panel';
    end
    uicontrol(gcf,'style','pushbutton','string',msg,'tag','clipboard',...
        'units','normalized','position',[xnow+2*width,ynow,2*width ht],'callback',...
        'plotimage3D(''clipboard'')','fontsize',fs,'tooltipstring',msg2);
    
    
    %add to group
    ynow=ynow-ht-.5*sep;
    uicontrol(gcf,'style','pushbutton','string','Add to Group','units',...
        'normalized','position',[xnow,ynow,1.3*width,ht],'callback',...
        'plotimage3D(''group'')','fontsize',fs,'tooltipstring',...
        'Include in group of linked plotimage3D figures','tag','group');
    %remove from group
    xnow=xnow+1.3*width;
    uicontrol(gcf,'style','pushbutton','string','Remove from Group','units',...
        'normalized','position',[xnow,ynow,1.3*width,ht],'callback',...
        'plotimage3D(''ungroup'')','fontsize',fs,'tooltipstring',...
        'Remove from group of linked plotimage3D figures','tag','ungroup');
    
    %group info
    xnow=xnot;
    ynow=ynow-ht-.5*sep;
    uicontrol(gcf,'style','pushbutton','string','Group info','units',...
        'normalized','position',[xnow,ynow,1.3*width,ht],'callback',...
        'plotimage3D(''groupinfo'');','fontsize',fs,'tooltipstring',...
        'Show group datasets and time shifts','tag','groupinfo');
    %clear group
    xnow=xnow+1.3*width;
    uicontrol(gcf,'style','pushbutton','string','Clear Group','units',...
        'normalized','position',[xnow,ynow,1.3*width,ht],'callback',...
        'plotimage3D(''cleargroup'')','fontsize',fs,'tooltipstring',...
        'Clear the group of linked figures to start a new group',...
        'tag','cleargroup');
    %group zoom
    xnow=xnow+1.3*width;
    uicontrol(gcf,'style','pushbutton','string','Group zoom','units',...
        'normalized','position',[xnow,ynow,1.3*width,ht],'callback',...
        'plotimage3D(''groupzoom'')','fontsize',fs,'tooltipstring',...
        'Zoom one member of a group and then push this to zoom all members',...
        'tag','cleargroup');
    
    %save views
    xnow=xnot;
    ynow=ynow-ht-.5*sep;
    uicontrol(gcf,'style','pushbutton','string','Save view','units',...
        'normalized','position',[xnow,ynow,width,ht],'callback',...
        'plotimage3D(''saveview'')','fontsize',fs,'tooltipstring',...
        'Save this view for easy return','tag','saveview');
    %forget views
    xnow=xnow+width;
    uicontrol(gcf,'style','pushbutton','string','Forget view','units',...
        'normalized','position',[xnow,ynow,width,ht],'callback',...
        'plotimage3D(''forgetview'')','fontsize',fs,'tooltipstring',...
        'Remove from list of saved views','tag','forgetview');
    %restore views
    xnow=xnow+width;
    uicontrol(gcf,'style','popupmenu','string',{'Saved views'},'units',...
        'normalized','position',[xnow,ynow,1.5*width,ht],'callback',...
        'plotimage3D(''restoreview'')','fontsize',fs,'tooltipstring',...
        'Remove from list of saved views','tag','savedviews');
    
    %cursor locate button
    xnow=xnot;
    ynow=ynow-ht-.5*sep;
    uicontrol(gcf,'style','pushbutton','string','cursor locate on','units',...
        'normalized','position',[xnow ynow 1.2*width ht],'tag','locate',...
        'fontsize',fs,'tooltipstring','Turn on location information at cursor',...
        'callback','plotimage3D(''locate'')','userdata',[]);
    %flipx button
    xnow=xnow+1.2*width;
    uicontrol(gcf,'style','pushbutton','string','flip xline','units',...
        'normalized','position',[xnow ynow width ht],'tag','locate',...
        'fontsize',fs,'tooltipstring','Reverse the x axis',...
        'callback','plotimage3D(''flipx'')','userdata',1,'tag','flipx');
    %userdata, 1 for normal -1 for reversed. This refers to the seismic
    %axis which uses the image convention for axis direction. For images, 
    %xdir normal increases to the right while ydir normal increases down.
    %For normal plots, xdir is the same but ydir increases up.
    %flipy button
    xnow=xnow+width;
    uicontrol(gcf,'style','pushbutton','string','flip inline','units',...
        'normalized','position',[xnow ynow width ht],'tag','locate',...
        'fontsize',fs,'tooltipstring','Reverse the y axis',...
        'callback','plotimage3D(''flipy'')','userdata',1,'tag','flipy');
    ynow=ynow-ht-.5*sep;
    xnow=xnot;
    uicontrol(gcf,'style','text','string','Associated figure windows:','units','normalized',...
        'position',[xnow,ynow-.5*nudge,2*width,ht],'fontsize',fs);
    xnow=xnow+2*width;
    uicontrol(gcf,'style','popupmenu','string',{'None'},'tag','windows',...
        'units','normalized','position',[xnow ynow 2*width ht],'fontsize',fs,...
        'callback','plotimage3D(''show_window'')');
    %userdata of this popup is a vector of windows handles
    
    set(gcf,'closerequestfcn','plotimage3D(''close'');');
    
    %make the seismic axes and show the first inline
    xnow=xnot+3*sep+5*width;
    width=.6;
    ht=.8;
    ynow=ynot;
    %hseismic=axes('position',[xnow,ynow,width,ht],'tag','seismic');
    axes('position',[xnow,ynow,width,ht],'tag','seismic');
    updateview;
    
elseif(strcmp(action,'inline'))
    %hfigs=getfigs;
    hthisfig=gcf;
    hinlinemode=findobj(hthisfig,'tag','inlinemode');
    udat=get(hinlinemode,'userdata');
    hxlinemode=udat(1);
    htslicemode=udat(2);
    hinline=udat(3);
    hxline=udat(4);
    htslice=udat(5);
    kol_off=.94*ones(1,3);
    kol_on=get(hinline,'userdata');
    tmp=get(hinline,'string');
    iline_chosen=str2double(tmp);
    if(isnan(iline_chosen))
        msgbox('You must enter an integer number to chose an inline',...
            'Ooops!');
        return;
    end
    hbmap=findobj(hthisfig,'tag','basemap');
    udat=get(hbmap,'userdata');
    iline=udat{4};
    imin=min(iline);imax=max(iline);
    if(iline_chosen<imin || iline_chosen>imax)
        msgbox(['Invalid inline number, must be between ' int2str(imin) ...
            ' and ' int2str(imax)],'Ooops!');
        return;
    end
    set(hinlinemode,'backgroundcolor',kol_on);
    set(hinline,'backgroundcolor',kol_on);
    set(hxlinemode,'backgroundcolor',kol_off);
    set(hxline,'backgroundcolor',kol_off,'string','all');
    set(htslicemode,'backgroundcolor',kol_off);
    set(htslice,'backgroundcolor',kol_off,'string','all');
    updateview;
    
elseif(strcmp(action,'xline'))
    %hfigs=getfigs;
    hthisfig=gcf;
    hinlinemode=findobj(hthisfig,'tag','inlinemode');
    udat=get(hinlinemode,'userdata');
    hxlinemode=udat(1);
    htslicemode=udat(2);
    hinline=udat(3);
    hxline=udat(4);
    htslice=udat(5);
    kol_off=.94*ones(1,3);
    kol_on=get(hinline,'userdata');
    tmp=get(hxline,'string');
    xline_chosen=str2double(tmp);
    if(isnan(xline_chosen))
        msgbox('You must enter an integer number to chose an xline',...
            'Ooops!');
        return;
    end
    hbmap=findobj(hthisfig,'tag','basemap');
    udat=get(hbmap,'userdata');
    xline=udat{3};
    xlmin=min(xline);xlmax=max(xline);
    if(xline_chosen<xlmin || xline_chosen>xlmax)
        msgbox(['Invalid xline number, must be between ' int2str(xlmin) ...
            ' and ' int2str(xlmax)],'Ooops!');
        return;
    end
    set(hinlinemode,'backgroundcolor',kol_off);
    set(hinline,'backgroundcolor',kol_off,'string','all');
    set(hxlinemode,'backgroundcolor',kol_on);
    set(hxline,'backgroundcolor',kol_on);
    set(htslicemode,'backgroundcolor',kol_off);
    set(htslice,'backgroundcolor',kol_off,'string','all');
    updateview;
    
elseif(strcmp(action,'tslice'))
    %hfigs=getfigs;
    hthisfig=gcf;
    hinlinemode=findobj(hthisfig,'tag','inlinemode');
    udat=get(hinlinemode,'userdata');
    hxlinemode=udat(1);
    htslicemode=udat(2);
    hinline=udat(3);
    hxline=udat(4);
    htslice=udat(5);
    kol_off=.94*ones(1,3);
    kol_on=get(hinline,'userdata');
    tmp=get(htslice,'string');
    tslice_chosen=str2double(tmp);
    if(isnan(tslice_chosen))
        msgbox('You must enter a valid time to chose a tslice',...
            'Ooops!');
        return;
    end
    hbmap=findobj(hthisfig,'tag','basemap');
    udat=get(hbmap,'userdata');
    t=udat{2};
    tmin=min(t);tmax=max(t);
    %check for milliseconds
    if(tslice_chosen>tmax)
        if(round(tslice_chosen)==tslice_chosen)
            tslice_chosen=tslice_chosen/1000;
            set(htslice,'string',num2str(tslice_chosen));
        end
    end
    if(tslice_chosen<tmin || tslice_chosen>tmax)
        msgbox(['Invalid time, must be between ' num2str(tmin) ...
            ' and ' num2str(tmax)],'Ooops!');
        return;
    end
    %make sure the time requested is one that we have
    it=near(t,tslice_chosen);
    tslice_chosen=t(it(1));
    set(htslice,'string',num2str(tslice_chosen));
    set(hinlinemode,'backgroundcolor',kol_off);
    set(hinline,'backgroundcolor',kol_off,'string','all');
    set(hxlinemode,'backgroundcolor',kol_off);
    set(hxline,'backgroundcolor',kol_off,'string','all');
    set(htslicemode,'backgroundcolor',kol_on);
    set(htslice,'backgroundcolor',kol_on);
    updateview;
    
elseif(strcmp(action,'previous'))
    hincr=findobj(gcf,'tag','increment');
    tmp=get(hincr,'string');
    inc=str2double(tmp);
    mode=determinemode;
    hseis=findobj(gcf,'tag','seismic');
    switch mode
        case 'inline'
            hinline=findobj(gcf,'tag','inlinebox');
            tmp=get(hinline,'string');
            inlinenow=str2double(tmp);
            inlinenext=inlinenow-inc;
            set(hinline,'string',int2str(inlinenext));
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            updateview;
            set(hseis,'xlim',xl,'ylim',yl);
            
        case 'xline'
            hxline=findobj(gcf,'tag','xlinebox');
            tmp=get(hxline,'string');
            xlinenow=str2double(tmp);
            xlinenext=xlinenow-inc;
            set(hxline,'string',int2str(xlinenext));
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            updateview;
            set(hseis,'xlim',xl,'ylim',yl);
            
        case 'tslice'
            htslice=findobj(gcf,'tag','tslicebox');
            dt=get(htslice,'userdata');
            inc=inc*dt;
            tmp=get(htslice,'string');
            tslicenow=str2double(tmp);
            tslicenext=tslicenow-inc;
            
            set(htslice,'string',num2str(tslicenext));
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            updateview;
            set(hseis,'xlim',xl,'ylim',yl);
    end
    
elseif(strcmp(action,'next'))
    hincr=findobj(gcf,'tag','increment');
    tmp=get(hincr,'string');
    inc=str2double(tmp);
    mode=determinemode;
    hseis=findobj(gcf,'tag','seismic');
    switch mode
        case 'inline'
            hinline=findobj(gcf,'tag','inlinebox');
            tmp=get(hinline,'string');
            inlinenow=str2double(tmp);
            inlinenext=inlinenow+inc;
            set(hinline,'string',int2str(inlinenext));
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            updateview;
            set(hseis,'xlim',xl,'ylim',yl);
        case 'xline'
            hxline=findobj(gcf,'tag','xlinebox');
            tmp=get(hxline,'string');
            xlinenow=str2double(tmp);
            xlinenext=xlinenow+inc;
            set(hxline,'string',int2str(xlinenext));
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            updateview;
            set(hseis,'xlim',xl,'ylim',yl);
            
        case 'tslice'
            htslice=findobj(gcf,'tag','tslicebox');
            dt=get(htslice,'userdata');
            inc=inc*dt;
            tmp=get(htslice,'string');
            tslicenow=str2double(tmp);
            tslicenext=tslicenow+inc;
            set(htslice,'string',num2str(tslicenext));
            xl=get(hseis,'xlim');
            yl=get(hseis,'ylim');
            updateview;
            set(hseis,'xlim',xl,'ylim',yl);
            
    end
    
elseif(strcmp(action,'clip'))
    hfigs=getfigs;
    hthisfig=gcf;
    clipnow=getclip;
    hampapply=findobj(hthisfig,'tag','ampapply');
    hampcontrols=get(hampapply,'userdata');
    if(length(clipnow)>1)
        set([hampapply hampcontrols],'visible','on');
        return;
    else
        set([hampapply hampcontrols],'visible','off');
    end
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    %get the dat amp values
    hbmap=findobj(hthisfig,'tag','basemap');
    udat=get(hbmap,'userdata');
    amp=udat{6};
    %update the values in the invisible max and min edit fields
    set(hampcontrols(1),'string',num2str(amp(1)+clipnow*amp(2)));%maximum
    set(hampcontrols(2),'string',num2str(amp(1)-clipnow*amp(2)));%minimum
    %get the seismic axes and update its clim property
    hseismic=findobj(hthisfig,'tag','seismic');
    clim=[amp(1)-clipnow*amp(2), amp(1)+clipnow*amp(2)];
    set(hseismic,'clim',clim);
    %process the other figs
    hclip=findobj(hthisfig,'tag','cliplevel');
    iclip=get(hclip,'value');
    for k=1:length(hotherfigs)
        hbmap=findobj(hotherfigs(k),'tag','basemap');
        udat=get(hbmap,'userdata');
        amp=udat{6};
        hampapply=findobj(hotherfigs(k),'tag','ampapply');
        hampcontrols=get(hampapply,'userdata');
        set(hampcontrols(1),'string',num2str(amp(1)+clipnow*amp(2)));%maximum
        set(hampcontrols(2),'string',num2str(amp(1)-clipnow*amp(2)));%minimum
        hcliplevel=findobj(hotherfigs(k),'tag','cliplevel');
        set(hcliplevel,'value',iclip);
        hseismic=findobj(hotherfigs(k),'tag','seismic');
        clim=[amp(1)-clipnow*amp(2), amp(1)+clipnow*amp(2)];
        set(hseismic,'clim',clim);
    end
    
elseif(strcmp(action,'manualclipping'))
    %hfigs=getfigs;
    hthisfig=gcf;
%     ind=hfigs~=hthisfig;
%     hotherfigs=hfigs(ind);
    clim=getclip;
    if(length(clim)~=2)
        error('logic failure');
    end
    hseismic=findobj(hthisfig,'tag','seismic');
    set(hseismic,'clim',clim);
    
elseif(strcmp(action,'colormap'))
    hfigs=getfigs;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    cmapnow=getcolormap;
    colormap(cmapnow);
    hcolormap=findobj(hthisfig,'tag','colormap');
    icolor=get(hcolormap,'value');
    cmap=get(hthisfig,'colormap');
    hbrighten=findobj(hthisfig,'tag','brighten');
    ibright=get(hbrighten,'value');
    hflip=findobj(hthisfig,'tag','flipcolormap');
    flip=get(hflip,'value');
    for k=1:length(hotherfigs)
        hcolormap=findobj(hotherfigs(k),'tag','colormap');
        set(hcolormap,'value',icolor);
        set(hotherfigs(k),'colormap',cmap);
        hbrighten=findobj(hotherfigs(k),'tag','brighten');
        set(hbrighten,'value',ibright)
        hflip=findobj(hotherfigs(k),'tag','flipcolormap');
        set(hflip,'value',flip);
    end
    
elseif(strcmp(action,'clipboard')||strcmp(action,'clipboardalt'))
    if(strcmp(action,'clipboardalt'))
        hidecontrols;
    end
    if(isunix)
        print -dtiff
        adon='tiff file in current directory';
    else
        fh=gcf;
        fh.Renderer='opengl';
        hgexport(fh,'-clipboard');
        adon='Windows clipboard';
    end
    if(strcmp(action,'clipboardalt'))
        restorecontrols;
    end
    msg=['Figure has been sent to ' adon];
    msgbox(msg,'Good news!');
    
elseif(strcmp(action,'grid'))
    hthisfig=gcf;
    hfigs=getfigs;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hgridopt=findobj(hthisfig,'tag','gridoptions');
    hgridcolor=findobj(hthisfig,'tag','gridcolors');
    hseismic=findobj(hthisfig,'tag','seismic');
    axes(hseismic);
    opt=get(hgridopt,'value');
    gridoptions=get(hgridopt,'string');
    iklr=get(hgridcolor,'value');
    klrs=get(hgridcolor,'userdata');
    klr=klrs{iklr};
    if(opt==1)
        grid off;
    else
        alpha=str2double(gridoptions{opt});
        grid on;
        set(hseismic,'gridalpha',alpha,'gridcolor',klr); 
    end
    %handle other figs
    for k=1:length(hotherfigs)
        hgridopt=findobj(hotherfigs(k),'tag','gridoptions');
        hgridcolor=findobj(hotherfigs(k),'tag','gridcolors');
        hseismic=findobj(hotherfigs(k),'tag','seismic');
        axes(hseismic); %#ok<LAXES>
        set(hgridopt,'value',opt);
        set(hgridcolor','value',iklr);
        if(opt==1)
            grid off;
        else
            alpha=str2double(gridoptions{opt});
            grid on;
            set(hseismic,'gridalpha',alpha,'gridcolor',klr);
        end
    end
    
elseif(strcmp(action,'group')||strcmp(action,'groupex'))
    %'groupex' is initiated by SANE to cause grouping.  
    if(~isempty(PLOTIMAGE3DDIFFDIAL))
        delete(PLOTIMAGE3DDIFFDIAL);
        PLOTIMAGE3DDIFFDIAL=[];
    end
    PLOTIMAGE3DMASTER=[];
    if(strcmp(action,'group'))
        hthisfig=gcf;
    else
        hthisfig=PLOTIMAGE3DTHISFIG;
    end
    hbmap=findobj(hthisfig,'tag','basemap');
    udat=get(hbmap,'userdata');
    if(~isempty(PLOTIMAGE3DFIGS))
        %check data size for compatibility
%         tmp=PLOTIMAGE3DDATASIZE;
%         if(length(udat{2})~=tmp(1)||length(udat{3})~=tmp(2)||length(udat{4})~=tmp(3))
%             msgbox('Cannot include this figure in group because data size is not compatible',...
%                 'Ooops!');
%             return;
%         end
        PLOTIMAGE3DDATASIZE=[PLOTIMAGE3DDATASIZE; [length(udat{2}) length(udat{3}) length(udat{4}) 0 0 0]];
    else
        %the first three numbers are nt nx and ny (dataset sizes) 
        %and the second three are dt dx and dy (units of seconds and line numbers) that give shifts
        %to match the master figure (first figure). Usually the master Figure will have zero shifts
        %unless the original master was removed from group. The basic rule relating time slice on
        %dataset k with that on dataset j is tk+dtk = tj+dtj (all values in seconds) and therefore
        %tj=tk+dtk-dtj.  If datasets do not have the same sample rate then this might not work out
        %to an exact sample.
        PLOTIMAGE3DDATASIZE=[length(udat{2}) length(udat{3}) length(udat{4}) 0 0 0];
    end
    ind=find(hthisfig==PLOTIMAGE3DFIGS,1);%see if we have already included it
    if(isempty(ind))
        if(~isempty(PLOTIMAGE3DFIGS))
            %equalize the displays
            hgroupfig=PLOTIMAGE3DFIGS(end);%previously the last member of the group
            
            h=findobj(hgroupfig,'tag','tmin');
            tmins=str2double(get(h,'string'));
            ival=get(h,'value');
            tmin=tmins(ival);
            h=findobj(hthisfig,'tag','tmin');
            tmins=str2double(get(h,'string'));
            ival=near(tmins,tmin);
            set(h,'value',ival);
            
            h=findobj(hgroupfig,'tag','tmax');
            tmaxs=str2double(get(h,'string'));
            ival=get(h,'value');
            tmax=tmaxs(ival);
            h=findobj(hthisfig,'tag','tmax');
            tmaxs=str2double(get(h,'string'));
            ival=near(tmaxs,tmax);
            set(h,'value',ival);
            
            h=findobj(hgroupfig,'tag','cliplevel');
            val=get(h,'value');
            h=findobj(hthisfig,'tag','cliplevel');
            set(h,'value',val);
            
            h=findobj(hgroupfig,'tag','colormap');
            val=get(h,'value');
            h=findobj(hthisfig,'tag','colormap');
            set(h,'value',val);
            
            h=findobj(hgroupfig,'tag','gridoptions');
            val=get(h,'value');
            h=findobj(hthisfig,'tag','gridoptions');
            set(h,'value',val);
            
            h=findobj(hgroupfig,'tag','brighten');
            val=get(h,'value');
            h=findobj(hthisfig,'tag','brighten');
            set(h,'value',val);
            
            h=findobj(hgroupfig,'tag','gridcolors');
            val=get(h,'value');
            h=findobj(hthisfig,'tag','gridcolors');
            set(h,'value',val);
            
            h=findobj(hgroupfig,'tag','flipcolormap');
            val=get(h,'value');
            h=findobj(hthisfig,'tag','flipcolormap');
            set(h,'value',val);
        end
        
        PLOTIMAGE3DFIGS=[PLOTIMAGE3DFIGS hthisfig];
        
        %send message to SANE
        if(issane)%issane is an internal function
            sane('pi3d:group',sanedata);%sanedata is an internal function
        end
        
        thisname=get(hthisfig,'name');
        msgbox([thisname ' added to group of linked figures'],'Done');
    end
    
elseif(strcmp(action,'ungroup')||strcmp(action,'ungroupex'))
    %'ungroupex' is initiated by SANE to cause ungrouping. 
    if(~isempty(PLOTIMAGE3DDIFFDIAL))
        delete(PLOTIMAGE3DDIFFDIAL);
        PLOTIMAGE3DDIFFDIAL=[];
    end
    if(strcmp(action,'ungroup'))
        hthisfig=gcf;
    else
        hthisfig=PLOTIMAGE3DTHISFIG;
    end
    ind=find(hthisfig==PLOTIMAGE3DFIGS,1);
    if(~isempty(ind))
        PLOTIMAGE3DFIGS(ind)=[];
        PLOTIMAGE3DDATASIZE(ind,:)=[];

        %send message to SANE
        if(issane)%issane is an internal function
            sane('pi3d:group',sanedata);%sanedata is an internal function
        end
        
        thisname=get(hthisfig,'name');
        msgbox([thisname ' removed from group of linked figures'],'Done');
    end
%     if(isempty(PLOTIMAGE3DFIGS))
%         PLOTIMAGE3DDATASIZE=[];
%     end
elseif(strcmp(action,'cleargroup'))
    if(~isempty(PLOTIMAGE3DDIFFDIAL))
        delete(PLOTIMAGE3DDIFFDIAL);
        PLOTIMAGE3DDIFFDIAL=[];
    end
    if(~isgraphics(PLOTIMAGE3DINFODIAL))
        delete(PLOTIMAGE3DINFODIAL)
    end
    PLOTIMAGE3DFIGS=[];
    PLOTIMAGE3DDATASIZE=[];
    PLOTIMAGE3DINFODIAL=[];
    
    %send message to SANE
    if(issane)%issane is an internal function
        sane('pi3d:group',sanedata);%sanedata is an internal function
    end
    
    msgbox(' Existing group cleared, start a new one! ','Done');
    
elseif(strcmp(action,'groupinfo'))
    if(isempty(PLOTIMAGE3DFIGS))
       msgbox('You have not yet created a group (or the group is empty) so there is nothing to show');
       return
    end
    hfigs=PLOTIMAGE3DFIGS;%get the grouped plotseis3D figures
    if(~isgraphics(PLOTIMAGE3DINFODIAL))
        delete(PLOTIMAGE3DINFODIAL)
    end
    hthisfig=gcf;
    tabledata=cell(length(hfigs),2);
    %get the time shifts
    figdata=PLOTIMAGE3DDATASIZE;
    %load up the tabledata
    for k=1:length(hfigs)
        hprevious=findobj(hfigs(k),'tag','previous');
        udat=get(hprevious,'userdata');
        name=strrep(udat{6},'plotimage3D ... ','');
        tabledata{k,1}=name;%names
        tabledata{k,2}=figdata(k,4);%time shifts
        %time shift rule tj+dtj=tk+dtk
    end

    %create a dialog
    pos=get(hthisfig,'position');
    width=pos(3)*.5;
    ht=pos(4)*.25;
    xnow=pos(1)+.5*(pos(3)-width);
    ynow=pos(2)+.5*(pos(4)-ht);
    hdial=figure('position',[xnow,ynow,width,ht]);
    pos=get(hdial,'position');
    width=.8;ht=.6;
    xnow=.1;ynow=.2;
    uitable(gcf,'units','normalized','position',[xnow,ynow,width,ht],'data',tabledata,...
        'columneditable',[false true],'columnname',{'Dataset name','Time shift (seconds)'},...
        'columnwidth',{pos(3)*width*.7,pos(3)*width*.2},'tag','table');
    
    %dismiss button
    xnow=.05;
    ynow=.1;
    width=.2;
    ht=.1;
    uicontrol(hdial,'style','pushbutton','string','Apply&Dismiss','tag','dismiss','units','normalized',...
        'position',[xnow,ynow,width,ht],'callback','plotimage3D(''dismissinfo'');',...
        'backgroundcolor','y','tooltipstring','Click to apply changes and dismiss this dialog',...
        'userdata',hthisfig);

    PLOTIMAGE3DINFODIAL=hdial;
    set(hdial,'closerequestfcn','plotimage3D(''dismissinfo'')','name','Plotimage3D Group Info');
elseif(strcmp(action,'groupzoom'))
    hfigs=getfigs;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hseis=findobj(hthisfig,'tag','seismic');
    xl=get(hseis,'xlim');
    yl=get(hseis,'ylim');
    set(hseis,'xlim',xl,'ylim',yl)
    for k=1:length(hotherfigs)
        hseis=findobj(hotherfigs(k),'tag','seismic');
        set(hseis,'xlim',xl,'ylim',yl);
    end
elseif(strcmp(action,'dismissinfo'))
    hdial=gcf;
    hfigs=getfigs;
    currentdata=PLOTIMAGE3DDATASIZE;
    htable=findobj(hdial,'tag','table');
    tabledata=get(htable,'data');
    newshifts=zeros(size(tabledata,1),1);
    for k=1:length(newshifts)
       newshifts(k)=tabledata{k,2};
    end
    currentshifts=currentdata(:,4);
    for k=1:length(newshifts)
        if(newshifts(k)~=currentshifts(k))
            %update view
            figure(hfigs(k))
            view=currentview;
            hbmap=findobj(gcf,'tag','basemap');
            udat=get(hbmap,'userdata');
            t=udat{2};
            dt=abs(t(2)-t(1));
            newshifts(k)=dt*round(newshifts(k)/dt);%round to nearest sample
            PLOTIMAGE3DDATASIZE(k,4)=newshifts(k);
%             if(strcmp(view{4},'tslice'))
%                 htslice=findobj(gcf,'tag','tslicebox');
%                 txt=get(htslice,'string');
%                 tnow=str2double(txt);
%                 
%                 tnow=tnow-currentshifts(k)+newshifts(k);
%                 set(htslice,'string',num2str(tnow));
%             end
            setview(view);
        end
    end
    delete(hdial);
    PLOTIMAGE3DINFODIAL=[];
    
elseif(strcmp(action,'saveview'))
    hfigs=getfigs;
    view=currentview;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hsavedviews=findobj(hthisfig,'tag','savedviews');
    viewlist=get(hsavedviews,'string');
    nviews=length(viewlist);
    switch view{4}
        case 'inline'
            viewlist{nviews+1}=[view{4} ': ' view{1}];
        case 'xline'
            viewlist{nviews+1}=[view{4} ': ' view{2}];
        case 'tslice'
            viewlist{nviews+1}=[view{4} ': ' view{3}];
    end
    %viewlist{nviews+1}=[view{4} ' inline: ' view{1} ' xline: ' view{2} ' tslice: ' view{3}];
    set(hsavedviews,'string',viewlist,'value',nviews+1);
    udat=get(hsavedviews,'userdata');
    nviews=length(udat);
    udat{nviews+1}=view;
    set(hsavedviews,'userdata',udat);
    %process the view menus of the other figs
    for k=1:length(hotherfigs)
        hsavedviews=findobj(hotherfigs(k),'tag','savedviews');
        viewlist=get(hsavedviews,'string');
        nviews=length(viewlist);
        switch view{4}
            case 'inline'
                viewlist{nviews+1}=[view{4} ': ' view{1}];
            case 'xline'
                viewlist{nviews+1}=[view{4} ': ' view{2}];
            case 'tslice'
                viewlist{nviews+1}=[view{4} ': ' view{3}];
        end
        %viewlist{nviews+1}=[view{4} ' inline: ' view{1} ' xline: ' view{2} ' tslice: ' view{3}];
        set(hsavedviews,'string',viewlist,'value',nviews+1);
        udat=get(hsavedviews,'userdata');
        nviews=length(udat);
        udat{nviews+1}=view;
        set(hsavedviews,'userdata',udat);
    end
    
elseif(strcmp(action,'forgetview'))
    hfigs=getfigs;
    view=currentview;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hsavedviews=findobj(hthisfig,'tag','savedviews');
    viewlist=get(hsavedviews,'userdata');
    %search for the current view among the saved views
    nviews=length(viewlist);
    iview=[];
    for k=1:nviews
        if(strcmp(view{1},viewlist{k}{1}))
            if(strcmp(view{2},viewlist{k}{2}))
                if(strcmp(view{3},viewlist{k}{3}))
                    if(strcmp(view{4},viewlist{k}{4}))
                        iview=k;
                    end
                end
            end
        end
    end
    if(isempty(iview))
%         msgbox('Current view is not in the saved list','Ooops!');
        return
    end
    viewlist(iview)=[];
    set(hsavedviews,'userdata',viewlist);
    viewnames=get(hsavedviews,'string');
    viewnames(iview+1)=[];
    set(hsavedviews,'string',viewnames,'value',1);
    %process other figs
    for k=1:length(hotherfigs)
        hsavedviews=findobj(hotherfigs(k),'tag','savedviews');
        viewlist=get(hsavedviews,'userdata');
        %search for the current view among the saved views
        nviews=length(viewlist);
        iview=[];
        for kk=1:nviews
            if(strcmp(view{1},viewlist{kk}{1}))
                if(strcmp(view{2},viewlist{kk}{2}))
                    if(strcmp(view{3},viewlist{kk}{3}))
                        if(strcmp(view{4},viewlist{kk}{4}))
                            iview=kk;
                        end
                    end
                end
            end
        end
        if(isempty(iview))
            %         msgbox('Current view is not in the saved list','Ooops!');
            return
        end
        viewlist(iview)=[];
        set(hsavedviews,'userdata',viewlist);
        viewnames=get(hsavedviews,'string');
        viewnames(iview+1)=[];
        set(hsavedviews,'string',viewnames,'value',1);
    end
    
elseif(strcmp(action,'restoreview'))
    hsavedviews=findobj(gcf,'tag','savedviews');
    views=get(hsavedviews,'userdata');
    desiredview=get(hsavedviews,'value')-1;
    if(desiredview>0)
        setview(views{desiredview});
    end
    
elseif(strcmp(action,'locate')||strcmp(action,'stoplocate'))
    hfigs=getfigs;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hlocate=findobj(gcf,'tag','locate');
    if(strcmp(action,'locate'))
        flag=1;
    else
        flag=0;
    end
    if(flag==1) %turn it on
        set(hlocate,'string','cursor locate off','callback','plotimage3D(''stoplocate'')');
        set(gcf,'windowbuttondownfcn','plotimage3D(''postlocation'')');
    else
        set(hlocate,'string','cursor locate on','callback','plotimage3D(''locate'')');
        set(gcf,'windowbuttondownfcn','');
        clearlocations;
    end
    
    %process other figs
    for k=1:length(hotherfigs)
        hlocate=findobj(hotherfigs(k),'tag','locate');
        if(flag==1) %turn it on
            set(hlocate,'string','cursor locate off','callback','plotimage3D(''stoplocate'')');
            set(hotherfigs(k),'windowbuttondownfcn','plotimage3D(''postlocation'')');
        else
            set(hlocate,'string','cursor locate on','callback','plotimage3D(''locate'')');
            set(hotherfigs(k),'windowbuttondownfcn','');
            clearlocations;
        end
    end
elseif(strcmp(action,'postlocation'))
    hlocate=findobj(gcf,'tag','locate');
    hfigs=getfigs;
    hthisfig=gcf;
    hseismic=findobj(hthisfig,'tag','seismic');
    existinglocations=get(hlocate,'userdata');
    currentpoint=get(hseismic,'currentpoint');
    xl=get(hseismic,'xlim');
    yl=get(hseismic,'ylim');
    if(currentpoint(1,1)<xl(1) || currentpoint(1,1)>xl(2))
        return;
    end
    if(currentpoint(1,2)<yl(1) || currentpoint(1,2)>yl(2))
        return;
    end
    % existinglocations are stored in a cell array, one entry per point. Each point is represented
    % by a two element vector of handles. The first is the handle to the line with the point marker
    % and the second is the text. check existing locations. If we have a match, then this is a
    % deletion
    if(isempty(currentpoint))
        return;
    end
    npts=length(existinglocations);
    badpoints=zeros(size(existinglocations));
    %small=10^12*eps;
    hclicked=gco;
    for k=1:npts
        thispoint=existinglocations{k};
        if(isempty(thispoint))
            badpoints(k)=1;
        else
            if(~isgraphics(thispoint(1)))
                badpoints(k)=1;
            else
%                 xpt=get(thispoint(1),'xdata');
%                 ypt=get(thispoint(1),'ydata');
%                 test=abs(currentpoint(1,1)-xpt)+abs(currentpoint(1,2)-ypt);
%                 if(test<small)
                if(hclicked==thispoint(1))
                    delete(thispoint);
                    return;
                end
            end
        end
    end
    ind=find(badpoints==1);
    if(~isempty(ind))
        existinglocations(ind)=[];
    end
    mode=determinemode;
    %use the gridcolor as the color
    kol=get(hseismic,'gridcolor');
    fs=9;mksize=6;
    if(strcmp(mode,'tslice'))
        newpoint(2)=text(currentpoint(1,1),currentpoint(1,2),...
            ['(' int2str(currentpoint(1,1)) ',' int2str(currentpoint(1,2)) ')'],...
            'fontsize',fs','color',kol);
    else
        newpoint(2)=text(currentpoint(1,1),currentpoint(1,2),...
            ['(' int2str(currentpoint(1,1)) ',' num2str(round(1000*currentpoint(1,2))/1000) ')'],...
            'fontsize',fs','color',kol);
    end
    newpoint(1)=line(currentpoint(1,1),currentpoint(1,2),'linestyle','none',...
        'marker','*','markersize',mksize,'color',kol);
    existinglocations{npts+1}=newpoint;
    set(hlocate,'userdata',existinglocations);
    
    %post the point in other figs
    ind= hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    for k=1:length(hotherfigs)
        figure(hotherfigs(k));
        hseismic=findobj(hotherfigs(k),'tag','seismic');
        hlocate=findobj(hotherfigs(k),'tag','locate');
        existinglocations=get(hlocate,'userdata');
        set(hotherfigs(k),'currentaxes',hseismic);
        if(strcmp(mode,'tslice'))
            newpoint(2)=text(currentpoint(1,1),currentpoint(1,2),...
                ['(' int2str(currentpoint(1,1)) ',' int2str(currentpoint(1,2)) ')'],...
                'fontsize',fs','color',kol);
        else
            newpoint(2)=text(currentpoint(1,1),currentpoint(1,2),...
                ['(' int2str(currentpoint(1,1)) ',' num2str(round(1000*currentpoint(1,2))/1000) ')'],...
                'fontsize',fs','color',kol);
        end
        newpoint(1)=line(currentpoint(1,1),currentpoint(1,2),'linestyle','none',...
            'marker','*','markersize',mksize,'color',kol);
        existinglocations{npts+1}=newpoint;
        set(hlocate,'userdata',existinglocations);
    end
    
    
elseif(strcmp(action,'flipx'))
    hfigs=getfigs;
    hthisfig=gcf;
    hflipx=findobj(hthisfig,'tag','flipx');
    dirflag=get(hflipx,'userdata');
    if(dirflag==1)
        set(hflipx,'userdata',-1);
    else
        set(hflipx,'userdata',1);
    end
    setaxesdir;
    
    %now handle other figs in group
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    for k=1:length(hotherfigs)
        figure(hotherfigs(k));
        hflipx=findobj(hotherfigs(k),'tag','flipx');
        if(dirflag==1)
            set(hflipx,'userdata',-1);
        else
            set(hflipx,'userdata',1);
        end
        setaxesdir;
    end
elseif(strcmp(action,'flipy'))
    hfigs=getfigs;
    hthisfig=gcf;
    hflipy=findobj(hthisfig,'tag','flipy');
    dirflag=get(hflipy,'userdata');
    if(dirflag==1)
        set(hflipy,'userdata',-1);
    else
        set(hflipy,'userdata',1);
    end
    setaxesdir;
    
    %now handle other figs in group
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    for k=1:length(hotherfigs)
        figure(hotherfigs(k));
        hflipy=findobj(hotherfigs(k),'tag','flipy');
        if(dirflag==1)
            set(hflipy,'userdata',-1);
        else
            set(hflipy,'userdata',1);
        end
        setaxesdir;
    end
    
elseif(strcmp(action,'increment'))
    hfigs=getfigs;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hinc=findobj(hthisfig,'tag','increment');
    val=get(hinc,'string');
    for k=1:length(hotherfigs)
       hinc=findobj(hotherfigs(k),'tag','increment');
       set(hinc,'string',val);
    end
    
elseif(strcmp(action,'tmin_tmax'))
    hfigs=getfigs;
    hthisfig=gcf;
    ind=hfigs~=hthisfig;
    hotherfigs=hfigs(ind);
    hobj=gcbo;
    htmin=findobj(gcf,'tag','tmin');
    htmax=findobj(gcf,'tag','tmax');
    tinc=get(htmin,'userdata');
    if(hobj==htmin)
        %tmin is changed
        tmin=get_tmin;%current tmin
        tmax=get_tmax;%current tmax
        tmaxlbl=get(htmax,'string');%current tmax labels
        Tmax=str2double(tmaxlbl{end});%largest possible tmax
        tmaxs=tmin+tinc:tinc:Tmax;%new tmax vector
        imax=near(tmaxs,tmax);
        tmaxlbl=num2strcell(tmaxs,-1);%new tmax labels
        set(htmax,'string',tmaxlbl,'value',imax);%reset tmax
        %process other figs
        for k=1:length(hotherfigs)
           htmin=findobj(hotherfigs(k),'tag','tmin');
           tmins=str2double(get(htmin,'string'));
           imin=near(tmins,tmin);
           set(htmin,'value',imin)
           htmax=findobj(hotherfigs(k),'tag','tmax');
           set(htmax,'string',tmaxlbl,'value',imax);%reset tmax
        end
    else
        %tmax is changed
        tmin=get_tmin;%current tmin
        tmax=get_tmax;%current tmax
        tminlbl=get(htmin,'string');%current tminx labels
        Tmin=str2double(tminlbl{1});%smallest possible tmin
        tmins=Tmin:tinc:tmax-tinc;%new tmin vector
        imin=near(tmins,tmin);
        tminlbl=num2strcell(tmins,-1);%new tmax labels
        set(htmin,'string',tminlbl,'value',imin);%reset tmax
        %process other figs
        for k=1:length(hotherfigs)
           htmax=findobj(hotherfigs(k),'tag','tmax');
           tmaxs=str2double(get(htmax,'string'));
           imax=near(tmaxs,tmax);
           set(htmax,'value',imax)
           htmin=findobj(hotherfigs(k),'tag','tmin');
           set(htmin,'string',tminlbl,'value',imin);%reset tmin
        end
    end
    updateview;
elseif(strcmp(action,'dismissdifference'))
    if(~isempty(PLOTIMAGE3DDIFFDIAL))
        PLOTIMAGE3DDIFFDIAL=[];
    end
    delete(gcf);
    
elseif(strcmp(action,'close'))
    button=questdlg('Are you sure you want to close this PLOTIMAGE3D dataset?','Just to be sure...');
    if(strcmp(button,'Yes'))
        hthisfig=gcbf;
        hwin=findobj(hthisfig,'tag','windows');
        hfigs=get(hwin,'userdata');
        for k=1:length(hfigs)
            if(isgraphics(hfigs(k)))
                delete(hfigs(k))
            end
        end
        delete(hthisfig);
    end
elseif(strcmp(action,'closesane'))
    hsane=gcbf;
    if(strcmp(get(hsane,'tag'),'sane'))
        hbut=gcbo;
        hpan=get(get(hbut,'parent'),'parent');
        idat=get(hpan,'userdata');%the is the number of the dataset whose window is being closed
        hfile=findobj(hsane,'tag','file');
        %hmsg=findobj(hsane,'tag','message');
        proj=get(hfile,'userdata');
        hthisfig=proj.pifigures{idat};
    else
        hthisfig=hsane;
    end
    hwin=findobj(hthisfig,'tag','windows');
    hfigs=get(hwin,'userdata');
    for k=1:length(hfigs)
        if(isgraphics(hfigs(k)))
            delete(hfigs(k))
        end
    end
    delete(hthisfig);
elseif(strcmp(action,'datanamechange'))
    %called by SANE to change the data name, arg2 (t) will be the Figure handle, and arg3 (x) will be
    %the new dataname
    hfig=t;
    dname=xline;
    hbmap=findobj(hfig,'tag','basemap');
    udat=get(hbmap,'userdata');
    oldname=udat{5};
    udat{5}=dname;
    set(hfig,'name',['plotimage3D ... ' dname]);
    hseismic=findobj(hfig,'tag','seismic');
    ht=get(hseismic,'title');
    ht.String=strrep(ht.String,oldname,dname);
    set(hbmap,'userdata',udat);
    hprevious=findobj(hfig,'tag','previous');
    udat=get(hprevious,'userdata');
    udat{6}=dname;
    set(hprevious,'userdata',udat)
elseif(strcmp(action,'closewindow'))
    hwin=gcf;
    hthisfig=get(hwin,'userdata');
    if(~isgraphics(hthisfig))
        delete(hwin);
        return
    end
    hwinlist=findobj(hthisfig,'tag','windows');
    winfigs=get(hwinlist,'userdata');
    winnames=get(hwinlist,'string');
    ival=get(hwinlist,'value');
    nwins=length(winfigs);
    if(nwins==0)
        return;%happens if string is 'None'
    end
    iwin=find(hwin==winfigs, 1);
    if(isempty(iwin))
        return;%not sure why this might happen
    end
    winnames(iwin)=[];
    winfigs(iwin)=[];
    if(ival==iwin)
        ival=min([length(winfigs) ival+1]);
    elseif(ival>iwin)
        ival=ival-1; 
    end
    if(isempty(winfigs))
        winnames{1}='None';
        ival=1;
    end
    set(hwinlist,'string',winnames,'value',ival,'userdata',winfigs);
    delete(hwin)
elseif(strcmp(action,'show_window'))
    hwin=gcbo;
    hfigs=get(hwin,'userdata');
    iwin=get(hwin,'value');
    if(isgraphics(hfigs(iwin)))
        figure(hfigs(iwin));
    else
        nwin=length(hfigs);
        fignames=get(hwin,'string');
        fignames(iwin)=[];
        figlist(iwin)=[];
        nwin=nwin-1;
        if(iwin>nwin)
            iwin=nwin;
        end
        if(isempty(figlist))
            fignames{1}='None';
            iwin=1;
        end
        set(hwin,'string',fignames,'value',iwin);
    end
end
end

function mode=determinemode
hthisfig=gcf;
hinlinemode=findobj(hthisfig,'tag','inlinemode');
udat=get(hinlinemode,'userdata');
% hxlinemode=udat(1);
% htslicemode=udat(2);
hinline=udat(3);
hxline=udat(4);
htslice=udat(5);
%determine mode
ykol=get(hinline,'userdata');%this is a flag to indicate mode
kol=get(hinline,'backgroundcolor');
if(kol==ykol)
    mode='inline';
end
kol=get(hxline,'backgroundcolor');
if(kol==ykol)
    mode='xline';
end
kol=get(htslice,'backgroundcolor');
if(kol==ykol)
    mode='tslice';
end
end

function updateview
global PLOTIMAGE3DMASTER PLOTIMAGE3DDIFFDIAL

if(~isempty(PLOTIMAGE3DDIFFDIAL))
    delete(PLOTIMAGE3DDIFFDIAL);
    PLOTIMAGE3DDIFFDIAL=[];
end

hfigs=getfigs;%get the grouped plotseis3D figures
hthisfig=gcf;
ind=hfigs~=hthisfig;
hotherfigs=hfigs(ind);%other figures in this group
%first deal with the main figure
hinlinemode=findobj(hthisfig,'tag','inlinemode');
udat=get(hinlinemode,'userdata');
%
hinline=udat(3);
hxline=udat(4);
htslice=udat(5);
%determine mode
mode=determinemode;
%get the data
hbmap=findobj(hthisfig,'tag','basemap');
udat=get(hbmap,'userdata');
seiss=udat{1};
t=udat{2};
dt=abs(t(2)-t(1));
xline=udat{3};
iline=udat{4};
dname=udat{5};
dname=strrep(dname,'plotimage3D ... ','');
amp=udat{6};
hseismic=findobj(gcf,'tag','seismic');
%get limits of current view to impose again if just updating
xl=get(hseismic,'xlim');
yl=get(hseismic,'ylim');
tmin=get_tmin;
tmax=get_tmax;
gridx=udat{7};
gridy=udat{8};
timeshift=gettimeshift;
switch mode
    case 'inline'
        t=t+timeshift;
        tmp=get(hinline,'string');
        iline_next=str2double(tmp);
        tmp=get(hxline,'string');
        if(strcmp(tmp,'all'))
            ixline=1:length(xline);
        else
            ind=strfind(tmp,':');
            ix1=str2double(tmp(1:ind-1));
            ix2=str2double(tmp(ind+1:end));
            ixline=ix1:ix2;
        end
        tmp=get(htslice,'string');
        if(strcmp(tmp,'all'))
            %it=1:length(t);
            it=near(t,tmin,tmax);
        else
            ind=strfind(tmp,':');
            it1=str2double(tmp(1:ind-1));
            it2=str2double(tmp(ind+1:end));
            it=it1:it2;
        end
        axes(hseismic)
        %find the image if it exists
        hi=findobj(gca,'type','image');
        inot=near(iline,iline_next);
        if(~isempty(hi))
            tag=get(hi,'tag');
            if(~strcmp(tag,'inline'))
                delete(hi);
                hi=[];
            end
        end
        if(isempty(hi))
            %make a new image
            %xx=1:length(xline);
            clipnow=getclip;
            if(length(clipnow)==1)
                clim=[amp(1)-clipnow*amp(2) amp(1)+clipnow*amp(2)];
            else
                clim=clipnow;
            end
            clearlocations;
            set(gcf,'nextplot','add');
            hi=imagesc(xline,t(it),squeeze(seiss(it,ixline,inot(1))),clim);
            set(gcf,'nextplot','new');
            %create a context menu
            hcm=uicontextmenu;
            uimenu(hcm,'label','Spectrum (2D)','callback',@show2dspectrum);
            uimenu(hcm,'label','Time-variant spectra','callback',@showtvspectrum);
            uimenu(hcm,'label','Time-variant relative phase','callback',@showtvphase);
            uimenu(hcm,'label','f-x phase','callback',@showfxphase);
            uimenu(hcm,'label','f-x amp','callback',@showfxamp);
            uimenu(hcm,'label','Amplitude histogram','callback',@amphist)
            uimenu(hcm,'label','SVD separation','callback',@showsvdsep);
            uimenu(hcm,'label','Difference plot','callback',@difference);
            uimenu(hcm,'label','Bandpass filter','callback',@filter);
            uimenu(hcm,'label','Spiking decon','callback',@deconvolution);
            set(hi,'uicontextmenu',hcm);
            %save the contextmenudata
            hprevious=findobj(gcf,'tag','previous');
            dname=get(gcf,'name');
            dname=strrep(dname,'plotimage3D ... ','');
                
            set(hprevious,'userdata',{squeeze(seiss(it,ixline,inot(1))),xline,t(it),'inline',iline(inot(1)),dname,gridx,dt});
            
            colorbar;
            set(hi,'tag','inline');
            %flipx;
            cmapnow=getcolormap;
            colormap(cmapnow);
            xlabel('xline number');
            ylabel('time (s)');
            ht=title([dname ' inline ' int2str(iline(inot))]);
            ht.Interpreter='none';
%             xt=get(hseismic,'xtick');
%             set(hseismic,'xticklabel',vector2textcell(xline(xt)));
            set(hseismic,'tag','seismic');
            setaxesdir;
            bigfont(hseismic,1.5,1);
            plotimage3D('grid');
        else
            %update the existing image
            set(hi,'cdata',squeeze(seiss(it,ixline,inot(1))),'ydata',t(it),'tag','inline');
            ylim([t(it(1)) t(it(end))])
            ht=title([dname ' inline ' int2str(iline(inot))]);
            ht.Interpreter='none';
            clipnow=getclip;
            if(length(clipnow)==1)
                clim=[amp(1)-clipnow*amp(2) amp(1)+clipnow*amp(2)];
            else
                clim=clipnow;
            end
            set(hseismic,'tag','seismic','clim',clim,'xlim',xl,'ylim',yl);
            %save the contextmenudata
            hprevious=findobj(gcf,'tag','previous');
            dname=get(gcf,'name');
            set(hprevious,'userdata',{squeeze(seiss(it,ixline,inot(1))),xline,t(it),'inline',iline(inot(1)),dname,gridx,dt});
            
        end
        %update the basemap
        axes(hbmap)
        hl=findobj(hbmap,'tag','currentline');
        ht=findobj(hbmap,'tag','currenttext');
        if(~isempty(hl)); delete(hl); end
        if(~isempty(ht)); delete(ht); end
        hnow=line(xline,iline(inot)*ones(size(xline)),'color','r');
        set(hnow,'tag','currentline');
        set(hbmap,'tag','basemap');
        nmid=round(length(xline)/2);
        text(xline(nmid),iline(inot),['inline ' int2str(iline(inot))],...
            'horizontalalignment','center','tag','currenttext');
    case 'xline'
        t=t+timeshift;
        tmp=get(hxline,'string');
        xline_next=str2double(tmp);
        tmp=get(hinline,'string');
        if(strcmp(tmp,'all'))
            iiline=1:length(iline);
        else
            ind=strfind(tmp,':');
            il1=str2double(tmp(1:ind-1));
            il2=str2double(tmp(ind+1:end));
            iiline=il1:il2;
        end
        tmp=get(htslice,'string');
        if(strcmp(tmp,'all'))
            %it=1:length(t);
            it=near(t,tmin,tmax);
        else
            ind=strfind(tmp,':');
            it1=str2double(tmp(1:ind-1));
            it2=str2double(tmp(ind+1:end));
            it=it1:it2;
        end
        axes(hseismic)
        %find the image if it exists
        hi=findobj(gca,'type','image');
        inot=near(xline,xline_next);
        if(~isempty(hi))
            tag=get(hi,'tag');
            if(~strcmp(tag,'xline'))
                delete(hi);
                hi=[];
            end
        end
        if(isempty(hi))
            %make a new image
            %xx=1:length(iline);
            clipnow=getclip;
            if(length(clipnow)==1)
                clim=[amp(1)-clipnow*amp(2) amp(1)+clipnow*amp(2)];
            else
                clim=clipnow;
            end
            clearlocations;
            set(gcf,'nextplot','add');
            hi=imagesc(iline,t(it),squeeze(seiss(it,inot(1),iiline)),clim);
            set(gcf,'nextplot','new');
            %create a context menu
            hcm=uicontextmenu;
            uimenu(hcm,'label','Spectrum (2D)','callback',@show2dspectrum);
            uimenu(hcm,'label','Time-variant spectra','callback',@showtvspectrum);
            uimenu(hcm,'label','Time-variant relative phase','callback',@showtvphase);
            uimenu(hcm,'label','f-x phase','callback',@showfxphase);
            uimenu(hcm,'label','f-x amp','callback',@showfxamp);
            uimenu(hcm,'label','Amplitude histogram','callback',@amphist)
            uimenu(hcm,'label','SVD separation','callback',@showsvdsep);
            uimenu(hcm,'label','Difference plot','callback',@difference);
            uimenu(hcm,'label','Bandpass filter','callback',@filter);
            uimenu(hcm,'label','Spiking decon','callback',@deconvolution);
            set(hi,'uicontextmenu',hcm);
            %save the contextmenudata
            hprevious=findobj(gcf,'tag','previous');
            dname=get(gcf,'name');
            dname=strrep(dname,'plotimage3D ... ','');
            set(hprevious,'userdata',{squeeze(seiss(it,inot(1),iiline)),iline,t(it),'xline',xline(inot(1)),dname,gridy,dt});
            
            colorbar;
            set(hi,'tag','xline');
            %flipx;
            cmapnow=getcolormap;
            colormap(cmapnow);
            xlabel('inline number');
            ylabel('time (s)');
            ht=title([dname ' xline ' int2str(xline(inot))]);
            ht.Interpreter='none';
%             xt=get(hseismic,'xtick');
%             set(hseismic,'xticklabel',vector2textcell(iline(xt)));
            set(hseismic,'tag','seismic');
            setaxesdir;            
            bigfont(hseismic,1.5,1);
            plotimage3D('grid');
        else
            %update the existing image, either clip or time range has changed
            set(hi,'cdata',squeeze(seiss(it,inot(1),iiline)),'ydata',t(it),'tag','xline');
            ylim([t(it(1)) t(it(end))])
            ht=title([dname ' xline ' int2str(xline(inot))]);
            ht.Interpreter='none';
            clipnow=getclip;
            if(length(clipnow)==1)
                clim=[amp(1)-clipnow*amp(2) amp(1)+clipnow*amp(2)];
            else
                clim=clipnow;
            end
            set(hseismic,'tag','seismic','clim',clim,'xlim',xl,'ylim',yl);
            %save the contextmenudata
            hprevious=findobj(gcf,'tag','previous');
            dname=get(gcf,'name');
            set(hprevious,'userdata',{squeeze(seiss(it,inot(1),iiline)),iline,t(it),'xline',xline(inot(1)),dname,gridy,dt});
        end
        %update the basemap
        axes(hbmap)
        hl=findobj(hbmap,'tag','currentline');
        ht=findobj(hbmap,'tag','currenttext');
        if(~isempty(hl)); delete(hl); end
        if(~isempty(ht)); delete(ht); end
        hnow=line(xline(inot)*ones(size(iline)),iline,'color','r');
        set(hnow,'tag','currentline');
        set(hbmap,'tag','basemap');
        nmid=round(length(iline)/2);
        text(xline(inot),iline(nmid),['xline ' int2str(xline(inot))],...
            'horizontalalignment','center','tag','currenttext');
        
        
    case 'tslice'
        t=t+timeshift;
        tmp=get(htslice,'string');
        tslice_next=str2double(tmp);
        tmp=get(hxline,'string');
        if(strcmp(tmp,'all'))
            ixline=1:length(xline);
        else
            ind=strfind(tmp,':');
            ix1=str2double(tmp(1:ind-1));
            ix2=str2double(tmp(ind+1:end));
            ixline=ix1:ix2;
        end
        tmp=get(hinline,'string');
        if(strcmp(tmp,'all'))
            iiline=1:length(iline);
        else
            ind=strfind(tmp,':');
            il1=str2double(tmp(1:ind-1));
            il2=str2double(tmp(ind+1:end));
            iiline=il1:il2;
        end
        axes(hseismic)
        %find the image if it exists
        hi=findobj(gca,'type','image');
        inot=near(t,tslice_next);
        if(~isempty(hi))
            tag=get(hi,'tag');
            if(~strcmp(tag,'tslice'))
                delete(hi);
                hi=[];
            end
        end
        if(isempty(hi))
            %make a new image
            %xx=1:length(xline);
            %yy=1:length(iline);
            clipnow=getclip;
            if(length(clipnow)==1)
                clim=[amp(1)-clipnow*amp(2) amp(1)+clipnow*amp(2)];
            else
                clim=clipnow;
            end
            clearlocations;
            set(gcf,'nextplot','add');
            hi=imagesc(xline,iline,squeeze(seiss(inot(1),ixline,iiline))',clim);
            set(gcf,'nextplot','new');
            %create a context menu
            hcm=uicontextmenu;
            uimenu(hcm,'label','Spectrum (2D)','callback',@show2dspectrum);
            uimenu(hcm,'label','SVD separation','callback',@showsvdsep);
            uimenu(hcm,'label','Footprint analysis','callback',@footprint);
            uimenu(hcm,'label','Surface plot','callback',@showsurf);
            uimenu(hcm,'label','Difference plots (Requires a Group)','callback',@difference);
            uimenu(hcm,'label','Amplitude histogram','callback',@amphist)
            
            set(hi,'uicontextmenu',hcm);
            %save the contextmenudata
            hprevious=findobj(gcf,'tag','previous');
            dname=get(gcf,'name');
            dname=strrep(dname,'plotimage3D ... ','');
            set(hprevious,'userdata',{squeeze(seiss(inot(1),ixline,iiline))',xline,iline,'tslice',t(inot(1)),dname,gridx,gridy});
            
            colorbar;
            set(hi,'tag','tslice');
            %flipx;
            cmapnow=getcolormap;
            colormap(cmapnow);
            xlabel('xline number');
            ylabel('inline number');
            ht=title([dname ' timeslice ' num2str(t(inot))]);
            ht.Interpreter='none';
%             xt=get(hseismic,'xtick');
%             set(hseismic,'xticklabel',vector2textcell(xline(xt)));
%             yt=get(hseismic,'ytick');
%             set(hseismic,'yticklabel',vector2textcell(iline(yt)));
            set(hseismic,'tag','seismic');
            setaxesdir;            
            bigfont(hseismic,1.5,1);
            plotimage3D('grid');
        else
            %update the existing image
            set(hi,'cdata',squeeze(seiss(inot(1),ixline,iiline))','tag','tslice');
            ht=title([dname ' timeslice ' num2str(t(inot))]);
            ht.Interpreter='none';
            clipnow=getclip;
            if(length(clipnow)==1)
                clim=[amp(1)-clipnow*amp(2) amp(1)+clipnow*amp(2)];
            else
                clim=clipnow;
            end
            set(hseismic,'tag','seismic','clim',clim,'xlim',xl,'ylim',yl);
            %save the contextmenudata
            hprevious=findobj(gcf,'tag','previous');
            dname=get(gcf,'name');
            dname=strrep(dname,'plotimage3D ... ','');
            set(hprevious,'userdata',{squeeze(seiss(inot(1),ixline,iiline))',xline,iline,'tslice',t(inot(1)),dname,gridx,gridy});
        end
        %update the basemap
        axes(hbmap)
        hl=findobj(hbmap,'tag','currentline');
        ht=findobj(hbmap,'tag','currenttext');
        if(~isempty(hl)); delete(hl); end
        if(~isempty(ht)); delete(ht); end
        xx=[xline(1) xline(end) xline(end) xline(1) xline(1)];
        yy=[iline(1) iline(1) iline(end) iline(end) iline(1)];
        hnow=line(xx,yy,'color','r');
        set(hnow,'tag','currentline');
        set(hbmap,'tag','basemap');
        xx=(xline(1)+xline(end))*.5;
        yy=(iline(1)+iline(end))*.5;
        text(xx,yy,['tslice ' num2str(t(inot))],...
            'horizontalalignment','center','tag','currenttext');
        
end

%now deal with the otherfigs
if(~isempty(hotherfigs))
    %this stuff with PLOTIMAGE3DMASTER is so that only one figure in the
    %group updates the others. Without this, they update each other
    %endlessly. So, the figure in which the click occurs calls the other
    %but the others just updatethemselves and do not call the others.
    xdir=get(hseismic,'xdir');
    ydir=get(hseismic,'ydir');
    if(isempty(PLOTIMAGE3DMASTER))
        PLOTIMAGE3DMASTER=1;
        view=currentview;
        for k=1:length(hotherfigs)
            figure(hotherfigs(k));
            setview(view);
            hseismic=findobj(hotherfigs(k),'tag','seismic');
            set(hseismic,'xdir',xdir,'ydir',ydir);
            PLOTIMAGE3DMASTER=PLOTIMAGE3DMASTER+1;%this seems useless but it makes the editor think PLOTIMAGE3DMASTER is in use
        end
        PLOTIMAGE3DMASTER=[];
    end
%     PLOTIMAGE3DMASTER=1;
%     view=currentview;
%     for k=1:length(hotherfigs)
%         figure(hotherfigs(k));
%         setview(view);
%         hseismic=findobj(hotherfigs(k),'tag','seismic');
%         set(hseismic,'xdir',xdir,'ydir',ydir);
%         PLOTIMAGE3DMASTER=PLOTIMAGE3DMASTER+1;%this seems useless but it makes the editor think PLOTIMAGE3DMASTER is in use
%     end
%     PLOTIMAGE3DMASTER=[];
end

end

function hfigs=getfigs
%check to see if current figure is the dismiss dialog
name=get(gcf,'name');
ind=strfind(name,'Group Info');
if(~isempty(ind))
    hbutt=findobj(gcf,'tag','dismiss');
    hthisfig=get(hbutt,'userdata');
else
    hthisfig=gcf;
end
global PLOTIMAGE3DFIGS
%PLOTIMAGE3DFIGS is an ordinary array of figure handles
if(isempty(PLOTIMAGE3DFIGS))
    hfigs=hthisfig;
    return;
else
    ind=isgraphics(PLOTIMAGE3DFIGS);
    PLOTIMAGE3DFIGS=PLOTIMAGE3DFIGS(ind);
    if(isempty(PLOTIMAGE3DFIGS))
        hfigs=hthisfig;
        return;
    end
    ind=find(hthisfig==PLOTIMAGE3DFIGS, 1);
    if(isempty(ind))
        hfigs=hthisfig;
    else
        hfigs=PLOTIMAGE3DFIGS;
    end
    return;
end
end

function timeshift=gettimeshift
hfigs=getfigs;
hthisfig=gcf;
global PLOTIMAGE3DDATASIZE
if(~isempty(PLOTIMAGE3DDATASIZE))
    ind= hthisfig==hfigs;
    figdata=PLOTIMAGE3DDATASIZE;
    if(length(figdata)>3)
        timeshift=figdata(ind,4);
    else
        timeshift=0;
    end
else
    timeshift=0;
end
end

function clipnow=getclip
%hfigs=getfigs;
hthisfig=gcf;
%ind=hfigs~=hthisfig;
%hotherfigs=hfigs(ind);
hclip=findobj(hthisfig,'tag','cliplevel');
cliplevels=get(hclip,'string');
iclip=get(hclip,'value');
if(strcmp(cliplevels{iclip},'manual'))
    hampapply=findobj(hthisfig,'tag','ampapply');
    hampcontrols=get(hampapply,'userdata');
    ampmax=str2double(get(hampcontrols(1),'string'));
    ampmin=str2double(get(hampcontrols(2),'string'));
    if(isnan(ampmax))
        msgbox('you have not entered a valid number for the maximum amplitude',...
            'Ooops!');
        return;
    end
    if(isnan(ampmin))
        msgbox('you have not entered a valid number for the minimum amplitude',...
            'Ooops!');
        return;
    end
    clipnow=[ampmin ampmax];
else
    clipnow=str2double(cliplevels{iclip});
end
end

function cmapnow=getcolormap 
%hfigs=getfigs;
hthisfig=gcf;
% ind=hfigs~=hthisfig;
% hotherfigs=hfigs(ind);
hcolormap=findobj(hthisfig,'tag','colormap');
colormaps=get(hcolormap,'string'); 
icolor=get(hcolormap,'value'); 
hflip=findobj(gcf,'tag','flipcolormap');
flip=get(hflip,'value');
hbrighten=findobj(gcf,'tag','brighten');
ibright=get(hbrighten,'value');
brightnesses=get(hbrighten,'string');
brightness=str2double(brightnesses{ibright});
m=128;
switch colormaps{icolor}
    case 'seisclrs'
        cmapnow=seisclrs(m);
    case 'parula'
        cmapnow=parula(m);
    case 'jet'
        cmapnow=jet(m);
    case 'redblue'
        cmapnow=redblue(m);
    case 'redblue2'
        cmapnow=redblue2(m);
    case 'redblue3'
        cmapnow=redblue3(m);
    case 'copper'
        cmapnow=copper(m);
    case 'blueblack'
        cmapnow=blueblack(m);
    case 'greenblack'
        cmapnow=greenblack(m);
    case 'bone'
        cmapnow=bone(m);
    case 'gray'
        cmapnow=gray(m);
    case 'bluebrown'
        cmapnow=bluebrown(m);
    case 'greenblue'
        cmapnow=greenblue(m);
    case 'winter'
        cmapnow=winter(m);
end
        
        
if(flip==1)
    cmapnow=flipud(cmapnow);
end
if(brightness~=0)
    cmapnow=brighten(cmapnow,brightness);
end

end

function hidecontrols
hthisfig=gcf;
hkids=get(hthisfig,'children');
hseismic=findobj(hthisfig,'tag','seismic');
hbmap=findobj(hthisfig,'tag','basemap');
vistate=cell(size(hkids));
xnot=.1;
for k=1:length(hkids)
    if(hkids(k)==hseismic)
        seisposn=get(hseismic,'position');
        set(hseismic,'position',[xnot seisposn(2) seisposn(3)+seisposn(1)-xnot seisposn(4)])
        vistate{k}='on';
    elseif(hkids(k)==hbmap)
        vistate{k}='on';
        hkidsb=get(hbmap,'children');
        visb=cell(size(hkidsb));
        for kk=1:length(hkidsb)
            visb{kk}=get(hkidsb(kk),'visible');
            set(hkidsb(kk),'visible','off');
        end
        set(hbmap,'visible','off');
    elseif(~strcmp(get(hkids(k),'type'),'colorbar'))
        vistate{k}=get(hkids(k),'visible');
        set(hkids(k),'visible','off'); 
    end
end
hcopyalt=findobj(hthisfig,'tag','clipboardalt');
set(hcopyalt,'userdata',{vistate visb seisposn});
end

function restorecontrols
hthisfig=gcf;
hkids=get(hthisfig,'children');
hseismic=findobj(hthisfig,'tag','seismic');
hbmap=findobj(hthisfig,'tag','basemap');
hcopyalt=findobj(hthisfig,'tag','clipboardalt');
udat=get(hcopyalt,'userdata');
vistate=udat{1};
visb=udat{2};
seisposn=udat{3};
for k=1:length(hkids)
    if(hkids(k)==hseismic)
        set(hseismic,'position',seisposn);
    elseif(hkids(k)==hbmap)
        hkidsb=get(hbmap,'children');
        for kk=1:length(hkidsb)
           set(hkidsb,'visible',visb{kk}); 
        end
        set(hbmap,'visible',vistate{k});
    elseif(~strcmp(get(hkids(k),'type'),'colorbar'))
        set(hkids(k),'visible',vistate{k});
    end
end
end

function view=currentview
hthisfig=gcf;
%get the 3 handles
hinline=findobj(hthisfig,'tag','inlinebox');
hxline=findobj(hthisfig,'tag','xlinebox');
htslice=findobj(hthisfig,'tag','tslicebox');
%get the text entries defining the view
itext=get(hinline,'string');
xtext=get(hxline,'string');
ttext=get(htslice,'string');
%determine the mode
mode=determinemode;
%define the view, it consists of text strings giving inline xline and time and the mode.
%mode is one of the three strings 'inline' 'xline' or 'tslice'
view={itext xtext ttext mode};
end

function setview(view)
hthisfig=gcf;
%get the 3 handles
hinline=findobj(hthisfig,'tag','inlinebox');
hxline=findobj(hthisfig,'tag','xlinebox');
htslice=findobj(hthisfig,'tag','tslicebox');
%set the view
set(hinline,'string',view{1});
set(hxline,'string',view{2});
set(htslice,'string',view{3});
%execute the proper callback
if(strcmp(view{4},'inline'))
    plotimage3D('inline');
elseif(strcmp(view{4},'xline'))
    plotimage3D('xline');
else
    plotimage3D('tslice');
end
end

function setaxesdir
hflipx=findobj(gcf,'tag','flipx');
hflipy=findobj(gcf,'tag','flipy');
hseismic=findobj(gcf,'tag','seismic');
hbmap=findobj(gcf,'tag','basemap');
dirflag=get(hflipx,'userdata');
mode=determinemode;
if(dirflag==1)
    set(hseismic,'xdir','normal');
    set(hbmap,'xdir','normal');
else
    set(hseismic,'xdir','reverse');
    set(hbmap,'xdir','reverse');
end
dirflag=get(hflipy,'userdata');
if(dirflag==1)
    if(strcmp(mode,'tslice'))
        set(hseismic,'ydir','reverse');
        set(hbmap,'ydir','reverse');
    else
        set(hseismic,'ydir','reverse');
        set(hbmap,'ydir','reverse');
    end
else
    if(strcmp(mode,'tslice'))
        set(hseismic,'ydir','normal');
        set(hbmap,'ydir','normal');
    else
        set(hseismic,'ydir','reverse');
        set(hbmap,'ydir','normal');
    end
end
end

function clearlocations
hfigs=getfigs;
hthisfig=gcf;
ind= hfigs~=hthisfig;
hotherfigs=hfigs(ind);
hlocate=findobj(hthisfig,'tag','locate');
existinglocations=get(hlocate,'userdata');
npts=length(existinglocations);
for k=1:npts
    thispoint=existinglocations{k};
    delete(thispoint);
end
set(hlocate,'userdata',[]);

%process the other figs
for k=1:length(hotherfigs)
    hlocate=findobj(hotherfigs(k),'tag','locate');
    existinglocations=get(hlocate,'userdata');
    npts=length(existinglocations);
    for kk=1:npts
        thispoint=existinglocations{kk};
        delete(thispoint);
    end
    set(hlocate,'userdata',[]);
end
end

function show2dspectrum(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};
gridx=udat{7};
gridy=udat{8};
cmap=get(hthisfig,'colormap');
hseis=findobj(hthisfig,'tag','seismic');
gridon=get(hseis,'xgrid');
gridcolor=get(hseis,'gridcolor');
gridalpha=get(hseis,'gridalpha');
switch mode
    case 'inline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];
        dt=y(2)-y(1);
        fmax=.25/dt;
        fnyq=.5/dt;
        dx=x(2)-x(1);
        kxnyq=.5/dx;
        pos=get(gcf,'position');
        datar=seisplotfk(data,y,x,dname2,fmax,gridy,gridx);
        if(strcmp(gridon,'on'))
            set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
            set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
        end
        colormap(cmap);
        set(gcf,'position',pos);
        axes(datar{1})
        xlabel('xline number')
        htit=get(gca,'title');
        str=get(htit,'string');
        str{2}=['x-t space dx=' num2str(dx) ', dt=' num2str(dt)];
        set(htit,'string',str);
        htit=get(datar{2},'title');
        str=get(htit,'string');
        str{2}=['kx-f space, kxnyq=' num2str(kxnyq) ', fnyq=' num2str(fnyq)];
        set(htit,'string',str);
        xdir=get(hseis,'xdir');
        ydir=get(hseis,'ydir');
        set(datar{1},'xdir',xdir,'ydir',ydir)
    case 'xline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];
        dt=y(2)-y(1);
        fmax=.25/dt;
        fnyq=.5/dt;
        dy=x(2)-x(1);
        kynyq=.5/dy;
        pos=get(gcf,'position');
        datar=seisplotfk(data,y,x,dname2,fmax,gridy,gridx);
        if(strcmp(gridon,'on'))
            set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
            set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
        end
        colormap(cmap)
        set(gcf,'position',pos);
        axes(datar{1})
        xlabel('inline number')
        htit=get(gca,'title');
        str=get(htit,'string');
        str{2}=['y-t space dy=' num2str(dy) ', dt=' num2str(dt)];
        set(htit,'string',str);
        htit=get(datar{2},'title');
        str=get(htit,'string');
        str{2}=['ky-f space, kynyq=' num2str(kynyq) ', fnyq=' num2str(fnyq)];
        set(htit,'string',str);
        xdir=get(hseis,'xdir');
        ydir=get(hseis,'ydir');
        set(datar{1},'xdir',xdir,'ydir',ydir)
    case 'tslice'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];
        dy=y(2)-y(1);
        dx=x(2)-x(1);
        kynyq=1/dy;
        kxnyq=1/dx;
        pos=get(gcf,'position');
        datar=seisplotfk(data,y,x,dname2,kynyq,gridy,gridx);
        colormap(cmap)
        if(strcmp(gridon,'on'))
            set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
            set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
        end
        set(gcf,'position',pos);
        axes(datar{1})
        xlabel('crossline number')
        ylabel('inline number')
        htit=get(gca,'title');
        str=get(htit,'string');
        str{2}=['x-y space dx=' num2str(dx) ', dy=' num2str(dy)];
        set(htit,'string',str);
        htit=get(datar{2},'title');
        str=get(htit,'string');
        str{2}=['kx-k_ space, kxnyq=' num2str(kxnyq) ', kynyq=' num2str(kynyq)];
        set(htit,'string',str);
        axes(datar{2})
        xlabel('Crossline wavenumber');
        ylabel('Inline wavenumber');
        %flip axes if needed to match main fig
        xdir=get(hseis,'xdir');
        ydir=get(hseis,'ydir');
        set(datar{1},'xdir',xdir,'ydir',ydir)
end

%Make entry in windows list and set closerequestfcn
winname=['2D spectrum ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function showsvdsep(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

dname2={strrep(dname,'plotimage3D ... ',''), [mode ' ' num2str(thisone)]};
cmap=get(hthisfig,'colormap');
hseis=findobj(hthisfig,'tag','seismic');
gridon=get(hseis,'xgrid');
gridcolor=get(hseis,'gridcolor');
gridalpha=get(hseis,'gridalpha');
pos=get(hthisfig,'position');
datar=seisplotsvd_sep(data,y,x,dname2);
colormap(cmap);
if(strcmp(gridon,'on'))
    set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{3},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
end
set(gcf,'position',pos);
%flip axes if needed to match main fig
xdir=get(hseis,'xdir');
ydir=get(hseis,'ydir');
set([datar{1} datar{2} datar{3}],'xdir',xdir,'ydir',ydir)

%Make entry in windows list and set closerequestfcn
winname=['SVD separation ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function showsurf(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

dname2=[strrep(dname,'plotimage3D ... ',''), mode ' ' num2str(thisone)];
cmap=get(hthisfig,'colormap');
hseis=findobj(hthisfig,'tag','seismic');
% gridon=get(hseis,'xgrid');
% gridcolor=get(hseis,'gridcolor');
% gridalpha=get(hseis,'gridalpha');
pos=get(hthisfig,'position');
% datar=seisplotsvd_sep(data,y,x,dname2);
figure
surf(x,y,data);shading flat
hsurf=gca;
colormap(cmap);
set(gcf,'name',['Surface plot ' dname2])
% if(strcmp(gridon,'on'))
%     set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
%     set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
%     set(datar{3},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
% end
set(gcf,'position',pos);
%flip axes if needed to match main fig
xdir=get(hseis,'xdir');
ydir=get(hseis,'ydir');
set(hsurf,'xdir',xdir,'ydir',ydir)

%Make entry in windows list and set closerequestfcn
winname=['Surface plot ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function footprint(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

hbmap=findobj(gcf,'tag','basemap');
udat=get(hbmap,'userdata');
dx=udat{7};
dy=udat{8};

dname2={strrep(dname,'plotimage3D ... ',''), [mode ' ' num2str(thisone)]};
cmap=get(hthisfig,'colormap');
hseis=findobj(hthisfig,'tag','seismic');
gridon=get(hseis,'xgrid');
gridcolor=get(hseis,'gridcolor');
gridalpha=get(hseis,'gridalpha');
pos=get(hthisfig,'position');
datar=seisplotsvd_foot(data,x,y,dx,dy,dname2);
colormap(cmap);
if(strcmp(gridon,'on'))
    set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{3},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{4},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{5},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{6},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
end
set(gcf,'position',pos);
%flip axes if needed to match main fig
xdir=get(hseis,'xdir');
ydir=get(hseis,'ydir');
set([datar{1} datar{2} datar{3} datar{4} datar{5} datar{6}],'xdir',xdir,'ydir',ydir)

%Make entry in windows list and set closerequestfcn
winname=['Footprint ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function difference(~,~)
global PLOTIMAGE3DDIFFDIAL
if(~isempty(PLOTIMAGE3DDIFFDIAL))
    delete(PLOTIMAGE3DDIFFDIAL);
    PLOTIMAGE3DDIFFDIAL=[];
end
hfigs=getfigs;%get the grouped plotseis3D figures
hthisfig=gcf;
ind=hfigs~=hthisfig;
hotherfigs=hfigs(ind);%other figures in this group
%get this figures info
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
% data=udat{1};
% x=udat{2};
% y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};
nameA=[dname ' ' mode ' ' num2str(thisone)];
%get other figures in group
if(isempty(hotherfigs))
    msgbox('Difference plots require that you have a ''group'' of at least two datasets.','Oops!');
    return
else
    namesB=cell(size(hotherfigs));
    for k=1:length(hotherfigs)
        hprevious=findobj(hotherfigs(k),'tag','previous');
        udat=get(hprevious,'userdata');
        % data=udat{1};
        % x=udat{2};
        % y=udat{3};
        mode=udat{4};
        thisone=udat{5};
        dname=udat{6};
        namesB{k}=[dname ' ' mode ' ' num2str(thisone)];
    end
    
end
%create a dialog
pos=get(hthisfig,'position');
width=pos(3)*.5;
ht=pos(4)*.25;
xnow=pos(1)+.5*(pos(3)-width);
ynow=pos(2)+.5*(pos(4)-ht);
hdial=figure('position',[xnow,ynow,width,ht]);
xnow=.05;ynow=.6;
width=.3;ht=.05;
uicontrol(hdial,'style','text','String',['A: ' nameA],'units','normalized','tag','namea',...
    'position',[xnow,ynow,width,ht],'userdata',hthisfig)
xnow=.4;ynow=.8;
width=.5;ht=.05;
uicontrol(hdial,'style','text','String','Choose dataset B','units','normalized','tag','nameb',...
    'position',[xnow,ynow,width,ht])
ynow=.4;ht=.4;
uicontrol(hdial,'style','listbox','String',namesB,'units','normalized','tag','namesb',...
    'position',[xnow,ynow,width,ht],'userdata',hotherfigs)

xnow=.55;ynow=.1;
width=.2;ht=.2;
hbutgrp=uibuttongroup(hdial,'units','normalized','position',[xnow,ynow,width,ht],'tag','option',...
    'title','Choose subtraction order');
uicontrol(hbutgrp,'style','radio','string','A - B','units','normalized','position',[0,.5,1,.5],...
    'enable','on','tag','AB','backgroundcolor','w');
uicontrol(hbutgrp,'style','radio','string','B - A','units','normalized','position',[0,0,1,.5],...
    'enable','on','tag','BA','backgroundcolor','w');
%do-it button
xnow=.05;
ynow=.2;
width=.1;
ht=.1;
uicontrol(hdial,'style','pushbutton','string','Do it','tag','doit','units','normalized',...
    'position',[xnow,ynow,width,ht],'callback',@differencedoit,...
    'backgroundcolor','c','tooltipstring','Click to create the difference plot of the selected datasets');
%dismiss button
xnow=.05;
ynow=.1;
width=.1;
ht=.1;
uicontrol(hdial,'style','pushbutton','string','Dismiss','tag','dismiss','units','normalized',...
    'position',[xnow,ynow,width,ht],'callback','plotimage3D(''dismissdifference'');',...
    'backgroundcolor','r','tooltipstring','Click to dismiss this dialog');

PLOTIMAGE3DDIFFDIAL=hdial;
set(hdial,'closerequestfcn','plotimage3D(''dismissdifference'')','name','Plotimage3D Difference Dialog');

end

function differencedoit(~,~)
global PLOTIMAGE3DDIFFDIAL
hdial=PLOTIMAGE3DDIFFDIAL;
%determine the A and B datasets
hnamea=findobj(hdial,'tag','namea');
tmp=get(hnamea,'string');
namea=tmp(4:end);
hA=get(hnamea,'userdata');%figure for dataset A
cmap=get(hA,'colormap');
hnamesb=findobj(hdial,'tag','namesb');
namesb=get(hnamesb,'string');
iB=get(hnamesb,'value');
hotherfigs=get(hnamesb,'userdata');
nameb=namesb{iB};
hB=hotherfigs(iB);
%get the two datasets
hprevious=findobj(hA,'tag','previous');
udat=get(hprevious,'userdata');
seisa=udat{1};
xa=udat{2};
ya=udat{3};
hprevious=findobj(hB,'tag','previous');
udat=get(hprevious,'userdata');
seisb=udat{1};
xb=udat{2};
yb=udat{3};

if(length(xa)~=length(xb) || length(ya)~=length(yb))
    msgbox(['The selected datasets have different sizes and a difference is not possible. '...
        'If the section view is inline or crossline, make sure the tmin and tmax settings are '...
        'the same for both datasets and try again.'],'Oops!');
    return
end

%determine the subtraction order
hab=findobj(hdial,'tag','AB');
order=get(hab,'value');

if(order==1)
    datar=seisplotdiff(seisa,seisb,ya,xa,namea,nameb);
else
    datar=seisplotdiff(seisb,seisa,ya,xa,nameb,namea);
end
colormap(cmap)
pos=get(hA,'position');
set(gcf,'position',pos)
%flip axes if needed to match main fig
hseis=findobj(hA,'tag','seismic');
xdir=get(hseis,'xdir');
ydir=get(hseis,'ydir');
set([datar{1} datar{2} datar{3}],'xdir',xdir,'ydir',ydir)

%Make entry in windows list and set closerequestfcn
winname=['difference ' mode ' ' num2str(thisone)];
hwin=findobj(hA,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function showtvspectrum(~,~)
hthisfig=gcf;
hprevious=findobj(gcf,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
t=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};
% get current clim to impose on new display
cl=get(gca,'clim');

switch mode
    case 'inline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];

        pos=get(gcf,'position');
        datar=seisplottvs(data,t,x,dname2,nan,nan);
        axes(datar{1})
        xlabel('xline number');
        set(gca,'clim',cl);
        set(gcf,'position',pos);
        
    case 'xline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];

        pos=get(gcf,'position');
        datar=seisplottvs(data,t,x,dname2,nan,nan);
        axes(datar{1})
        xlabel('inline number');
        set(gca,'clim',cl);
        set(gcf,'position',pos);
    otherwise
        return;
        
end

%Make entry in windows list and set closerequestfcn
winname=['TVS ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);
end

function showtvphase(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
t=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

switch mode
    case 'inline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];

        pos=get(gcf,'position');
        haxes=seisplotphase(data,t,x,nan,nan,nan,20,dname2);
        axes(haxes{1})
        xlabel('xline number');
        axes(haxes{4})
        xlabel('xline number');
        set(gcf,'position',pos);
        
    case 'xline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];

        pos=get(gcf,'position');
        haxes=seisplotphase(data,t,x,nan,nan,nan,20,dname2);
        axes(haxes{1})
        xlabel('inline number');
        axes(haxes{4})
        xlabel('inline number');
        set(gcf,'position',pos);
    otherwise
        return;
        
end
%Make entry in windows list and set closerequestfcn
winname=['tv phase ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);
end

function showfxphase(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
thistag=get(gcf,'tag');
udat=get(hprevious,'userdata');
seis=udat{1};
x=udat{2};
t=udat{3};
mode=udat{4};
thisone=udat{5};%line number
dname=udat{6};
% get current clim to impose on new display
cl=get(gca,'clim');
fmax=nan;
%cmap=get(gcf,'colormap');
switch mode
    case 'inline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];
        pos=get(gcf,'position');
        flag=1;
        xname='xline';
        datar=seisplotfx(seis,t,x,dname2,nan,nan,fmax,xname,flag);
        axes(datar{1});
        xlabel('xline number');
        set(gca,'clim',cl);
        set(gcf,'position',pos);
        
    case 'xline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];

        pos=get(gcf,'position');
        flag=1;
        xname='iline';
        datar=seisplotfx(seis,t,x,dname2,nan,nan,fmax,xname,flag);
        axes(datar{1})
        xlabel('inline number');
        set(gca,'clim',cl);
        set(gcf,'position',pos);
        %colormap(cmap);
        
    otherwise
        return;
        
end

set(gcf,'tag',thistag);%this means that if pi3D window is 'fromsane' then so will be this window

%Make entry in windows list and set closerequestfcn
winname=['fx phase ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function showfxamp(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
thistag=get(gcf,'tag');
udat=get(hprevious,'userdata');
seis=udat{1};
x=udat{2};
t=udat{3};
mode=udat{4};
thisone=udat{5};%line number
dname=udat{6};
% get current clim to impose on new display
cl=get(gca,'clim');
fmax=nan;
%cmap=get(gcf,'colormap');
switch mode
    case 'inline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];
        pos=get(gcf,'position');
        flag=0;
        xname='xline';
        datar=seisplotfx(seis,t,x,dname2,nan,nan,fmax,xname,flag);
        axes(datar{1});
        xlabel('xline number');
        set(gca,'clim',cl);
        set(gcf,'position',pos);
        %colormap(cmap)
    case 'xline'
        dname2=[strrep(dname,'plotimage3D ... ','') ' ' mode ' ' num2str(thisone)];

        pos=get(gcf,'position');
        flag=0;
        xname='iline';
        datar=seisplotfx(seis,t,x,dname2,nan,nan,fmax,xname,flag);
        axes(datar{1})
        xlabel('inline number');
        set(gca,'clim',cl);
        set(gcf,'position',pos);
        %colormap(cmap)
    otherwise
        return;
        
end

set(gcf,'tag',thistag);%this means that if pi3D window is 'fromsane' then so will be this window

%Make entry in windows list and set closerequestfcn
winname=['fx amplitude ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function deconvolution(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

dname2=[strrep(dname,'plotimage3D ... ',''), [' ' mode ' ' num2str(thisone)]];
cmap=get(hthisfig,'colormap');
hseis=findobj(hthisfig,'tag','seismic');
gridon=get(hseis,'xgrid');
gridcolor=get(hseis,'gridcolor');
gridalpha=get(hseis,'gridalpha');
pos=get(hthisfig,'position');
datar=seisplotdecon(data,y,x,dname2);
colormap(cmap);
if(strcmp(gridon,'on'))
    set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
end
set(gcf,'position',pos);
%flip axes if needed to match main fig
xdir=get(hseis,'xdir');
ydir=get(hseis,'ydir');
set([datar{1} datar{2}],'xdir',xdir,'ydir',ydir)

%Make entry in windows list and set closerequestfcn
winname=['Spiking decon ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
crf=get(hfig,'closerequestfcn');
set(hfig,'closerequestfcn',[crf 'plotimage3D(''closewindow'');'],'userdata',hthisfig);
end

function filter(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
x=udat{2};
y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

dname2=[strrep(dname,'plotimage3D ... ',''), [' ' mode ' ' num2str(thisone)]];
cmap=get(hthisfig,'colormap');
hseis=findobj(hthisfig,'tag','seismic');
gridon=get(hseis,'xgrid');
gridcolor=get(hseis,'gridcolor');
gridalpha=get(hseis,'gridalpha');
pos=get(hthisfig,'position');
datar=seisplotfilt(data,y,x,dname2);
colormap(cmap);
if(strcmp(gridon,'on'))
    set(datar{1},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
    set(datar{2},'xgrid',gridon,'ygrid',gridon,'gridcolor',gridcolor,'gridalpha',gridalpha);
end
set(gcf,'position',pos);
%flip axes if needed to match main fig
xdir=get(hseis,'xdir');
ydir=get(hseis,'ydir');
set([datar{1} datar{2}],'xdir',xdir,'ydir',ydir)

%Make entry in windows list and set closerequestfcn
winname=['Bandpass filter ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function amphist(~,~)
hthisfig=gcf;
hprevious=findobj(hthisfig,'tag','previous');
udat=get(hprevious,'userdata');
data=udat{1};
% x=udat{2};
% y=udat{3};
mode=udat{4};
thisone=udat{5};
dname=udat{6};

dname2=[strrep(dname,'plotimage3D ... ',''), [' ' mode ' ' num2str(thisone)]];

pos=get(hthisfig,'position');

inonzero= data~=0.0;

figure('position',[pos(1:2) .5*pos(3) .5*pos(4)]);
hist(data(inonzero),200);
xlabel('amplitude');ylabel('number of samples');
ht=title(['Amplitude histogram for ' dname2]);
ht.Interpreter='none';
grid

%Make entry in windows list and set closerequestfcn
winname=['Amp histogram ' mode ' ' num2str(thisone)];
hwin=findobj(hthisfig,'tag','windows');
hfig=gcf;
currentwindows=get(hwin,'string');
currentfigs=get(hwin,'userdata');

nwin=length(currentwindows);
if(nwin==1)
   if(strcmp(currentwindows{1},'None'))
       currentwindows{1}=winname;
       currentfigs(1)=hfig;
       nwin=0;
   else
       currentwindows{2}=winname;
       currentfigs(2)=hfig;
   end
else
    currentwindows{nwin+1}=winname;
    currentfigs(nwin+1)=hfig;
end
set(hwin,'string',currentwindows,'value',nwin+1,'userdata',currentfigs)
set(hfig,'closerequestfcn','plotimage3D(''closewindow'')','userdata',hthisfig);

end

function tmin=get_tmin
htmin=findobj(gcf,'tag','tmin');
ival=get(htmin,'value');
tmins=get(htmin,'string');
tmin=str2double(tmins{ival});
end

function tmax=get_tmax
htmax=findobj(gcf,'tag','tmax');
ival=get(htmax,'value');
tmaxs=get(htmax,'string');
if(ival<=length(tmaxs))
    tmax=str2double(tmaxs{ival});
else
    tmax=str2double(tmaxs{end});
end
end

function flag=issane
%This function determines if the current figure was created by sane or not and returns a logical flag 
tag=get(gcf,'tag');
flag=false;
if(strcmp(tag,'fromsane'))
    udat=get(gcf,'userdata');
    if(iscell(udat))
        if(length(udat)==2)
            if(isgraphics(udat{2}))
                flag=true;
            end
        end
    end
end

end

function sd=sanedata
%this is used to communicate with SANE. When sane launches a plotimage3d
%window, it puts information in the userdata of the plotimage3D figure that
%is required for sane to know which dataset is sending the message. This
%function retrieves that data to send back to sane.
sd=get(gcf,'userdata');
end