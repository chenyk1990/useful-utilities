function datar=seisplotphase(seis,t,x,twin,t1s,xref,delx,dname,fmin,fmax)
% SEISPLOTPHASE: plots a seismic gather and reletive phase profiles in 3 time windows
%
% datar=seisplotphase(seis,t,x,twin,t1s,xref,delx,dname,fmin,fmax)
%
% A new figure is created and divided into two same-sized axes (side-by-side). The seismic gather is
% plotted as an image in the left-hand-side and relative phase profiles in different time
% windows are plotted in the right-hand-side. A relative phase profile is created by choosing a
% reference trace and comparing it to all traces in the gather. For each comparison the best
% constant phase rotation is found that best matches the given trace to the reference trace.
%
% seis ... input seismic matrix
% t ... time coordinate vector for seis. This is the row coordinate of seis. 
% x ... space coordinate vector for seis
% twin ... window length (seconds) (nan gets default)
% ********** default = 1.0 *************
% t1s ... vector of 3 window center times (nan gets default)
% ********** default = [t(1)+.5*twin(1) t(1)+twin(1)+.5*twin(2) t(1)+twin(1)+twin(2)+.5*twin(3)] *********
% xref ... x coordinate of the reference trace  (nan gets default)
% ********** default .5*(x(1)+x(end));
% dname ... text string giving a name for the dataset that will annotate
%       the plots.
% ************ default dname =[] ************
% fmin ... frequency of high pass filter to be applied before measurements. Enter [] for no filter.
% ************ default = fmin=10 ***********
% fmax ... frequency of low-pass filter to be applied before measurements. Enter [] for no filter.
% ************ default fmax=60 **************
%
% datar ... Return data which is a length 2 cell array containing
%           data{1} ... handle of the seismic axes
%           data{2} ... handle of the phase axes
%
% G.F. Margrave, Devon, 2017
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

global DRAGLINE_MOTION DRAGLINE_XLIMS DRAGLINE_YLIMS DRAGLINE_SHOWPOSN DRAGLINE_CALLBACK DRAGLINE_MOTIONCALLBACK DRAGLINE_PAIRED

if(~ischar(seis))
    action='init';
else
    action=seis;
end

datar=[];%initialize return data to null

if(strcmp(action,'init'))
    
    
    if(length(t)~=size(seis,1))
        error('time coordinate vector does not match seismic');
    end
    if(length(x)~=size(seis,2))
        error('space coordinate vector does not match seismic');
    end
    if(nargin<10)
        fmax=60;
    end
    if(nargin<9)
        fmin=10;
    end
    if(nargin<8)
        dname=[];
    end
    if(nargin<7)
        delx=20;
    end
    if(nargin<4)
        twin=nan;
    end
    if(nargin<5)
        t1s=nan;
    end
    if(nargin<6)
        xref=nan;
    end
   
    if(isnan(xref))
        xref=.5*(x(1)+x(end));
    end
    
    if(isnan(twin))
        twin=1.0;
    end
    
    if(isnan(t1s))
        t1s=t(1)+[.5*twin, 1.5*twin, 2.5*twin];
    end
    
    if(length(t1s)~=3 )
        error('t1s must be length 3');
    end
    
    xmin=min(x);
    xmax=max(x);
    xs=[xmin xmax];
    t1s=[t1s(:) t1s(:)];
    
    xwid=.35;
    yht=.75;
    xsep=.1;
    xnot=.1;
    ynot=.1;
    

    figure
    hax1=subplot('position',[xnot ynot xwid yht]);

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis);
    if(iclip==1)
        clim=[-amax amax];
    else
        clim=[am-clip*sigma am+clip*sigma];
    end
        
    imagesc(x,t,seis,clim);colormap(seisclrs)
%     brighten(.5);
    grid
    title(dname ,'interpreter','none')
    
    %draw window start times
    
    klrs=get(hax1,'colororder');
    lw=1;
    line(xs,t1s(1,:),'color',klrs(2,:),'linestyle','--','buttondownfcn','seisplotphase(''dragline'');','tag','1','linewidth',lw);
    %line(xs,[t1s(1)+twins(1) t1s(1)+twins(1)],'color',klrs(2,:),'linestyle',':','buttondownfcn','seisplotphase(''dragline'');','tag','1b','linewidth',lw);
    line(xs,t1s(2,:),'color',klrs(3,:),'linestyle','--','buttondownfcn','seisplotphase(''dragline'');','tag','2','linewidth',lw);
    %line(xs,[t1s(2)+twins(2) t1s(2)+twins(2)],'color',klrs(3,:),'linestyle',':','buttondownfcn','seisplotphase(''dragline'');','tag','2b','linewidth',lw);
    line(xs,t1s(3,:),'color',klrs(4,:),'linestyle','--','buttondownfcn','seisplotphase(''dragline'');','tag','3','linewidth',lw);
    %line(xs,[t1s(3)+twins(3) t1s(3)+twins(3)],'color',klrs(4,:),'linestyle',':','buttondownfcn','seisplotphase(''dragline'');','tag','3b','linewidth',lw);
    
    %reference line
    yl=get(gca,'ylim');
    line([xref xref],yl,'color','k','linestyle','--','buttondownfcn','seisplotphase(''dragline'');','tag','4','linewidth',lw);
    
    maxmeters=7000;
    
    if(max(t)<10)
        ylabel('time (s)')
    elseif(max(t)<maxmeters)
        ylabel('depth (m)')
    else
        ylabel('depth (ft)')
    end
    if(max(x)<maxmeters)
        xlabel('distance (m)')
    else
        xlabel('distance (ft)')
    end
    
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid yht]);
    pos=get(hax2,'position');
    
    
    %make a clip control
    xnow=xnot+xwid;
    wid=.055;ht=.05;sep=.005;
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clipxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotphase(''clipxt'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1,t1s,twin,dname,xs,klrs(2:4,:),pos,x,xref,delx,fmin,fmax},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brightenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotphase(''brightenxt'');',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darkenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotphase(''brightenxt'');',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightnessxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Info','tag','info','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotphase(''info'');',...
        'tooltipstring','Click for gate adjustment instructions','userdata',0,'backgroundcolor','y');
    
    set(hax1,'tag','seis');
    

    haxes=tvphasepro(seis,t,x,xref,t1s,xs,twin,hax2,klrs(2:4,:),delx);
    
    %recompute button
    xnow=xnot+2*xwid+xsep;
    ht=.05;
    ynow=ynot+yht;
    uicontrol(gcf,'style','pushbutton','string','recompute','tag','recompute','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotphase(''recompute'');',...
        'tooltipstring','recompute the phases','userdata',haxes);
    %other phase controls
    ynow=ynow-.5*ht-sep;
    uicontrol(gcf,'style','text','string','Ensemble half width','units','normalized',...
        'position',[xnow,ynow,wid,.5*ht]);
    ynow=ynow-ht-sep;
    delxvec=[0,1,3,6,10,20,30,40,60,100];
    ind=find(delx==delxvec, 1);
    if(isempty(ind))
        tmp=[delxvec,delx];
        delxvec=sort(tmp);
    end
    delxstr=num2strcell(delxvec);
    idel=find(delxvec==delx);
    uicontrol(gcf,'style','popupmenu','string',delxstr,'tag','delx','units','normalized',...
        'position',[xnow+.25*wid,ynow,.5*wid,ht],'value',idel,...
        'tooltipstring',...
        'This is the half-width in traces of the ensemble to be compared to the reference trace');
    
    ynow=ynow-sep;
    uicontrol(gcf,'style','text','string','Temporal window','units','normalized',...
        'position',[xnow,ynow,wid,.5*ht]);
    ynow=ynow-ht-sep;
    tmax=min([t(end)-t(1), 2]);
    twinvec=.2:.2:tmax;
    ind=find(twin==twinvec, 1);
    if(isempty(ind))
        tmp=[twinvec,twin];
        twinvec=sort(tmp);
    end
    twinstr=num2strcell(twinvec);
    iwin=find(twinvec==twin);
    uicontrol(gcf,'style','popupmenu','string',twinstr,'tag','twin','units','normalized',...
        'position',[xnow+.25*wid,ynow,.5*wid,ht],'value',iwin,...
        'tooltipstring',...
        'This is the width in seconds of the temporal window at each analysis time');
    
    ynow=ynow-sep;
    uicontrol(gcf,'style','text','string','Fmin:','units','normalized',...
        'position',[xnow,ynow,.3*wid,.5*ht]);
    uicontrol(gcf,'style','edit','string',num2str(fmin),'tag','fmin','units','normalized',...
        'position',[xnow+.3*wid+sep,ynow,.5*wid,.5*ht],...
        'tooltipstring',...
        'Low-cut frequency (Hz), leave blank for no filter');
     ynow=ynow-.5*ht-sep;
    uicontrol(gcf,'style','text','string','Fmin:','units','normalized',...
        'position',[xnow,ynow,.3*wid,.5*ht]);
    uicontrol(gcf,'style','edit','string',num2str(fmax),'tag','fmax','units','normalized',...
        'position',[xnow+.3*wid+sep,ynow,.5*wid,.5*ht],...
        'tooltipstring',...
        'High-cut frequency (Hz), leave blank for no filter');
    
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
    set(gcf,'name',['Phase profiles for ' dname]);
    if(nargout>0)
        datar=cell(1,4);
        datar{1}=hax1;
        datar{2}=haxes(1);
        datar{3}=haxes(2);
        datar{4}=haxes(3);
    end
elseif(strcmp(action,'clipxt'))
    hclip=findobj(gcf,'tag','clipxt');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
   % amin=udat{5};
    sigma=udat{3};
    hax=udat{6};
    if(iclip==1)
        clim=[-amax amax];
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set(hax,'clim',clim);

elseif(strcmp(action,'brightenxt'))
    hbut=gcbo;
    hbright=findobj(gcf,'tag','brightenxt');
    if(hbut==hbright)
        inc=.1;
    else
        inc=-.1;
    end
    brighten(inc);
    hbrightness=findobj(gcf,'tag','brightnessxt');
    brightlvl=get(hbrightness,'userdata');
    brightlvl=brightlvl+inc;
    if(abs(brightlvl)<.01)
        brightlvl=0;
    end
    set(hbrightness,'string',['lvl ' num2str(brightlvl)],'userdata',brightlvl)

elseif(strcmp(action,'dragline'))
    hnow=gcbo;
    hclipxt=findobj(gcf,'tag','clipxt');
    udat=get(hclipxt,'userdata');
    haxe=udat{6};
    
    h1=findobj(haxe,'tag','1');
    yy=get(h1,'ydata');
    t1=yy(1);
   
    h2=findobj(haxe,'tag','2');
    yy=get(h2,'ydata');
    t2=yy(2);

    h3=findobj(haxe,'tag','3');
    yy=get(h3,'ydata');
    t3=yy(1);
    
    h4=findobj(haxe,'tag','4');

    
    hi=findobj(haxe,'type','image');
    t=get(hi,'ydata');
    tmin=t(1);tmax=t(end);
    DRAGLINE_SHOWPOSN='on';
    DRAGLINE_CALLBACK='';
    DRAGLINE_MOTIONCALLBACK='';
    if(hnow==h1)
        %clicked on t1
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[tmin t2];
        %DRAGLINE_PAIRED=h1b;
    elseif(hnow==h2)
        %clicked on t2
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[t1 t3];
        %DRAGLINE_PAIRED=h2b;
    elseif(hnow==h3)
        %clicked on t3
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[t2 tmax];
        %DRAGLINE_PAIRED=h3b;
    elseif(hnow==h4)
        DRAGLINE_MOTION='xonly';
        DRAGLINE_YLIMS=[];
        DRAGLINE_XLIMS=[];
    end
    
    dragline('click')
elseif(strcmp(action,'recompute'))
    hrecomp=gcbo;
    haxes=get(hrecomp,'userdata');
    hclipxt=findobj(gcf,'tag','clipxt');
    udat=get(hclipxt,'userdata');
    hax1=udat{6};
    t1s=udat{7};
    twin=udat{8};
    xs=udat{10};
    klrs=udat{11};
    pos=udat{12};
    x=udat{13};
    
    h1=findobj(hax1,'tag','1');
    yy=get(h1,'ydata');
    t1s(1)=yy(1);
    
    h2=findobj(hax1,'tag','2');
    yy=get(h2,'ydata');
    t1s(2)=yy(1);
    
    h3=findobj(hax1,'tag','3');
    yy=get(h3,'ydata');
    t1s(3)=yy(1);
    
    h4=findobj(hax1,'tag','4');
    xx=get(h4,'xdata');
    xref=xx(1);
    
    hdelx=findobj(gcf,'tag','delx');
    idel=get(hdelx,'value');
    delstr=get(hdelx,'string');
    delx=str2double(delstr{idel});
        
    htwin=findobj(gcf,'tag','twin');
    iwin=get(htwin,'value');
    twinstr=get(htwin,'string');
    twin=str2double(twinstr{iwin});
    
    udat{7}=t1s;
    udat{8}=twin;
    udat{14}=xref;
    udat{15}=delx;
    set(hclipxt,'userdata',udat);
    
    hi=findobj(hax1,'type','image');
    seis=get(hi,'cdata');
    t=get(hi,'ydata');
    delete(haxes);
    haxes=tvphasepro(seis,t,x,xref,t1s,xs,twin,pos,klrs,delx);
    boldlines(haxes,4);
    bigfont(haxes,1.5,1);
    set(hrecomp,'userdata',haxes);
    hsep=findobj(gcf,'tag','separate');
    set(hsep,'string','separate spectra','userdata',0)
elseif(strcmp(action,'info'))
    hclipxt=findobj(gcf,'tag','clipxt');
    udat=get(hclipxt,'userdata');
    twin=udat{8};
    delx=udat{15};
    msg=['Phase analysis times are shown on the seismic section with colors to match the ',...
        'resulting curves. There are 3 analysis times shown as dashed lines. ',...
        ['At each analysis time, a window of size ' num2str(twin) ' seconds centered on the time is defined.'],...
        'With the left mouse button, click on any of these lines ',...
        'and drag it to a new position. ',...
        'After adjusting the times, click recompute to generate new estimates. ',...
        'The phase estimates require the simultaneous estimate of time shift (or delay) ',...
        'and cross correlation. These results are also shown.'...
        ' Phase, delay, and lag are all computed relative to a reference trace which is indicated ',...
        [' by the black dashed line. Estmates were computed every delx=' int2str(delx) ' traces. '],...
        ['At each estimation location, an ensemble of ' int2str(2*delx+1) ' traces were compared to ']...
        'the reference trace to calculate the three statistics. ',...
        'The ensemble size is indicated by the gray patch in the phase axes.'];
    pos=get(gcf,'position');
    msgbox(msg,'TVPhase gate adjustment instructions');
    pos2=get(gcf,'position');
    xc=pos(1)+.5*pos(3);
    yc=pos(2)+.5*pos(4);
    x1=xc-.5*pos2(3);
    y1=yc-.5*pos2(4);
    set(gcf,'position',[x1 y1 pos2(3:4)]);
    
    
end
end

function [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(data)
% data ... input data
%
% 
% clips ... determined clip levels
% clipstr ... cell array of strings for each clip level for use in popup menu
% clip ... starting clip level
% iclip ... index into clips where clip is found
% sigma ... standard deviation of data
% am ... mean of data
% amax ... max of data
% amin ... min of data

sigma=std(data(:));
am=mean(data(:));
amin=min(data(:));
amax=max(data(:));
nsigma=ceil((amax-amin)/sigma);%number of sigmas that span the data

clips=[20 15 10 8 6 4 3 2 1 .1 .01 .001 .0001]';
if(nsigma<clips(1))
    ind= clips<nsigma;
    clips=[nsigma;clips(ind)];
else
    n=floor(log10(nsigma/clips(1))/log10(2));
    newclips=zeros(n,1);
    newclips(1)=nsigma;
    for k=n:-1:2
        newclips(k)=2^(n+1-k)*clips(1);
    end
    clips=[newclips;clips];
end

clipstr=cell(size(clips));
nclips=length(clips);
clipstr{1}='none';
for k=2:nclips
    clipstr{k}=['clip= ' num2str(sigfig(clips(k),3))];
end
iclip=near(clips,3);
clip=clips(iclip);

end