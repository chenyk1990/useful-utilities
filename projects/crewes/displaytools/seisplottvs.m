function datar=seisplottvs(seis,t,x,dname,t1s,twins,fmax)
% SEISPLOTTVS: plots a seismic gather and its frequency spectrum in time windows
%
% datar=seisplottvs(seis,t,x,dname,t1s,twins,fmax)
%
% A new figure is created and divided into two same-sized axes (side-by-side). The seismic gather is
% plotted as an image in the left-hand-side and its temporal amplitude spectra in different time
% windows are plotted in the right-hand-side. Controls are provided to adjust the clipping and to
% brighten or darken the image plots.
%
% seis ... input seismic matrix
% t ... time coordinate vector for seis. This is the row coordinate of seis. 
% x ... space coordinate vector for seis
% dname ... text string giving a name for the dataset that will annotate
%       the plots.
% ************ default dname =[] ************
% t1s ... vector of 3 window start times (nan gets default)
% ********** default = [t(1) t(1)+twin t(2)+2*twin] where twin=(t(end)-t(1))/3 *********
% twins ... vector of 3 window lengths (nan gets default)
% ********** default = [twin twin twin] *************
% fmax ... maximum frequency to include on the frequency axis.
% ************ default = .5/(t(2)-t(1)) which is Nyquist ***********
%
% datar ... Return data which is a length 2 cell array containing
%           data{1} ... handle of the seismic axes
%           data{2} ... handle of the spectral axes
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
global SANE_TIMEWINDOWS
global FMAX DBLIM

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
    
    if(nargin<4)
        dname=[];
    end
    if(nargin<5)
        t1s=nan;
    end
    if(nargin<6)
        twins=nan;
    end
    
    if(isnan(t1s) && isnan(twins))
        if(~isempty(SANE_TIMEWINDOWS))
            t1s=SANE_TIMEWINDOWS(:,1);
            t2s=SANE_TIMEWINDOWS(:,2);
            twins=t2s-t1s;
        else
            twin=(t(end)-t(1))/3;
            t1s=[t(1)+.05*twin t(1)+twin t(1)+1.95*twin];
            twins=twin*ones(1,3);
        end
    end
    
    if(isnan(t1s))
        twin=(t(end)-t(1))/3;
        t1s=[t(1)+.05*twin t(1)+twin t(1)+1.95*twin];
    end
    if(isnan(twins))
        twin=(t(end)-t(1))/3;
        twins=twin*ones(1,3);
    end
    
    if(length(t1s)~=3 || length(twins)~=3)
        error('t1s and twins must be length 3');
    end
    
    fnyq=.5/(t(2)-t(1));
    
    if(nargin<7)
        if(isempty(FMAX))
            fmax=fnyq;
        else
            fmax=FMAX;
        end
    end
    
    if(fmax>fnyq)
        fmax=fnyq;
    end
    
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
    xmin=min(x);
    xmax=max(x);
    klrs=get(hax1,'colororder');
    lw=1;
    line([xmin xmax],[t1s(1) t1s(1)],'color',klrs(2,:),'linestyle','--','buttondownfcn','seisplottvs(''dragline'');','tag','1','linewidth',lw);
    line([xmin xmax],[t1s(1)+twins(1) t1s(1)+twins(1)],'color',klrs(2,:),'linestyle',':','buttondownfcn','seisplottvs(''dragline'');','tag','1b','linewidth',lw);
    line([xmin xmax],[t1s(2) t1s(2)],'color',klrs(3,:),'linestyle','--','buttondownfcn','seisplottvs(''dragline'');','tag','2','linewidth',lw);
    line([xmin xmax],[t1s(2)+twins(2) t1s(2)+twins(2)],'color',klrs(3,:),'linestyle',':','buttondownfcn','seisplottvs(''dragline'');','tag','2b','linewidth',lw);
    line([xmin xmax],[t1s(3) t1s(3)],'color',klrs(4,:),'linestyle','--','buttondownfcn','seisplottvs(''dragline'');','tag','3','linewidth',lw);
    line([xmin xmax],[t1s(3)+twins(3) t1s(3)+twins(3)],'color',klrs(4,:),'linestyle',':','buttondownfcn','seisplottvs(''dragline'');','tag','3b','linewidth',lw);
    
    %x boundary lines
    tmin=min(t);
    tmax=max(t);
    xdel=.01*(xmax-xmin);
    line([xmin+xdel xmin+xdel],[tmin tmax],'color','k','linestyle','--','buttondownfcn','seisplottvs(''dragline'');','tag','x1','linewidth',lw);
    line([xmax-xdel xmax-xdel],[tmin tmax],'color','k','linestyle','--','buttondownfcn','seisplottvs(''dragline'');','tag','x2','linewidth',lw);
    
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
    %make a button to reset time windows to the global values
    xnow=xnot+xwid;
    wid=.055;ht=.05;sep=.005;
    ynow=ynot+yht+sep;
    uicontrol(gcf,'style','pushbutton','string','Reset windows to globals','units','normalized',...
        'position',[xnow,ynow,1.5*wid,.5*ht],'callback','seisplottvs(''resetwindows'')','tag','resetwin',...
        'tooltipstring','Resets windows to the most recent published values');
    %make a clip control
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clipxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplottvs(''clipxt'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1,t1s,twins,fmax,dname},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brightenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplottvs(''brightenxt'');',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darkenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplottvs(''brightenxt'');',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightnessxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Info','tag','info','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplottvs(''info'');',...
        'tooltipstring','Click for gate adjustment instructions','userdata',0);
    
    set(hax1,'tag','seis');
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid yht]);
    tpad=2*max(twins);
    tvdbspec(t,seis,t1s,twins,tpad,dname,hax2);  
    yl=get(gca,'ylim');
    dblimmin=yl(1);
    if(~isempty(DBLIM))
        dblim=DBLIM;
    else
        dblim=dblimmin;
    end
    xlim([0 fmax])
    ylim([dblim 0])
    %make a clip control
    

    xnow=xnot+2*xwid+xsep;
    ht=.025;
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','pushbutton','string','recompute','tag','recompute','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplottvs(''recompute'');',...
        'tooltipstring','recompute the spectra');
    ynow=ynow-ht;
    uicontrol(gcf,'style','pushbutton','string','separate spectra','tag','separate','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplottvs(''separate'');',...
        'tooltipstring','separate the spectra for easier viewing','userdata',0);
     ynow=ynow-ht;
    uicontrol(gcf,'style','text','string','Fmax:','units','normalized',...
        'position',[xnow,ynow,.5*wid,ht],'tooltipstring','The maximum frequency to show');
    uicontrol(gcf,'style','edit','string',num2str(fmax),'units','normalized','tag','fmax',...
        'position',[xnow+.5*wid,ynow,.5*wid,ht],'tooltipstring','Enter a value in Hz.',...
        'callback','seisplottvs(''setlims'');','userdata',fnyq);
    ynow=ynow-ht;
    uicontrol(gcf,'style','text','string','db limit:','units','normalized',...
        'position',[xnow,ynow,.5*wid,ht],'tooltipstring','The minimum decibel level to show');
    uicontrol(gcf,'style','edit','string',num2str(dblim),'units','normalized','tag','dblim',...
        'position',[xnow+.5*wid,ynow,.5*wid,ht],'tooltipstring','Enter a negative number',...
        'callback','seisplottvs(''setlims'');','userdata',dblimmin);
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
    set(hax2,'tag','tvs');
    
    set(gcf,'name',['TVS analysis for ' dname]);
    if(nargout>0)
        datar=cell(1,2);
        datar{1}=hax1;
        datar{2}=hax2;
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
    t1s=udat{7};
    twins=udat{8};
    twin=.25*min(twins);
    
    h1=findobj(haxe,'tag','1');
    yy=get(h1,'ydata');
    t1=yy(1);
    h1b=findobj(haxe,'tag','1b');
    yy=get(h1b,'ydata');
    t1b=yy(1);
    h2=findobj(haxe,'tag','2');
    yy=get(h2,'ydata');
    t2=yy(2);
    h2b=findobj(haxe,'tag','2b');
    yy=get(h2b,'ydata');
    t2b=yy(2);
    h3=findobj(haxe,'tag','3');
    yy=get(h3,'ydata');
    t3=yy(1);
    h3b=findobj(haxe,'tag','3b');
    yy=get(h3b,'ydata');
    t3b=yy(1);
    hx1=findobj(haxe,'tag','x1');
    xx=get(hx1,'xdata');
    x1=xx(1);
    hx2=findobj(haxe,'tag','x2');
    xx=get(hx2,'xdata');
    x2=xx(1);
    
    hi=findobj(haxe,'type','image');
    t=get(hi,'ydata');
    x=get(hi,'xdata');
    tmin=t(1);tmax=t(end);
    xmin=min(x);xmax=max(x);
    xdel=.01*(xmax-xmin);
    DRAGLINE_SHOWPOSN='on';
    DRAGLINE_CALLBACK='';
    DRAGLINE_MOTIONCALLBACK='';
    if(hnow==h1)
        %clicked on t1
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[tmin t1b];
        DRAGLINE_PAIRED=h1b;
    elseif(hnow==h2)
        %clicked on t2
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[tmin t2b];
        DRAGLINE_PAIRED=h2b;
    elseif(hnow==h3)
        %clicked on t3
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[tmin t3b];
        DRAGLINE_PAIRED=h3b;
    elseif(hnow==h1b)
        %clicked on t1b
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[t1 tmax-twin];
        DRAGLINE_PAIRED=h1;
    elseif(hnow==h2b)
        %clicked on t2b
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[t2 tmax-twin];
        DRAGLINE_PAIRED=h2;
    elseif(hnow==h3b)
        %clicked on t3b
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[t3 tmax-twin];
        DRAGLINE_PAIRED=h3;
    elseif(hnow==hx1)
        %clicked on x1
        DRAGLINE_MOTION='xonly';
        DRAGLINE_XLIMS=[xmin+xdel x2-xdel];
        DRAGLINE_PAIRED=hx2;
    elseif(hnow==hx2)
        %clicked on x2
        DRAGLINE_MOTION='xonly';
        DRAGLINE_XLIMS=[x1+xdel xmax-xdel];
        DRAGLINE_PAIRED=hx1;
    end
    
    dragline('click')
    
elseif(strcmp(action,'resetwindows'))
    hclipxt=findobj(gcf,'tag','clipxt');
    udat=get(hclipxt,'userdata');
    haxe=udat{6};
    
    tglobal=SANE_TIMEWINDOWS;
    t1s=tglobal(:,1);
    t2s=tglobal(:,2);
    
    h1=findobj(haxe,'tag','1');
    set(h1,'ydata',[t1s(1) t1s(1)]);
    h1b=findobj(haxe,'tag','1b');
    set(h1b,'ydata',[t2s(1) t2s(1)]);
    h2=findobj(haxe,'tag','2');
    set(h2,'ydata',[t1s(2) t1s(2)]);
    h2b=findobj(haxe,'tag','2b');
    set(h2b,'ydata',[t2s(2) t2s(2)]);
    h3=findobj(haxe,'tag','3');
    set(h3,'ydata',[t1s(3) t1s(3)]);
    h3b=findobj(haxe,'tag','3b');
    set(h3b,'ydata',[t2s(3) t2s(3)]);
    
elseif(strcmp(action,'setlims'))
    hfmax=findobj(gcf,'tag','fmax');
    hdblim=findobj(gcf,'tag','dblim');
    tmp=get(hfmax,'string');
    fmax=str2double(tmp);
    fnyq=get(hfmax,'userdata');
    if(isnan(fmax) || fmax>fnyq || fmax<0)
        fmax=fnyq;
        set(hfmax,'string',num2str(fmax));
    end
    tmp=get(hdblim,'string');
    dblim=str2double(tmp);
    if(isnan(dblim))
        dblim=get(hdblim,'userdata');
        set(hdblim,'string',num2str(dblim));
    end
    if(dblim>0)
        dblim=-dblim;
        set(hdblim,'string',num2str(dblim));
    end
    htvs=findobj(gcf,'tag','tvs');
    axes(htvs);
    xlim([0 fmax]);
    ylim([dblim 0]);
    
    FMAX=fmax;
    DBLIM=dblim;
    
elseif(strcmp(action,'recompute'))
    hclipxt=findobj(gcf,'tag','clipxt');
    udat=get(hclipxt,'userdata');
    hax1=udat{6};
    t1s=udat{7};
    twins=udat{8};
    fmax=udat{9};
    dname=udat{10};
    dbflag=1;
    
    h1=findobj(hax1,'tag','1');
    yy=get(h1,'ydata');
    t1s(1)=yy(1);
    h1b=findobj(hax1,'tag','1b');
    yy=get(h1b,'ydata');
    twins(1)=yy(1)-t1s(1);
    
    h2=findobj(hax1,'tag','2');
    yy=get(h2,'ydata');
    t1s(2)=yy(1);
    h2b=findobj(hax1,'tag','2b');
    yy=get(h2b,'ydata');
    twins(2)=yy(1) -t1s(2);
    
    h3=findobj(hax1,'tag','3');
    yy=get(h3,'ydata');
    t1s(3)=yy(1);
    h3b=findobj(hax1,'tag','3b');
    yy=get(h3b,'ydata');
    twins(3)=yy(1)-t1s(3);
    
    hx1=findobj(hax1,'tag','x1');
    xx=get(hx1,'xdata');
    x1=xx(1);
    hx2=findobj(hax1,'tag','x2');
    xx=get(hx2,'xdata');
    x2=xx(1);
    
    
    udat{7}=t1s;
    udat{8}=twins;
    set(hclipxt,'userdata',udat);
    
    t2s=t1s+twins;
    
    SANE_TIMEWINDOWS=[t1s(:) t2s(:)];
    
    hi=findobj(hax1,'type','image');
    seis=get(hi,'cdata');
    t=get(hi,'ydata');
    x=get(hi,'xdata');
    indx=near(x,x1,x2);
    
    hax2=findobj(gcf,'tag','tvs');
    tpad=2*max(twins);
    tvdbspec(t,seis,t1s,twins,tpad,dname,hax2,dbflag,indx);
    set(hax2,'tag','tvs');
    boldlines(hax2,4);
    bigfont(hax2,1.5,1);
    xlim([0 fmax])
    
    hsep=findobj(gcf,'tag','separate');
    set(hsep,'string','separate spectra','userdata',0)
elseif(strcmp(action,'info'))
    msg=['Spectral analysis gates are shown on the seismic section with colors to match the ',...
        'spectral curves. There are 3 analysis gates with the gate top shown as a dashed line ',...
        'and the gate bottom as a dotted line. With the left mouse button, click on either the ',...
        'top or bottom of a gate and drag it to a new position. (The bottom cannot be dragged ',...
        'past the top and vice versa.) With the right mouse button, click and drag either the ',...
        'gate top or bottom and they will move together keeping the gate with constant. ',...
        'After adjusting the gates, click recompute to see the spectra. Avoid very narrow ',...
        'gates. Gates may overlap.'];
    pos=get(gcf,'position');
    msgbox(msg,'TVS gate adjustment instructions');
    pos2=get(gcf,'position');
    xc=pos(1)+.5*pos(3);
    yc=pos(2)+.5*pos(4);
    x1=xc-.5*pos2(3);
    y1=yc-.5*pos2(4);
    set(gcf,'position',[x1 y1 pos2(3:4)]);
elseif(strcmp(action,'separate'))
    hsep=gcbo;
    sep=get(hsep,'userdata');
    if(sep==0)
        %we are separating
        hax2=findobj(gcf,'tag','tvs');
        yl=get(hax2,'ylim');
        sep=round(abs(diff(yl))/10);
        hl=findobj(hax2,'type','line');
        yl=get(hl(1),'ydata');
        set(hl(1),'ydata',yl-3*sep);
        yl=get(hl(2),'ydata');
        set(hl(2),'ydata',yl-2*sep);
        yl=get(hl(3),'ydata');
        set(hl(3),'ydata',yl-sep);
        set(hsep,'userdata',sep);
        set(hsep,'string','combine spectra')
    else
        %we are un-separating
        hax2=findobj(gcf,'tag','tvs');
        hl=findobj(hax2,'type','line');
        yl=get(hl(1),'ydata');
        set(hl(1),'ydata',yl+3*sep);
        yl=get(hl(2),'ydata');
        set(hl(2),'ydata',yl+2*sep);
        yl=get(hl(3),'ydata');
        set(hl(3),'ydata',yl+sep);
        set(hsep,'userdata',0);
        set(hsep,'string','separate spectra')
    end
    
    
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