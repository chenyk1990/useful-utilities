function datar=seisplotdecon(seis1,t1,x1,dname1)
% SEISPLOTDECON: Interactive deconvolution of a seismic stack or gather
%
% datar=seisplotdecon(seis,t,x,dname)
%
% A new figure is created and divided into two same-sized axes (side-by-side). The input seismic
% gather is platted as an image in the left-hand-side and a deconvolved and bandpass filtered gather
% is plotted as an image in the right-hand-side. Initial display uses default parameters which will
% probably please no one. Controls are provided to adjust the deconvolution and filter and re-apply.
% The data should be regularly sampled in both t and x.
%
% seis ... seismic matrix
% t ... time coordinate vector for seis
% x ... space coordinate vector for seis
%   *********** default = 1:number_of_traces ************
% dname ... text string nameing the seismic matrix.
%   *********** default = 'Input data' **************
%
% datar ... Return data which is a length 2 cell array containing
%           data{1} ... handle of the input seismic axes
%           data{2} ... handle of the filter seismic axes
% These return data are provided to simplify plotting additional lines and
% text in either axes.
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
global DECONGATE_TOP DECONGATE_BOT DECON_OP DECON_STAB DECON_FMIN DECON_FMAX

if(~ischar(seis1))
    action='init';
else
    action=seis1;
end

datar=[];%initialize return data to null

if(strcmp(action,'init'))
    
    if(nargin<2)
        error('at least 3 inputs are required');
    end
    if(nargin<3)
        x1=1:size(seis1,2);
    end
    if(nargin<4)
        dname1='Input data';
    end
    
    x2=x1;
    t2=t1;
    dt=t1(2)-t1(1);
    fnyq=.5/dt;
    if(isempty(DECON_FMAX))
        fmax=round(.25*fnyq);
    else
        fmax=DECON_FMAX;
    end
    if(isempty(DECON_FMIN))
        fmin=10;
    else
        fmin=DECON_FMIN;
    end
    if(isempty(DECON_OP))
        top=.1;
    else
        top=DECON_OP;
    end
    if(isempty(DECON_STAB))
        stab=.001;
    else
        stab=DECON_STAB;
    end

    
    %seis2=filter_stack(seis1,t1,fmin,fmax,'method','filtf');
    %dname2=[dname1 ' filtered, fmin=' num2str(fmin) ', fmax=' num2str(fmax)];
    seis2=seis1;
    
    if(length(t1)~=size(seis1,1))
        error('time coordinate vector does not match first seismic matrix');
    end
    if(length(x1)~=size(seis1,2))
        error('space coordinate vector does not match first seismic matrix');
    end
    if(length(t2)~=size(seis2,1))
        error('time coordinate vector does not match second seismic matrix');
    end
    if(length(x2)~=size(seis2,2))
        error('space coordinate vector does not match second seismic matrix');
    end
    
    if(iscell(dname1))
        dname1=dname1{1};
    end

    xwid=.35;
    yht=.8;
    xsep=.1;
    xnot=.1;
    ynot=.1;
    

    figure
    hax1=subplot('position',[xnot ynot xwid yht]);

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis1);
    clim=[am-clip*sigma am+clip*sigma];
        
    imagesc(x1,t1,seis1,clim);colormap(seisclrs)
    brighten(.5);
    grid
    ht=title(dname1);
    ht.Interpreter='none';
    maxmeters=7000;
    if(max(t1)<10)
        ylabel('time (s)')
    elseif(max(t1)<maxmeters)
        ylabel('depth (m)')
    else
        ylabel('depth (ft)')
    end
    xlabel('line coordinate')
    
    %draw decon gate
    if(~isempty(DECONGATE_TOP))
        ttop=DECONGATE_TOP*ones(1,2);
        tbot=DECONGATE_BOT*ones(1,2);
    else
        trange=t1(end)-t1(1);
        ttop=t1(1)+.25*trange*ones(1,2);
        tbot=ttop+.5*trange;
    end
    lw=.5;
    xs=[x1(1) x1(end)];
    h1=line(xs,ttop,'color','r','linestyle','-','buttondownfcn','seisplotdecon(''dragline'');','tag','ttop','linewidth',lw);
    h2=line(xs,tbot,'color','r','linestyle','--','buttondownfcn','seisplotdecon(''dragline'');','tag','tbot','linewidth',lw);
    
    legend([h1 h2],'design gate top','design gate bottom','location','southeast')
    
    %set gates to published value
    xnow=xnot+xwid;
    ynow=ynot+yht;
    wid=.055;ht=.05;sep=.005;
    uicontrol(gcf,'style','pushbutton','string','Use published gate','tag','setgate','units','normalized',...
        'position',[xnow ynow wid .5*ht],'callback','seisplotdecon(''setgate'');',...
        'tooltipstring','Sets the decon gate to the last published value');
    
    %make a clip control
    xnow=xnot+xwid;
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip1','units','normalized',...
        'position',[xnow,ynow,1.2*wid,ht],'callback','seisplotdecon(''clip1'')','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brighten','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''brighten'')',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darken','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''brighten'')',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightness','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    
    set(hax1,'tag','seis1');
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid yht]);

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis2);
    clim=[am-clip*sigma am+clip*sigma];
        
    imagesc(x2,t2,seis2,clim);colormap(seisclrs)
    brighten(.5);
    grid
%     dname2=dname1;
%     ht=title(dname2);
%     ht.Interpreter='none';
    
    if(max(t2)<10)
        ylabel('time (s)')
    elseif(max(t2)<maxmeters)
        ylabel('depth (m)')
    else
        ylabel('(depth (ft)')
    end
%     if(max(x2)<maxmeters)
%         xlabel('distance (m)')
%     else
%         xlabel('distance (ft)')
%     end
    xlabel('line coordinate')
    %make a clip control

    xnow=xnot+2*xwid+xsep;
    ht=.05;
    ynow=ynot+yht-ht;
    %wid=.045;sep=.005;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip2','units','normalized',...
        'position',[xnow,ynow,1.2*wid,ht],'callback','seisplotdecon(''clip2'')','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax2},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    %decon parameters
    ht=.025;
    ynow=ynow-ht-sep;
    xnow=xnow+sep;
    uicontrol(gcf,'style','text','string','Decon parameters:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring','These are for spiking decon (Wiener, aka deconw)');
    ynow=ynow-ht-sep;
    wid=wid*.5;
    uicontrol(gcf,'style','text','string','oplen:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring','decon operator length in seconds');
    uicontrol(gcf,'style','edit','string',num2str(top),'units','normalized','tag','oplen',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring','Enter a value in seconds between 0 and 1');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','stab:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring','stability or white noise constant');
    uicontrol(gcf,'style','edit','string',num2str(stab),'units','normalized','tag','stab',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring','Enter a value between 0 and 1');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Apply Decon','units','normalized',...
        'position',[xnow,ynow,2*wid,ht],'callback','seisplotdecon(''applydecon'');',...
        'tooltipstring','Apply current decon and filter specs','tag','deconbutton');
    
    %filter parameters
    ynow=ynow-2*ht-sep;
    xnow=xnow+sep;
    wid=wid*2;
    uicontrol(gcf,'style','text','string','Filter parameters:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring','These are for a post-decon bandpass');
    ynow=ynow-ht-sep;
    wid=wid*.3;
    uicontrol(gcf,'style','text','string','Fmin:','units','normalized',...
        'position',[xnow,ynow,wid,ht],'tooltipstring',...
        'This is the minimum frequency (Hz) to pass, enter zero for a lowpass filter');
    uicontrol(gcf,'style','edit','string',num2str(fmin),'units','normalized','tag','fmin',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring',['Enter a value in Hz between 0 and ' num2str(fnyq)]);
     uicontrol(gcf,'style','text','string','dFmn:','units','normalized',...
        'position',[xnow+2*(wid+sep),ynow,wid,ht],'tooltipstring',...
        'This is the rolloff width on the lowend. Leave blank for the default which is .5*Fmin');
    uicontrol(gcf,'style','edit','string','','units','normalized','tag','dfmin',...
        'position',[xnow+3*(wid+sep),ynow,wid,ht],'tooltipstring','Enter a value in Hz between 0 and Fmin');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','Fmax:','units','normalized',...
        'position',[xnow,ynow,wid,ht]);
    uicontrol(gcf,'style','edit','string',num2str(fmax),'units','normalized','tag','fmax',...
        'position',[xnow+wid+sep,ynow,wid,ht],'tooltipstring',['Enter a value in Hz between 0 and ' num2str(fnyq)]);
    uicontrol(gcf,'style','text','string','dFmx:','units','normalized',...
        'position',[xnow+2*(wid+sep),ynow,wid,ht],'tooltipstring',...
        'This is the rolloff width on the high end. Leave blank for the default which is 10 Hz');
    uicontrol(gcf,'style','edit','string','','units','normalized','tag','dfmax',...
        'position',[xnow+3*(wid+sep),ynow,wid,ht],'tooltipstring','Enter a value in Hz between 0 and Fnyq-Fmax');
    ynow=ynow-ht-sep;
    wid=0.03;
    uicontrol(gcf,'style','text','string','Phase:','units','normalized',...
        'position',[xnow,ynow,wid,ht]);
    uicontrol(gcf,'style','popupmenu','string',{'zero','minimum'},'units','normalized','tag','phase',...
        'position',[xnow+wid+sep,ynow,1.3*wid,ht],'tooltipstring','Usually choose zero');
    ynow=ynow-ht-sep;
    wid=0.055;
    uicontrol(gcf,'style','pushbutton','string','Apply Filter','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''applyfilter'');',...
        'tooltipstring','Apply current filter specs');
    %spectra
    ynow=ynow-2*ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Show spectra','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''spectra'');',...
        'tooltipstring','Show spectra in separate window','tag','spectra','userdata',[]);
    
    
    ynow=ynow-2*ht-sep;
     uicontrol(gcf,'style','text','string','Compute performace:','units','normalized',...
        'position',[xnow,ynow,1.5*wid,ht],'tooltipstring','For decon only');
    ynow=ynow-ht-sep;
     uicontrol(gcf,'style','text','string','','units','normalized','tag','performance',...
        'position',[xnow,ynow,1.5*wid,ht]);
    
    %zoom buttons
    wid=.1;
    pos=get(hax1,'position');
    xnow=pos(1)+.5*pos(3)-.5*wid;
    ynow=.97;
    uicontrol(gcf,'style','pushbutton','string','Zoom #1 like #2','units','normalized',...
        'position',[xnow ynow wid ht],'tag','1like2','callback','seisplotdecon(''equalzoom'');');
    
    pos=get(hax2,'position');
    xnow=pos(1)+.5*pos(3)-.5*wid;
    uicontrol(gcf,'style','pushbutton','string','Zoom #2 like #1','units','normalized',...
        'position',[xnow ynow wid ht],'tag','2like1','callback','seisplotdecon(''equalzoom'');');
    
    %results popup
    xnow=pos(1);
    ynow=pos(2)+pos(4)-ht;
    wid=pos(3);
    ht=3*ht;
    fs=14;
    uicontrol(gcf,'style','popupmenu','string','Diddley','units','normalized','tag','results',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''select'');','fontsize',fs,...
        'fontweight','bold')
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
    set(hax2,'tag','seis2');
    seisplotdecon('applydecon');
%     if(iscell(dname2))
%         dn2=dname2{1};
%     else
%         dn2=dname2;
%     end
    set(gcf,'name',['Spiking decon analysis for ' dname1],'closerequestfcn','seisplotdecon(''close'');');
    if(nargout>0)
        datar=cell(1,2);
        datar{1}=hax1;
        datar{2}=hax2;
    end
elseif(strcmp(action,'clip1'))
    hclip=findobj(gcf,'tag','clip1');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
    amin=udat{5};
    sigma=udat{3};
    hax=udat{6};
    if(iclip==1)
        clim=[amin amax];
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set(hax,'clim',clim);
elseif(strcmp(action,'clip2'))
    hclip=findobj(gcf,'tag','clip2');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
    amin=udat{5};
    sigma=udat{3};
    hax=udat{6};
    if(iclip==1)
        %clim=[amin amax];
        clim=[amin amax];
    else
        clip=clips(iclip-1);
        clim=[am-clip*sigma,am+clip*sigma];
        %clim=[amin am+clip*sigma];
    end
    set(hax,'clim',clim);
elseif(strcmp(action,'brighten'))
    hbut=gcbo;
    hbright=findobj(gcf,'tag','brighten');
    if(hbut==hbright)
        inc=.1;
    else
        inc=-.1;
    end
    brighten(inc);
    hbrightness=findobj(gcf,'tag','brightness');
    brightlvl=get(hbrightness,'userdata');
    brightlvl=brightlvl+inc;
    if(abs(brightlvl)<.01)
        brightlvl=0;
    end
    set(hbrightness,'string',['lvl ' num2str(brightlvl)],'userdata',brightlvl)
elseif(strcmp(action,'equalzoom'))
    hbut=gcbo;
    hseis1=findobj(gcf,'tag','seis1');
    hseis2=findobj(gcf,'tag','seis2');
    tag=get(hbut,'tag');
    switch tag
        case '1like2'
            xl=get(hseis2,'xlim');
            yl=get(hseis2,'ylim');
            set(hseis1,'xlim',xl,'ylim',yl);
            
        case '2like1'
            xl=get(hseis1,'xlim');
            yl=get(hseis1,'ylim');
            set(hseis2,'xlim',xl,'ylim',yl);
    end
elseif(strcmp(action,'dragline'))
    hnow=gcbo;
%     hclipxt=findobj(gcf,'tag','clipxt');
%     udat=get(hclipxt,'userdata');
%     haxe=udat{6};
%     t1s=udat{7};
%     twins=udat{8};
%     twin=.25*min(twins);
    
    hseis1=findobj(gcf,'tag','seis1');

    h1=findobj(hseis1,'tag','ttop');
    yy=get(h1,'ydata');
    ttop=yy(1);
   
    h2=findobj(hseis1,'tag','tbot');
    yy=get(h2,'ydata');
    tbot=yy(2);

    
    hi=findobj(hseis1,'type','image');
    t=get(hi,'ydata');
    tmin=t(1);tmax=t(end);
    DRAGLINE_SHOWPOSN='on';
    DRAGLINE_CALLBACK='';
    DRAGLINE_MOTIONCALLBACK='';
    if(hnow==h1)
        %clicked on ttop
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[tmin tbot];
        %DRAGLINE_PAIRED=h1b;
    elseif(hnow==h2)
        %clicked on tbot
        DRAGLINE_MOTION='yonly';
        DRAGLINE_YLIMS=[ttop tmax];
        %DRAGLINE_PAIRED=h2b;
    end
    
    dragline('click')
elseif(strcmp(action,'applydecon'))
    %plan: apply the decon parameters and update the performace label. Then put the result in
    %userdata of the decon button and call 'apply filter'. Apply filter will produce the label and
    %the saved result. The most-recent decon without a filter remains in the button's user data so
    %that a different filter can be applied. Save results will always have both decon and filter
    hseis1=findobj(gcf,'tag','seis1');
    hi=findobj(hseis1,'type','image');
    seis=get(hi,'cdata');
    t=get(hi,'ydata');
    %get the design gate
    htop=findobj(hseis1,'tag','ttop');
    yy=get(htop,'ydata');
    ttop=yy(1);
    hbot=findobj(hseis1,'tag','tbot');
    yy=get(hbot,'ydata');
    tbot=yy(1);
    idesign=near(t,ttop,tbot);
    %get the operator length 
    hop=findobj(gcf,'tag','oplen');
    val=get(hop,'string');
    top=str2double(val);
    if(isnan(top))
        msgbox('oplen is not recognized as a number','Oh oh ...');
        return;
    end
    if(top<0 || top>1)
        msgbox('oplen is unreasonable, enter a value in seconds');
        return;
    end
    %get the stab 
    hstab=findobj(gcf,'tag','stab');
    val=get(hstab,'string');
    stab=str2double(val);
    if(isnan(stab))
        msgbox('stab is not recognized as a number','Oh oh ...');
        return;
    end
    if(stab<0 || stab>1)
        msgbox('stab is unreasonable, enter a value between 0 and 1');
        return;
    end
    %deconvolve
    t1=clock;
    seisd=deconw_stack(seis,t,0,ttop,tbot,1,top,stab);
    DECONGATE_TOP=ttop;
    DECONGATE_BOT=tbot;
    DECON_OP=top;
    DECON_STAB=stab;
    t2=clock;
    timepertrace=round(100000*etime(t2,t1)/size(seis,2))/1000;
    hperf=findobj(gcf,'tag','performance');
    set(hperf,'string',[num2str(timepertrace) ' ms/trace'])
    hdbut=findobj(gcf,'tag','deconbutton');
    set(hdbut,'userdata',{seisd,ttop,tbot,top,stab});
    seisplotdecon('applyfilter');
elseif(strcmp(action,'applyfilter'))
    hdbut=findobj(gcf,'tag','deconbutton');
    udat=get(hdbut,'userdata');
    seisd=udat{1};
    ttop=udat{2};
    tbot=udat{3};
    top=udat{4};
    stab=udat{5};
    hseis2=findobj(gcf,'tag','seis2');
    hi=findobj(hseis2,'type','image');
    t=get(hi,'ydata');
    fnyq=.5/(t(2)-t(1));
    hobj=findobj(gcf,'tag','fmin');
    val=get(hobj,'string');
    fmin=str2double(val);
    if(isnan(fmin))
        msgbox('Fmin is not recognized as a number','Oh oh ...');
        return;
    end
    if(fmin<0 || fmin>fnyq)
        msgbox(['Fmin must be greater than 0 and less than ' num2str(fnyq)],'Oh oh ...');
        return;
    end
    hobj=findobj(gcf,'tag','dfmin');
    val=get(hobj,'string');
    if(~isempty(val))
        dfmin=str2double(val);
        if(isnan(dfmin))
            msgbox('dFmin is not recognized as a number','Oh oh ...');
            return;
        end
        if(dfmin<0 || dfmin>fmin)
            msgbox(['dFmin must be greater than 0 and less than ' num2str(fmin)],'Oh oh ...');
            return;
        end
    else
        dfmin=.5*fmin;
    end
    hobj=findobj(gcf,'tag','fmax');
    val=get(hobj,'string');
    fmax=str2double(val);
    if(isnan(fmax))
        msgbox('Fmax is not recognized as a number','Oh oh ...');
        return;
    end
    if(fmax<0 || fmax>fnyq)
        msgbox(['Fmax must be greater than 0 and less than ' num2str(fnyq)],'Oh oh ...');
        return;
    end
    if(fmax<=fmin && fmax~=0)
        msgbox('Fmax must be greater than Fmin','Oh oh ...');
        return;
    end
    hobj=findobj(gcf,'tag','dfmax');
    val=get(hobj,'string');
    if(~isempty(val))
        dfmax=str2double(val);
        if(isnan(dfmax))
            msgbox('dFmax is not recognized as a number','Oh oh ...');
            return;
        end
        if(dfmax<0 || dfmax>fnyq-fmax)
            msgbox(['dFmax must be greater than 0 and less than ' num2str(fnyq-fmax)],'Oh oh ...');
            return;
        end
    else
        dfmax=10;
    end
    hobj=findobj(gcf,'tag','phase');
    ival=get(hobj,'value');
    phase=ival-1;
    
    seis2=filter_stack(seisd,t,fmin,fmax,'method','filtf','phase',phase,'dflow',dfmin,'dfhigh',dfmax);
    DECON_FMIN=fmin;
    DECON_FMAX=fmax;
    set(hi,'cdata',seis2);
    axes(hseis2);
    dname=['Decon oplen=' num2str(top) ', stab=' num2str(stab) ', gate ' time2str(ttop) '-' time2str(tbot)];
    name=[dname ', & [' num2str(fmin) ',' num2str(dfmin) ']-[' num2str(fmax) ',' num2str(dfmax) '] filter'];
    %update clipping
    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis2);
    clim=[am-clip*sigma am+clip*sigma];
    hclip2=findobj(gcf,'tag','clip2');
    set(hclip2,'string',clipstr','value',iclip,'userdata',{clips,am,sigma,amax,amin,hseis2});
    set(hseis2,'clim',clim);
    seisplotdecon('clip2');
    %save the results and update hresults
    hresults=findobj(gcf,'tag','results');
    results=get(hresults,'userdata');
    if(isempty(results))
        nresults=1;
        results.names={name};
        results.data={seis2};
        results.top={top};
        results.ttop={ttop};
        results.tbot={tbot};
        results.stab={stab};
        results.fmins={fmin};
        results.dfmins={dfmin};
        results.fmaxs={fmax};
        results.dfmaxs={dfmax};
        results.phases={phase};
    else
        nresults=length(results.names)+1;
        results.names{nresults}=name;
        results.data{nresults}=seis2;
        results.top{nresults}=top;
        results.ttop{nresults}=ttop;
        results.tbot{nresults}=tbot;
        results.stab{nresults}=stab;
        results.fmins{nresults}=fmin;
        results.dfmins{nresults}=dfmin;
        results.fmaxs{nresults}=fmax;
        results.dfmaxs{nresults}=dfmax;
        results.phases{nresults}=phase;
    end
    set(hresults,'string',results.names,'value',nresults,'userdata',results)
    
    %see if spectra window is open
    hspec=findobj(gcf,'tag','spectra');
    hspecwin=get(hspec,'userdata');
    if(isgraphics(hspecwin))
        seisplotdecon('spectra');
    end
    
elseif(strcmp(action,'spectra'))
    hfig=gcf;
    name=get(hfig,'name');
    ind=strfind(name,'Spectral display');
    if(isempty(ind))
        hmaster=hfig;
    else
        hmaster=get(hfig,'userdata');
    end
    hseis1=findobj(hmaster,'tag','seis1');
    hseis2=findobj(hmaster,'tag','seis2');
    hi=findobj(hseis1,'type','image');
    seis1=get(hi,'cdata');
    hi=findobj(hseis2,'type','image');
    seis2=get(hi,'cdata');
    t=get(hi,'ydata');
    hspec=findobj(hmaster,'tag','spectra');
    hspecwin=get(hspec,'userdata');
    if(isempty(hspecwin))
        %make the spectral window if it does not already exist
        pos=get(hmaster,'position');
        wid=pos(3)*.5;ht=pos(4)*.5;
        x0=pos(1)+pos(3)-wid;y0=pos(2);
        hspecwin=figure('position',[x0,y0,wid,ht],'closerequestfcn','seisplotdecon(''closespec'');','userdata',hmaster);
        set(hspecwin,'name','Spectral display window')
        
        whitefig;
        x0=.1;y0=.1;awid=.7;aht=.8;
        subplot('position',[x0,y0,awid,aht]);
        sep=.01;
        ht=.05;wid=.075;
        ynow=y0+aht-ht;
        xnow=x0+awid+sep;
        uicontrol(gcf,'style','text','string','tmin:','units','normalized',...
            'position',[xnow,ynow,wid,ht])
        ntimes=10;
        tinc=round(10*(t(end)-t(1))/ntimes)/10;
        %times=[fliplr(0:-tinc:t(1)) tinc:tinc:t(end)-tinc];
        times=t(1):tinc:t(end)-tinc;
        %times=t(1):tinc:t(end)-tinc;
        stimes=num2strcell(times);
        ynow=ynow-ht-sep;
        uicontrol(gcf,'style','popupmenu','string',stimes,'units','normalized','tag','tmin',...
            'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''spectra'');','userdata',times);
        ynow=ynow-ht-sep;
        uicontrol(gcf,'style','text','string','tmax:','units','normalized',...
            'position',[xnow,ynow,wid,ht])
        times=t(end):-tinc:tinc;
        stimes=num2strcell(times);
        ynow=ynow-ht-sep;
        uicontrol(gcf,'style','popupmenu','string',stimes,'units','normalized','tag','tmax',...
            'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''spectra'');','userdata',times);
        ynow=ynow-ht-sep;
        uicontrol(gcf,'style','text','string','db range:','units','normalized',...
            'position',[xnow,ynow,wid,ht])
        db=-20:-20:-160;
        idb=near(db,-100);
        dbs=num2strcell(db);
        ynow=ynow-ht-sep;
        uicontrol(gcf,'style','popupmenu','string',dbs,'units','normalized','tag','db','value',idb,...
            'position',[xnow,ynow,wid,ht],'callback','seisplotdecon(''spectra'');','userdata',db);
        set(hspec,'userdata',hspecwin);
    else
        figure(hspecwin);
    end
    htmin=findobj(gcf,'tag','tmin');
    times=get(htmin,'userdata');
    it=get(htmin,'value');
    tmin=times(it);
    htmax=findobj(gcf,'tag','tmax');
    times=get(htmax,'userdata');
    it=get(htmax,'value');
    tmax=times(it);
    ind=near(t,tmin,tmax);
    hdb=findobj(gcf,'tag','db');
    db=get(hdb,'userdata');
    dbmin=db(get(hdb,'value'));
    pct=10;
    [S1,f]=fftrl(seis1(ind,:),t(ind),pct);
    S2=fftrl(seis2(ind,:),t(ind),pct);
    A1=mean(abs(S1),2);
    A2=mean(abs(S2),2);
    hh=plot(f,todb(A1),f,todb(A2));
    set(hh,'linewidth',2)
    xlabel('Frequency (Hz)')
    ylabel('decibels');
    ylim([dbmin 0])
    grid on
    legend('Input','Decon+filter'); 
    title(['Average ampltude spectra, tmin=' time2str(tmin) ', tmax=' time2str(tmax)]);
elseif(strcmp(action,'closespec'))
    hfig=gcf;
    hdaddy=get(hfig,'userdata');
    hspec=findobj(hdaddy,'tag','spectra');
    set(hspec,'userdata',[]);
    delete(hfig);
    if(isgraphics(hdaddy))
        figure(hdaddy);
    end
elseif(strcmp(action,'select'))
    hfig=gcf;
    hresults=findobj(hfig,'tag','results');
    results=get(hresults,'userdata');
    iresult=get(hresults,'value');
    hseis2=findobj(hfig,'tag','seis2');
    hi=findobj(hseis2,'type','image');
    set(hi,'cdata',results.data{iresult});
    hop=findobj(hfig,'tag','oplen');
    set(hop,'string',num2str(results.top{iresult}));
    hstab=findobj(hfig,'tag','stab');
    set(hstab,'string',num2str(results.stab{iresult}));
    httop=findobj(hfig,'tag','ttop');
    set(httop,'ydata',ones(1,2)*results.ttop{iresult});
    htbot=findobj(hfig,'tag','tbot');
    set(htbot,'ydata',ones(1,2)*results.tbot{iresult});
    hfmin=findobj(hfig,'tag','fmin');
    set(hfmin,'string',num2str(results.fmins{iresult}));
    hdfmin=findobj(hfig,'tag','dfmin');
    set(hdfmin,'string',num2str(results.dfmins{iresult}));
    hfmax=findobj(hfig,'tag','fmax');
    set(hfmax,'string',num2str(results.fmaxs{iresult}));
    hdfmax=findobj(hfig,'tag','dfmax');
    set(hdfmax,'string',num2str(results.dfmaxs{iresult}));
    hphase=findobj(hfig,'tag','phase');
    set(hphase,'value',results.phases{iresult}+1);
    %update clipping
    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(results.data{iresult});
    clim=[am-clip*sigma am+clip*sigma];
    hclip2=findobj(gcf,'tag','clip2');
    set(hclip2,'string',clipstr','value',iclip,'userdata',{clips,am,sigma,amax,amin,hseis2});
    set(hseis2,'clim',clim);
    seisplotdecon('clip2')
    %see if spectra window is open
    hspec=findobj(hfig,'tag','spectra');
    hspecwin=get(hspec,'userdata');
    if(isgraphics(hspecwin))
        seisplotdecon('spectra');
    end
elseif(strcmp(action,'setgate'))
    if(isempty(DECONGATE_TOP))
        return
    end
    ttop=DECONGATE_TOP;
    tbot=DECONGATE_BOT;
    %set the design gate
    hseis1=findobj(gcf,'tag','seis1');
    htop=findobj(hseis1,'tag','ttop');
    set(htop,'ydata',ttop*ones(1,2));
    hbot=findobj(hseis1,'tag','tbot');
    set(hbot,'ydata',tbot*ones(1,2));
elseif(strcmp(action,'close'))
    hspec=findobj(gcf,'tag','spectra');
    hspecwin=get(hspec,'userdata');
    delete(hspecwin);
    crf=get(gcf,'closerequestfcn');
    ind=strfind(crf,';');
    if(ind(1)==length(crf))
        delete(gcf);
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

%clips=linspace(nsigma,1,nclips)';
clips=[20 15 10 8 6 4 3 2 1 .5 .25 .1 .075 .05 .025 .01 .005 .001 .0001]';
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