function seispick(seis,t,x,dname)
% SEISPLOT ... create a sesimic image plot in a new figure window
%
% seisplot(seis,t,x,dname)
%
% This function provides a similar display to that of plotimage but without the extra features (such
% as picking). This is suitable when you want a quick seismic display in a new figure but don't nead
% the extra features of plotimage. Seisplot also gives interactive ability to choose the clip level.
% Use seisplota for a similar display in the current axes. You might also consider seisplotfk and
% seisplottwo. The former shows a seismic section and its f-k spectrum while the latter shows two
% seismic sections in a single window for comparison.
%
% seis ... seismic matrix (gather).One trace per column
% t ... time coordinate vector. length(t) must equal size(seis,1)
% x ... space coordinate vector. length(x) must equal size(seis,2)
% dname ... dataset name (string). Used to title the plot and to label the figure
% ************** default = [] ****************
%
% NOTE: The image is plotted with Matlab's imagesc. This function only annotates the axes
% precisely correctly if the x and t vectors are regularly sampled. This is usually the case
% with t but less often so with x. For precisly annotated tick marks when x is not regular, the
% only current option is to uses plotseis or plotseismic which both plot wiggle traces not
% images.
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

global ZOOM_VALUE

if(~ischar(seis))
    action='init';
else
    action=seis;
end

if(strcmp(action,'init'))
    
    if(nargin<4)
        dname=[];
    end
    xname='';
    yname='';
    if(nargin<3)
        x=1:size(seis,2);
        xname='column number';
    end
    if(nargin<2)
        t=(1:size(seis,1))';
        yname='row number';
    end
    
    if(length(t)~=size(seis,1))
        error('time coordinate vector does not match seismic');
    end
    if(length(x)~=size(seis,2))
        error('space coordinate vector does not match seismic');
    end
    
    figure
    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis);
    clim=[am-clip*sigma am+clip*sigma];
        
    hi=imagesc(x,t,seis,clim);colormap(seisclrs);
%     brighten(.5);
    grid
    hcm=uicontextmenu;
    uimenu(hcm,'label','Time-variant spectrum','callback',@showtvspectrum);
    uimenu(hcm,'label','2D spectrum','callback',@showfkspectrum);
    set(hi,'uicontextmenu',hcm);
    
    if(~isempty(dname))
        ht=title(dname);
        set(ht,'interpreter','none');
        set(gcf,'name',dname);
    end
    maxmeters=7000;
    if(isempty(yname))
        if(max(t)<10)
            yname='time (s)';
        elseif(max(t)<maxmeters)
            yname='depth (m)';
        else
            yname='(depth (ft)';
        end
    end
    ylabel(yname);
    if(isempty(xname))
        if(max(x)<maxmeters)
            xname='distance (m)';
        else
            xname='distance (ft)';
        end
    end
    xlabel(xname);
    
    %compute the average Hilbert envelope
    aveenv=zeros(size(seis,1),1);
    for k=1:length(x)
        aveenv=aveenv+env(double(seis(:,k)));
    end
    aveenv=aveenv/length(x);
    
    %make a clip control
    xnow=.05;ynow=.875;wid=.055;ht=.05;sep=.005;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''clip'')','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,dname,aveenv,yname,xname},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped');
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brighten','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''brighten'')',...
        'tooltipstring','push once or multiple times to brighten the image');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darken','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''brighten'')',...
        'tooltipstring','push once or multiple times to darken the image');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightness','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','current image brightness','userdata',0);
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Amp histogram','tag','histogram','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''histogram'')',...
        'tooltipstring','show amplitude histogram');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Publish zoom','tag','histogram','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''publish'')',...
        'tooltipstring','Publish current zoom limits');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','Match zoom','tag','histogram','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''match'')',...
        'tooltipstring','Match zoom to published limits');
    ynow=ynow-ht-sep;
    xl=get(gca,'xlim');
    yl=get(gca,'ylim');
    uicontrol(gcf,'style','pushbutton','string','Un-zoom','tag','histogram','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplot(''unzoom'')',...
        'tooltipstring','Un-zoom to original view','userdata',{xl yl});
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
elseif(strcmp(action,'clip'))
    hclip=findobj(gcf,'tag','clip');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
    sigma=udat{3};
    if(iclip==1)
        clim=[-amax amax];
    else
        clip=clips(iclip-1);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set(gca,'clim',clim);
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
elseif(strcmp(action,'histogram'))
    p=get(gcf,'position');
    hi=findobj(gca,'type','image');
    seis=get(hi,'cdata');
    t=get(hi,'ydata');
    nsamps=numel(seis);
    hclip=findobj(gcf,'tag','clip');
    udat=get(hclip,'userdata');
    am=udat{2};
    amax=udat{4};
    amin=udat{5};
    dname=udat{6};
    aveenv=udat{7};
    yname=udat{8};
    sigma=udat{3};
    figure
    q3=.5*p(3);
    q4=.5*p(4);
    q1=p(1)+.5*(p(3)-q3);
    q2=p(2)+.5*(p(4)-q4);
    set(gcf,'position',[q1 q2 q3 q4],'name',['Amplitude histogram ' dname])
    nbins=nsamps/500;
    if(nbins<100); nbins=100; end
    if(nbins>1000);nbins=1000;end
    subplot(2,1,1)
    hist(seis(:),nbins);
    ht=title({dname 'Amplitude histogram'});
    set(ht,'interpreter','none');
    xlabel('Amplitude');
    ylabel('Number of samples');
    xl=get(gca,'xlim');
    yl=get(gca,'ylim');
    yinc=diff(yl)/10;
    xinc=diff(xl)/20;
    fs=9;
    ynow=yl(2)-yinc;
    text(xl(1)+xinc,ynow,['Max amp= ' num2str(amax)],'fontsize',fs)
    ynow=ynow-yinc;
    text(xl(1)+xinc,ynow,['Min amp= ' num2str(amin)],'fontsize',fs)
    ynow=ynow-yinc;
    text(xl(1)+xinc,ynow,['Mean amp= ' num2str(am)],'fontsize',fs)
    ynow=ynow-yinc;
    text(xl(1)+xinc,ynow,['Std dev= ' num2str(sigma)],'fontsize',fs)
    subplot(2,1,2)
    %fit exponential
    nt2=round(length(t)/1);
    ind=near(t,t(1),t(nt2));
    Emax=max(aveenv);
    p=polyfit(t(ind),log(aveenv(ind)+.0001*Emax),2);
    expfit=exp(polyval(p,t(ind)));
    plot(t,aveenv,t(ind),expfit(ind),'r')
    legend('average envelope',['log(env)= ' num2str(sigfig(p(1),2)) 't^2 + ' num2str(sigfig(p(2),2)) 't + ' num2str(sigfig(p(3),2))])
    title('Average trace envelope')
    xlabel(yname)
    ylabel('Amplitude')
elseif(strcmp(action,'publish'))
    yl=get(gca,'ylim');
    xl=get(gca,'xlim');
    ZOOM_VALUE{1}=xl;
    ZOOM_VALUE{2}=yl;
elseif(strcmp(action,'match'))
    if(isempty(ZOOM_VALUE))
        return;
    end
    xl=ZOOM_VALUE{1};
    yl=ZOOM_VALUE{2};
    hi=findobj(gca,'type','image');
    x=get(hi,'xdata');
    y=get(hi,'ydata');
    fudge=diff(xl)*.1;
    x1=min(x)-fudge;
    x2=max(x)+fudge;
    if(~between(x1,x2,xl(1),2) && ~between(x1,x2,xl(2),2))
        msgbox('Published zoom limits are incompatible with this data');
        return;
    end
    fudge=diff(yl)*.1;
    y1=min(y)-fudge;
    y2=max(y)+fudge;
    if(~between(y1,y2,yl(1),2) && ~between(y1,y2,yl(2),2))
        msgbox('Published zoom limits are incompatible with this data');
        return;
    end
    set(gca,'xlim',xl,'ylim',yl)
elseif(strcmp(action,'unzoom'))
    udat=get(gco,'userdata');
    xl=udat{1};
    yl=udat{2};
    set(gca,'xlim',xl,'ylim',yl);
    
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

function showtvspectrum(~,~)
%get the data
hi=findobj(gca,'type','image');
%x=get(hi,'xdata');
t=get(hi,'ydata');
seis=get(hi,'cdata');

tmin=min(t);
tmax=max(t);

dname=get(gcf,'name');

%divide into 3 zones
tdel=(tmax-tmin)/3;
tnots=[tmin tmin+tdel tmin+2*tdel];
twins=tdel*ones(size(tdel));
tpad=2*tdel;

pos=get(gcf,'position');
figure
tvdbspec(t,seis,tnots,twins,tpad,dname,gca);
prepfiga
set(gcf,'position',pos);

end

function showfkspectrum(~,~)
%get the data
hi=findobj(gca,'type','image');
x=get(hi,'xdata');
t=get(hi,'ydata');
seis=get(hi,'cdata');

dname=get(gcf,'name');

time=1;
if(max(t)>30)
    time=0;
end

[fk,f,k]=fktran(seis,t,x);
seisplot(abs(fk),f,k,[dname ' 2D spectrum'])
if(time==1)
    xlabel('Wavenumber');ylabel('Frequency (Hz)');
else
    xlabel('Wavenumber');ylabel('Wavenumber');
end

end