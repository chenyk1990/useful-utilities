function datar=seisplotsvd_sep(seis,t,x,dname)
% seisplotsvd_sep: plots a seismic gather and its SVD separation into Gross and Detail side-by-side
%
% datar=seisplotsvd_sep(seis,t,x,dname)
%
% A new figure is created and divided into three same-sized axes (side-by-side). The seismic gather
% (matrix) is plotted as an image in the left-hand-side and its SVD separation into Gross and Detail
% (see svd_sep) cousin are plotted as images other two axes. Controls are provided to adjust the
% singular value cuttoff and to brighten or darken the image plots.
%
% seis ... input seismic matrix
% t ... time coordinate vector for seis (y or row coordinate). Only used
%       for plotting
% *********** default 1:nrows ***************
% x ... space coordinate vector for seis (x or column coordinate)
% *********** default 1:ncols ***************
% dname ... text string giving a name for the dataset that will annotate
%       the plots.
% ************ default dname =[] ************
%
% datar ... Return data which is a length 4 cell array containing
%           data{1} ... handle of the seismic axes
%           data{2} ... handle of the svd Gross axes
%           data{3} ... handle of the svd Detail axes
%           data{4} ... singular values (vector)
% These return data are provided to simplify plotting additional lines and
% text in either axes.
%
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
    
    [nrows,ncols]=size(seis);
    if(nargin<3)
        x=1:ncols;
    end
    if(nargin<2)
        t=(1:nrows)';
    end
    if(length(t)~=nrows)
        error('time coordinate vector does not match seismic');
    end
    if(length(x)~=ncols)
        error('space coordinate vector does not match seismic');
    end
    
    if(nargin<4)
        dname=[];
    end
    
    xwid=.25;
    yht=.75;
    xsep=.05;
    xnot=.05;
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
    
    %process dname if it is too long
    toolong=50;
    if(iscell(dname))
        if(length(dname{1})>toolong)
            str1=dname{1};
            str2=dname{2};
            ind=isspace(str1);
            ind2=find(ind>0);%points to word breaks
            ind3=find(ind2<toolong);
            if(~isempty(ind3))
               str1a=str1(1:ind2(ind3(end)));
               str2a=[str1(ind2(ind3(end))+1:end) ' ' str2];
               dname{1}=str1a;
               dname{2}=str2a;
            end
        end
    end
            
    ht=title(dname);
    ht.Interpreter='none';
    
    maxmeters=7000;
    
    if(max(t)<10)
        ylabel('time (s)')
    elseif(max(t)<maxmeters)
        ylabel('distance (m)')
    else
        ylabel('distance (ft)')
    end
    if(max(x)<maxmeters)
        xlabel('distance (m)')
    else
        xlabel('distance (ft)')
    end
    %make a clip control

    xnow=xnot+xwid;
    wid=.04;ht=.05;sep=.005;
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clipxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_sep(''clipxt'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brightenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_sep(''brightenxt'');',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darkenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_sep(''brightenxt'');',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightnessxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    
    set(hax1,'tag','seis');
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid yht]);
    hax3=subplot('position',[xnot+2*xwid+1.5*xsep ynot xwid yht]);
    [U,S,V]=svd(seis);
    singvals=diag(S);
    
    %make the singcut axes
    hax4=subplot('position',[xnot+3*xwid+1.55*xsep ynot .25*xwid .8*yht]);
    
    nsing=length(singvals);
    hh=semilogy(singvals,1:nsing);flipy
    set(hh,'tag','singvals');
    %ylabel('singular values')
    set(hax4,'xtick',[],'yaxislocation','right','nextplot','add')
    xlim([floor(singvals(end)) ceil(singvals(1))])
    ylim([1 length(singvals)])
    singcut=10;
    line([singvals(1) singvals(end)],singcut*ones(1,2),'color','r',...
        'buttondownfcn','seisplotsvd_sep(''dragline'');','tag','singcut');
    set(hax4,'tag','singcutaxe','userdata',{U,singvals,V,singcut,dname});
    if(iscell(dname))
        set(gcf,'name',['SVD analysis for ' dname{1}]);
    else
        set(gcf,'name',['SVD analysis for ' dname]);
    end
    
    [gross,detail,singvalsg,singvalsd]=reconstruct(singcut);
    
    hhg=semilogy(singvalsg,1:nsing,'r.');
    set(hhg,'tag','singvalsg');
    hhd=semilogy(singvalsd,1:nsing,'g.');
    set(hhd,'tag','singvalsd');
    legend([hh hhg hhd],'All','Gross','Detail','location','southeast');
    title('SingVals')
    
    axes(hax2);

    %make a clip control
    xnow=xnot+3*xwid+1.5*xsep;
    ht=.05;
    ynow=ynot+yht-ht;
    %wid=.045;sep=.005;
    hclipsvd=uicontrol(gcf,'style','popupmenu','string','xxx','tag','clipsvd','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_sep(''clipsvd'');','value',1,...
        'userdata',[],'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped');
    if(max(t)>10)
        uicontrol(gcf,'style','pushbutton','string','Grid CC','tag','gridcc','units','normalized',...
            'position',[xnow+wid,ynow+.5*ht,wid,.5*ht],'callback','seisplotsvd_sep(''gridcc'');',...
            'userdata',[],'tooltipstring','Show correlations with the acquisition grid.');
    end
    ynow2=ynow-ht-sep;
    hbutgrp=uibuttongroup(gcf,'units','normalized','position',[xnow,ynow2,1.5*wid,1.5*ht],'tag','cutoff');
%     uicontrol(hbutgrp,'style','radio','string','Gauss cutoff','units','normalized','position',[0,.5,1,.5],...
%         'enable','on','tag','gauss','backgroundcolor','w');
%     uicontrol(hbutgrp,'style','radio','string','Sharp cutoff','units','normalized','position',[0,0,1,.5],...
%         'enable','on','tag','sharp','backgroundcolor','w');
    uicontrol(hbutgrp,'style','radio','string','Gauss cutoff','units','normalized','position',[0,.5,1,.5],...
        'enable','on','tag','gauss','backgroundcolor','w','callback','seisplotsvd_sep(''singcutold'');');
    uicontrol(hbutgrp,'style','radio','string','Sharp cutoff','units','normalized','position',[0,0,1,.5],...
        'enable','on','tag','sharp','backgroundcolor','w','callback','seisplotsvd_sep(''singcutold'');');
    
    
    %make an info button
    msg=['The seismic matrix on the left is shown separated into its Gross structure and its Detail. ',...
        'This separation is done using SVD (singular-value decomposition) and is controlled by the parameter singcut. ',...
        'There are typically hundreds to thousands of singular values in an image where, sorted from largest to smallest, ',...
        'the first few control the Gross structure and the remainder control the detail. These are plotted in the ',...
        'axis on the right. You can adjust the cutoff singular value (singcut) by clicking and dragging the horizontal ',...
        'red line in the SingVals axis. When you release the line a new separation will be displayed.',...
        'The singular values defining the Gross structure are the original singular values multiplied by a Gaussian ',...
        'centered on the largest singular value and whose standard deviation is singcut. The Detail singular values are the original singular values ',...
        'minus those that define the Gross structure. In this way the original seismic matrix is always equal ',...
        'to the sum of Gross and Detail.'];
    uicontrol(gcf,'style','pushbutton','string','Info','tag','info','units','normalized',...
        'userdata',msg,'position',[xnow,ynow+2*ht,wid,ht],'callback','seisplotsvd_sep(''info'');',...
        'backgroundcolor','y','fontsize',16);


    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(gross);
    set(hclipsvd,'userdata',{clips,am,sigma,amax,amin,[hax2 hax3]});
    set(hclipsvd,'string',clipstr,'value',iclip);
    
    %clim=[amin am+clip*sigma];
    if(iclip==1)
        clim=[-amax amax];
    else
        clim=[am-clip*sigma am+clip*sigma];
    end
        
    imagesc(x,t,gross,clim);colormap(seisclrs);
%     brighten(.5);
    grid
    
    title({'SVD gross',['singcut= ' int2str(singcut)]});
    
    set(gca,'yticklabel','');
    if(max(x)<maxmeters)
        xlabel('distance (m)')
    else
        xlabel('distance (ft)')
    end
    axes(hax3);
    imagesc(x,t,detail,clim);colormap(seisclrs);
%     brighten(.5);
    grid
    
    title({'SVD detail',['singcut= ' int2str(singcut)]});
    
    set(gca,'yticklabel','');
    if(max(x)<maxmeters)
        xlabel('distance (m)')
    else
        xlabel('distance (ft)')
    end
    
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
    set(hax2,'tag','svd1');
    set(hax3,'tag','svd2');
    
    if(nargout>0)
        datar=cell(1,4);
        datar{1}=hax1;
        datar{2}=hax2;
        datar{3}=hax3;
        datar{4}=singvals;
    end
elseif(strcmp(action,'dragline'))
    %hthisline=gcbo;
    %h1=findobj(gcf,'tag','singcut');
    hax4=findobj(gcf,'tag','singcutaxe');
    yl=get(hax4,'ylim');
    nsing=yl(2);
    xl=get(hax4,'xlim');
    %tmp=get(h1,'ydata');
    %singcut=round(tmp(1));
    DRAGLINE_MOTION='yonly';
    DRAGLINE_XLIMS=xl;
    DRAGLINE_YLIMS=[1 nsing];
    DRAGLINE_SHOWPOSN='on';
    DRAGLINE_CALLBACK='seisplotsvd_sep(''singcutold'');';
    DRAGLINE_MOTIONCALLBACK='';
    DRAGLINE_PAIRED=[];
    dragline('click');
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
elseif(strcmp(action,'clipsvd'))
    hclip=findobj(gcf,'tag','clipsvd');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
    %amin=udat{5};
    sigma=udat{3};
    haxes=udat{6};
    if(iclip==1)
        %clim=[amin amax];
        clim=[-amax amax];
    else
        clip=clips(iclip-1);
        clim=[am-clip*sigma,am+clip*sigma];
        %clim=[amin am+clip*sigma];
    end
    set(haxes(1),'clim',clim);
    set(haxes(2),'clim',clim);
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
elseif(strcmp(action,'brightensvd'))
    hbut=gcbo;
    hbright=findobj(gcf,'tag','brightensvd');
    if(hbut==hbright)
        inc=.1;
    else
        inc=-.1;
    end
    brighten(inc);
    hbrightness=findobj(gcf,'tag','brightnesssvd');
    brightlvl=get(hbrightness,'userdata');
    brightlvl=brightlvl+inc;
    if(abs(brightlvl)<.01)
        brightlvl=0;
    end
    set(hbrightness,'string',['lvl ' num2str(brightlvl)],'userdata',brightlvl)
elseif(strcmp(action,'singcutold'))
    hsingcutaxe=findobj(gcf,'tag','singcutaxe');
    udat=get(hsingcutaxe,'userdata');
    singvals=udat{2};
    nsing=length(singvals);
    h1=findobj(gcf,'tag','singcut');
    val=get(h1,'ydata');
    singcut=round(val(1));
    if(singcut<1); singcut=1; end
    if(singcut>nsing); singcut=nsing; end
    %dname=udat{5};
    
    [gross,detail,singvalsg,singvalsd]=reconstruct(singcut);
    haxesvd1=findobj(gcf,'tag','svd1');
    axes(haxesvd1);
    hi=findobj(haxesvd1,'type','image');
    set(hi,'cdata',gross);
    title({'SVD Gross',['singcut= ' int2str(singcut)]});
    haxesvd2=findobj(gcf,'tag','svd2');
    axes(haxesvd2);
    hi=findobj(haxesvd2,'type','image');
    set(hi,'cdata',detail);
    hhg=findobj(gcf,'tag','singvalsg');
    set(hhg,'xdata',singvalsg);
    hhg=findobj(gcf,'tag','singvalsd');
    set(hhg,'xdata',singvalsd);
    title({'SVD Detail',['singcut= ' int2str(singcut)]});
elseif(strcmp(action,'gridcc')||strcmp(action,'gridcc2'))
    hcc=findobj(gcf,'tag','gridcc');
    if(strcmp(action,'gridcc2'))
       a=askthingsfini;
       if(a==-1)
           return;
       end
       dx=str2double(a(1,:));
       if(isnan(dx))
           msgbox('Bad numerica value!! Try again');
           return
       end
       dy=str2double(a(2,:));
       if(isnan(dy))
           msgbox('Bad numerica value!! Try again');
           return
       end
       dxaq=str2double(a(3,:));
       if(isnan(dxaq))
           msgbox('Bad numerica value!! Try again');
           return
       end
       dyaq=str2double(a(4,:));
       if(isnan(dyaq))
           msgbox('Bad numerica value!! Try again');
           return
       end
       set(hcc,'userdata',[dx dy dxaq dyaq]);
    end
    %get the grid info if it exists
    udat=get(hcc,'userdata');
    if(isempty(udat))
       %put up the askthings dialog
       %the transfer
       transfer='seisplotsvd_sep(''gridcc2'');';
       %the questions
       q1='Physical distance between xlines';
       q2='Physical distance between inlines';
       q3='Acquisition line spacing in xline direction';
       q4='Acquisition line spacing in inline direction';
       q=char(q1,q2,q3,q4);
       askthingsinit(transfer,q);
       return;
    end
    hthisfig=gcf;
    dx=udat(1);
    dy=udat(2);
    dxaq=udat(3);
    dyaq=udat(4);
    %now get each image and compute its correlations
    %the slice
    hax=findobj(hthisfig,'tag','seis');
    hi=findobj(hax,'type','image');
    seis=get(hi,'cdata');
    x=get(hi,'xdata');
    y=get(hi,'ydata');
    [ccx1,xlags1,ccy1,ylags1,xaqgrid,yaqgrid]=ccfoot(seis,x*dx,y*dy,dxaq,dyaq);
    %the gross
    hax=findobj(hthisfig,'tag','svd1');
    hi=findobj(hax,'type','image');
    seis=get(hi,'cdata');
    x=get(hi,'xdata');
    y=get(hi,'ydata');
    [ccx2,xlags2,ccy2,ylags2]=ccfoot(seis,x*dx,y*dy,dxaq,dyaq);
    %the detail
    hax=findobj(hthisfig,'tag','svd2');
    hi=findobj(hax,'type','image');
    seis=get(hi,'cdata');
    x=get(hi,'xdata');
    y=get(hi,'ydata');
    [ccx3,xlags3,ccy3,ylags3]=ccfoot(seis,x*dx,y*dy,dxaq,dyaq);
    %now put up a new figure and plot
    pos=get(hthisfig,'position');
    hccfig=figure;
    subplot(1,3,1)
    plot(xlags1,ccx1,ylags1,ccy1);
    y0=min([ccx1(:);ccy1(:)]);
    y0=max([0, floor(y0/.1)*.1-.1]);
    y1=max([ccx1(:);ccy1(:)]);
    y1=min([1, ceil(y1/.1)*.1+.1]);
    ylim([y0 y1]);
    xlabel('lag');
    ylabel('correlation')
    title('Full seismic');
    legend('x grid correlation','y grid correlation','location','south');
    %legend('source grid','receiver grid','location','south');
    xl=get(gca,'xlim');
    delx=xl(2)/2;
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    subplot(1,3,2)
    plot(xlags2,ccx2,ylags2,ccy2);
    ylim([y0 y1]);
    xlabel('lag');
    ylabel('correlation')
    title('Gross');
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    subplot(1,3,3)
    plot(xlags3,ccx3,ylags3,ccy3);
    ylim([y0 y1]);
    xlabel('lag');
    ylabel('correlation')
    title('Detail');
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    boldlines;
    bigfont(gcf,1.5,1);
    wid=pos(3)*.5;
    ht=pos(4)*.5;
    hax4=findobj(hthisfig,'tag','singcutaxe');
    udat=get(hax4,'userdata');
    dname=udat{5};
    if(iscell(dname))
        dname= dname{1};
    end
    set(hccfig,'position',[pos(1)+.5*(pos(3)-wid),pos(2)+.5*(pos(4)-ht),wid,ht],'name',['Grid correlations for ' dname])
    
    %showgridsbutton
    uicontrol(hccfig,'style','pushbutton','string','Show grids','units','normalized',...
        'position',[.91,.85,.08,.05],'callback','seisplotsvd_sep(''showgrids'');',...
        'userdata',{xaqgrid yaqgrid x y dname},'tag','showgrids');
    
%     figure
%     plot(xlags1,ccx1,ylags1,ccy1);
%     y0=min([ccx1(:);ccy1(:)]);
%     y0=max([0, floor(y0/.1)*.1-.1]);
%     y1=max([ccx1(:);ccy1(:)]);
%     y1=min([1, ceil(y1/.1)*.1+.1]);
%     ylim([y0 y1]);
%     xlabel('lag');
%     ylabel('correlation')
%     title('Full seismic');
%     grid
%     boldlines;
%     bigfont(gcf,1.5,1);
%     legend('source grid','receiver grid','location','south');
%     pos=get(gcf,'position');
%     set(gcf,'position',[200 200 600 800]);
    
elseif(strcmp(action,'showgrids'))
    hshowgrids=findobj(gcf,'tag','showgrids');
    udat=get(hshowgrids,'userdata');
    xaqgrid=udat{1};
    yaqgrid=udat{2};
    x=udat{3};
    y=udat{4};
    dname=udat{5};
    datar=seisplottwo(xaqgrid,y,x(1:end-1),'X grid',yaqgrid,y(1:end-1),x,'Y grid');
    axes(datar{1}); xlabel('');ylabel('');
    axes(datar{2}); xlabel('');ylabel('');
    pos=get(gcf,'position');
    set(gcf,'position',[pos(1:2) .67*pos(3) pos(4)],'name',['Aquisition grids for ' dname]);
    
    
elseif(strcmp(action,'info'))
    hinfo=findobj(gcf,'tag','info');
    udat=get(hinfo,'userdata');
    if(iscell(udat))
        msg=udat{1};
        h=udat{2};
        if(isgraphics(h))
            delete(h);
        end
    else
        msg=udat;
    end
    h=msgbox(msg,'SVD Separation','help');
    set(hinfo,'userdata',{msg,h});
end
end

function [gross,detail,singvalsg,singvalsd]=reconstruct(singcut)
hsingcutaxe=findobj(gcf,'tag','singcutaxe');
udat=get(hsingcutaxe,'userdata');
hgauss=findobj(gcf,'tag','gauss');
if(~isempty(hgauss))
    flag=get(hgauss,'value');
else
    flag=1;
end
U=udat{1};
singvals=udat{2};
V=udat{3};
m=size(U,1);
n=size(V,1);
nsing=length(singvals);
j=1:nsing;
if(flag==1)
    g=exp(-(j-1).^2/singcut^2)';
else
    g=zeros(size(j))';
    ind= j<=singcut;
    g(ind)=1;
end
h=1-g;
singvalsg=singvals.*g;
tmp=diag(singvalsg);
if(m>n)
    Sg=[tmp;zeros(m-n,n)];
elseif(n>m)
    Sg=[tmp zeros(m,n-m)];
else
    Sg=tmp;
end
singvalsd=singvals.*h;
tmp=diag(singvalsd);
if(m>n)
    Sd=[tmp;zeros(m-n,n)];
elseif(n>m)
    Sd=[tmp zeros(m,n-m)];
else
    Sd=tmp;
end
gross=U*Sg*V';
detail=U*Sd*V';

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