function datar=seisplotsvd_foot(slice,x,y,dx,dy,dname)
% seisplotsvd_foot: plots a seismic gather and its SVD separation into Gross and Detail side-by-side
%
% datar=seisplotsvd_foot(slice,x,y,dx,dy,dname)
%
% A new figure is created and divided into three same-sized axes (side-by-side). The seismic gather
% (matrix) is plotted as an image in the left-hand-side and its SVD separation into Gross and Detail
% (see svd_sep) cousin are plotted as images other two axes. Controls are provided to adjust the
% singular value cuttoff and to brighten or darken the image plots.
%
% slice ... input seismic matrix
% x ... column coordinate vector for slice. This is usually crossline number.
% y ... row coordinate vector for slice. This will typically be inline number.
% dx ... physical distance (ft or meters) between crosslines.
% dy ... physical distance (ft or meters) between inlines.
% dname ... text string giving a name for the dataset that will annotate
%       the plots.
% ************ default dname =[] ************
%
% datar ... Return data which is a length 6 cell array containing
%           data{1} ... handle of the seismic axes
%           data{2} ... handle of the svd Gross axes
%           data{3} ... handle of the svd Detail axes
%           data{4} ... handle of the result axes
%           data{5} ... handle of the filtered Gross axes
%           data{6} ... handle of the difference axes
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

if(~ischar(slice))
    action='init';
else
    action=slice;
end

datar=[];%initialize return data to null

if(strcmp(action,'init'))
    
    [nrows,ncols]=size(slice);
    if(nargin<3)
        x=1:ncols;
    end
    if(nargin<2)
        y=(1:nrows)';
    end
    if(length(y)~=nrows)
        error('y coordinate vector does not match seismic');
    end
    if(length(x)~=ncols)
        error('space coordinate vector does not match seismic');
    end
    
    if(nargin<4)
        dname=[];
    end
    
    xwid=.25;
    yht=.35;
    xsep=.05;
    xnot=.05;
    ynot=.55;
    

    figure
    hax1=subplot('position',[xnot ynot xwid yht]);

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(slice);
    if(iclip==1)
        clim=[-amax amax];
    else
        clim=[am-clip*sigma am+clip*sigma];
    end
        
    imagesc(x,y,slice,clim);colormap(seisclrs)
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
        dname{1}=['Input: ' dname{1}];
    else
        dname=['Input: ' dname];
    end        
    ht=title(dname);
    ht.Interpreter='none';
    ylabel('inline number');
    set(gca,'xticklabel','');
    
    %make a clip control

    xnow=xnot+xwid;
    wid=.04;ht=.05;sep=.005;
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clipxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_foot(''clipxt'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brightenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_foot(''brightenxt'');',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darkenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_foot(''brightenxt'');',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightnessxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid yht]);
    hax3=subplot('position',[xnot+2*xwid+1.5*xsep ynot xwid yht]);
    [U,S,V]=svd(slice);
    singvals=diag(S);
    
    %make the singcut axes
    ysing=.1;
    hax4=subplot('position',[xnot+3*xwid+1.55*xsep ysing .25*xwid 1.5*yht]);
    
    nsing=length(singvals);
    hh=semilogy(singvals,1:nsing);flipy
    set(hh,'tag','singvals');
    %ylabel('singular values')
    set(hax4,'xtick',[],'yaxislocation','right','nextplot','add')
    xlim([floor(singvals(end)) ceil(singvals(1))])
    ylim([1 length(singvals)])
    singcut=10;
    line([singvals(1) singvals(end)],singcut*ones(1,2),'color','r',...
        'buttondownfcn','seisplotsvd_foot(''dragline'');','tag','singcut');
    set(hax4,'tag','singcutaxe','userdata',{U,singvals,V,singcut,dname});
    if(iscell(dname))
        set(gcf,'name',['Footprint analysis for ' dname{1}]);
    else
        set(gcf,'name',['Footprint analysis for ' dname]);
    end
    
    [gross,detail,singvalsg,singvalsd]=decompose(singcut);
    
    hhg=semilogy(singvalsg,1:nsing,'r.');
    set(hhg,'tag','singvalsg');
    hhd=semilogy(singvalsd,1:nsing,'g.');
    set(hhd,'tag','singvalsd');
    legend([hh hhg hhd],'All','Gross','Detail','location','southeast');
    title('SingVals')
    legendfontsize(.9)
    
    axes(hax2);

    %make a clip control
    xnow=xnot+3*xwid+1.5*xsep;
    ht=.05;
    ynow=ynot+yht-ht;
    %wid=.045;sep=.005;
    hclipsvd=uicontrol(gcf,'style','popupmenu','string','xxx','tag','clipsvd','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_foot(''clipsvd'');','value',1,...
        'userdata',[],'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped');
    uicontrol(gcf,'style','pushbutton','string','Grid CC','tag','gridcc','units','normalized',...
        'position',[xnow+wid,ynow+.5*ht,wid,.5*ht],'callback','seisplotsvd_foot(''gridcc'');',...
        'userdata',[],'tooltipstring','Show correlations with the acquisition grid.');
    ynow2=ynow-.25*ht;
    hbutgrp=uibuttongroup(gcf,'units','normalized','position',[xnow,ynow2,1.5*wid,.75*ht],'tag','cutoff');
    uicontrol(hbutgrp,'style','radio','string','Gauss cutoff','units','normalized','position',[0,.5,1,.5],...
        'enable','on','tag','gauss','backgroundcolor','w');
    uicontrol(hbutgrp,'style','radio','string','Sharp cutoff','units','normalized','position',[0,0,1,.5],...
        'enable','on','tag','sharp','backgroundcolor','w');
    
    %wavenumber controls
    sigmaxvals={'1.0','0.75' '0.5' '0.375' '0.25' '0.125' '0.0625' '0.0313' '0.0156' '0.0078','0.0039','0.0020'};
    sigmayvals=sigmaxvals;
    isigG=6;
    ynow2=ynow2-ht;
    nudge=.1*ht;
    uicontrol(gcf,'style','text','string','Sigmax(G):','units','normalized','position',[xnow,ynow2-nudge,.75*wid,ht],...
        'tooltipstring','Gross filtering: the smaller the number the greater the filtering','backgroundcolor','w');
    uicontrol(gcf,'style','popupmenu','string',sigmaxvals,'units','normalized','position',[xnow+.75*wid,ynow2,.75*wid,ht],...
        'tooltipstring','Values are fractions of crossline wavenumber Nyquist','value',isigG,'tag','sigmaxG',...
        'callback','seisplotsvd_foot(''sigmachange'');','userdata',{dx dy});
    ynow2=ynow2-.5*ht;
    uicontrol(gcf,'style','text','string','Sigmay(G):','units','normalized','position',[xnow,ynow2-nudge,.75*wid,ht],...
        'tooltipstring','Gross filtering: the smaller the number the greater the filtering','backgroundcolor','w');
    hsigmay=uicontrol(gcf,'style','popupmenu','string',sigmayvals,'units','normalized','position',[xnow+.75*wid,ynow2,.75*wid,ht],...
        'tooltipstring','Values are fractions of inline wavenumber Nyquist','value',isigG,'tag','sigmayG',...
        'callback','seisplotsvd_foot(''sigmachange'');');
    isigD=1;
    ynow2=ynow2-.5*ht;
    nudge=.1*ht;
    uicontrol(gcf,'style','text','string','Sigmax(D):','units','normalized','position',[xnow,ynow2-nudge,.75*wid,ht],...
        'tooltipstring','Detail filtering: the smaller the number the greater the filtering','backgroundcolor','w');
    uicontrol(gcf,'style','popupmenu','string',sigmaxvals,'units','normalized','position',[xnow+.75*wid,ynow2,.75*wid,ht],...
        'tooltipstring','Values are fractions of crossline wavenumber Nyquist','value',isigD,'tag','sigmaxD',...
        'callback','seisplotsvd_foot(''sigmachange'');');
    ynow2=ynow2-.5*ht;
    uicontrol(gcf,'style','text','string','Sigmay(D):','units','normalized','position',[xnow,ynow2-nudge,.75*wid,ht],...
        'tooltipstring','Detail filtering: the smaller the number the greater the filtering','backgroundcolor','w');
    uicontrol(gcf,'style','popupmenu','string',sigmayvals,'units','normalized','position',[xnow+.75*wid,ynow2,.75*wid,ht],...
        'tooltipstring','Values are fractions of inline wavenumber Nyquist','value',isigD,'tag','sigmayD',...
        'callback','seisplotsvd_foot(''sigmachange'');');
    %zoom buttons
    ynow2=ynow2-.5*ht;
    uicontrol(gcf,'style','pushbutton','string','Equal zoom','units','normalized','position',[xnow,ynow2,wid,.7*ht],...
        'callback','seisplotsvd_foot(''eqzoom'');','tag','eqzoom',...
        'tooltipstring','Zoom one axis then push this to equalize all axes');
    
    uicontrol(gcf,'style','pushbutton','string','Unzoom all','units','normalized','position',[xnow+wid,ynow2,wid,.7*ht],...
        'callback','seisplotsvd_foot(''unzoom'');','tag','unzoom',...
        'tooltipstring','Restore all axes to original ranges');
    
    
    %make an info button
    msg=['The seismic matrix on the top-left is shown separated into its Gross structure (top center) and its Detail (top right). ',...
        'This separation is done using SVD (singular-value decomposition) and is controlled by the parameter singcut. ',...
        'There are typically hundreds to thousands of singular values in an image where, sorted from largest to smallest, ',...
        'the first few control the Gross structure and the remainder control the detail. These are plotted in the ',...
        'axis on the right. You can adjust the cutoff singular value (singcut) by clicking and dragging the horizontal ',...
        'red line in the SingVals axis. When you release the line a new separation will be displayed.',...
        'The singular values defining the Gross structure are the original singular values multiplied by a Gaussian ',...
        'centered on the largest singular value and whose standard deviation is singcut. The Detail singular values are the original singular values ',...
        'minus those that define the Gross structure. In this way the original seismic matrix is always equal ',...
        'to the sum of Gross and Detail. On the bottom row are the results of using this separation to suppress footprint or other problems.',...
        'The idea is to find a separation where the footprint is mostly in the Gross panel. Then both panels are wavenumber filtered ',...
        'to suppress higher wavenumbers but the Gross panel is usually more harshly filtered than the Detail. ',...
        'Filtering is done by multiplying the wavenumber spectrum by a Gaussian whose standard ',...
        'deviations in x and y are specified by the Sigmax(G), Sigmay(G), Sigmax(D), Sigmay(D) ',...
        'controls. These sigma values are specified as fractions of the corresponding spatial Nyquist.',...
        ' Thus sigmax=0.25 means one-quarter of the Nyquist wavenumber in x. So the smaller the',...
        ' sigma value the greater the filtering. A sigma value of 1 is essentially no filtering.',...
        'The filtered Gross and filtered Detail and summed to produce the output (Result).'];
    uicontrol(gcf,'style','pushbutton','string','Info','tag','info','units','normalized',...
        'userdata',msg,'position',[xnow,ynow+1.5*ht,wid,.8*ht],'callback','seisplotsvd_foot(''info'');',...
        'backgroundcolor','y','fontsize',12);


    set(hclipsvd,'userdata',{clips,am,sigma,amax,amin});
    set(hclipsvd,'string',clipstr,'value',iclip);
    
        
    imagesc(x,y,gross,clim);
    grid
    
    title({'SVD Gross',['singcut= ' int2str(singcut)]});
    
    set(gca,'yticklabel','','xticklabel','');

    axes(hax3);
    imagesc(x,y,detail,clim);
%     brighten(.5);
    grid
    
    title({'SVD Detail',['singcut= ' int2str(singcut)]});
    
    set(gca,'yticklabel','','xticklabel','');

    
    %now do a wavenumber filter with default parameters
    sigmaxG=str2double(sigmaxvals{isigG});
    sigmayG=str2double(sigmayvals{isigG});
    xx=dx*(1:length(x));
    yy=dy*(1:length(y))';
    [grossa,aG]=wavenumberfilt(gross,xx,yy,sigmaxG,sigmayG);
    sigmaxD=str2double(sigmaxvals{isigD});
    sigmayD=str2double(sigmayvals{isigD});
    [detaila,aD]=wavenumberfilt(detail,xx,yy,sigmaxD,sigmayD);
    
    slicea=grossa+detaila;
    sliced=slice-slicea;
    ynot=.1;
    
    hax4=subplot('position',[xnot,ynot,xwid,yht]);
    imagesc(x,y,slicea,clim);
    grid
    title('Result: after footprint suppression')   
    xlabel('xline number');
    ylabel('inline number');
    
    %toggle buttons
    ynow2=ynot+yht+.6*ht;
    xnow=xnot+xwid-.5*wid;
    hbutgrp2=uibuttongroup(gcf,'units','normalized','position',[xnow,ynow2,2*wid,1.2*ht],'tag','IOD','userdata',{slice slicea sliced dname});
    uicontrol(hbutgrp2,'style','radio','string','Input & Result','units','normalized','position',[0,.6667,1,.34],...
        'enable','on','tag','in&out','backgroundcolor','w','callback','seisplotsvd_foot(''IOD'');');
    uicontrol(hbutgrp2,'style','radio','string','Input & Difference','units','normalized','position',[0,.3333,1,.34],...
        'enable','on','tag','in&diff','backgroundcolor','w','callback','seisplotsvd_foot(''IOD'');');
    uicontrol(hbutgrp2,'style','radio','string','Difference & Result','units','normalized','position',[0,0,1,.34],...
        'enable','on','tag','out&diff','backgroundcolor','w','callback','seisplotsvd_foot(''IOD'');');  
    
    %Spectra toggle buttons
    ynow2=ynot+yht+.6*ht;
    xnow=xnot+2*xwid+xsep-.5*wid;
    hbutgrp3=uibuttongroup(gcf,'units','normalized','position',[xnow,ynow2,1.5*wid,ht],'tag','specdat',...
        'userdata',{gross grossa detail detaila [aG aD] singcut [sigmaxG sigmayG sigmaxD sigmayD] x y 1},...
        'backgroundcolor','w');
    %the last entry in userdata is a 0-1 flag. 1 means data are shown 0 means spectra
    uicontrol(hbutgrp3,'style','radio','string','Show data','units','normalized','position',[0,.5,1,.5],...
        'enable','on','tag','data','backgroundcolor','w','callback','seisplotsvd_foot(''specdat'');');
    uicontrol(hbutgrp3,'style','radio','string','Show spectra','units','normalized','position',[0,0,1,.5],...
        'enable','on','tag','spectra','backgroundcolor','w','callback','seisplotsvd_foot(''specdat'');');
   
    
    hax5=subplot('position',[xnot+xwid+xsep ynot xwid yht]);
    
    imagesc(x,y,grossa,clim);
    title({['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]});
    grid
    set(gca,'yticklabel','');
    xlabel('xline number');
    
    hax6=subplot('position',[xnot+2*xwid+1.5*xsep ynot xwid yht]);
    
    imagesc(x,y,detaila,clim);
    title({['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]});
    grid
    set(gca,'yticklabel','');
    xlabel('xline number');
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
    set(hax1,'tag','slice');
    set(hax2,'tag','gross');
    set(hax3,'tag','detail');
    set(hax4,'tag','result');
    set(hax5,'tag','grossa');
    set(hax6,'tag','detaila');
    
    set(hsigmay,'userdata',{hax1 hax2 hax3 hax4 hax5 hax6});
    
    if(nargout>0)
        datar=cell(1,4);
        datar{1}=hax1;
        datar{2}=hax2;
        datar{3}=hax3;
        datar{4}=hax4;
        datar{5}=hax5;
        datar{6}=hax6;
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
    DRAGLINE_CALLBACK='seisplotsvd_foot(''singcutchange'');';
    DRAGLINE_MOTIONCALLBACK='';
    DRAGLINE_PAIRED=[];
    dragline('click');
elseif(strcmp(action,'clipxt'))
    %determine if data or spectra
    hspecdat=findobj(gcf,'tag','specdat');
    udatspecdat=get(hspecdat,'userdata');
    displayflag=udatspecdat{10};%1 for data 0 for spectra
    hclip=findobj(gcf,'tag','clipxt');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    if(displayflag==1)
        am=udat{2};
        amax=udat{4};
        sigma=udat{3};
    else
       am=udat{6};
       sigma=udat{7};
       amax=udat{8};
       amin=udat{9};
    end
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    if(iclip==1)
        if(displayflag==1)
            clim=[-amax amax];
        else
            clim=[amin amax];
        end
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set([haxes{1} haxes{4}],'clim',clim);
elseif(strcmp(action,'clipsvd'))
    %determine if data or spectra
    hspecdat=findobj(gcf,'tag','specdat');
    udatspecdat=get(hspecdat,'userdata');
    displayflag=udatspecdat{10};%1 for data 0 for spectra
    hclip=findobj(gcf,'tag','clipsvd');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    if(displayflag==1)
        am=udat{2};
        amax=udat{4};
        sigma=udat{3};
    else
       am=udat{6};
       sigma=udat{7};
       amax=udat{8};
       amin=udat{9};
    end
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    if(iclip==1)
        if(displayflag==1)
            clim=[-amax amax];
        else
            clim=[amin amax];
        end
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set([haxes{2} haxes{3} haxes{5} haxes{6}],'clim',clim);
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
elseif(strcmp(action,'singcutchange'))
    %determine singcut
    hsingcutaxe=findobj(gcf,'tag','singcutaxe');
    udat=get(hsingcutaxe,'userdata');
    singvals=udat{2};
    nsing=length(singvals);
    h1=findobj(gcf,'tag','singcut');
    val=get(h1,'ydata');
    singcut=round(val(1));
    if(singcut<1); singcut=1; end
    if(singcut>nsing); singcut=nsing; end
    %decompose
    [gross,detail,singvalsg,singvalsd]=decompose(singcut);
    
    %get the axes
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    %refresh the displays
    axes(haxes{2});
    hi=findobj(haxes{2},'type','image');
    set(hi,'cdata',gross);
    title({'SVD Gross',['singcut= ' int2str(singcut)]});
    axes(haxes{3});
    hi=findobj(haxes{3},'type','image');
    set(hi,'cdata',detail);
    hhg=findobj(gcf,'tag','singvalsg');
    set(hhg,'xdata',singvalsg);
    hhg=findobj(gcf,'tag','singvalsd');
    set(hhg,'xdata',singvalsd);
    title({'SVD Detail',['singcut= ' int2str(singcut)]});
    %update userdat in hspecdat
    hspecdat=findobj(gcf,'tag','specdat');
    udat=get(hspecdat,'userdata');
    udat{1}=gross;
    udat{3}=detail;
    udat{6}=singcut;
    set(hspecdat,'userdata',udat);
    %call sigmachange to force recompute
    seisplotsvd_foot('sigmachange');

elseif(strcmp(action,'sigmachange'))
    %determine the sigma values
    hsigmax=findobj(gcf,'tag','sigmaxG');
    vals=get(hsigmax,'string');
    ival=get(hsigmax,'value');
    sigmaxG=str2double(vals{ival});
    hsigmay=findobj(gcf,'tag','sigmayG');
    vals=get(hsigmay,'string');
    ival=get(hsigmay,'value');
    sigmayG=str2double(vals{ival});
    hsigmaxD=findobj(gcf,'tag','sigmaxD');
    vals=get(hsigmaxD,'string');
    ival=get(hsigmaxD,'value');
    sigmaxD=str2double(vals{ival});
    hsigmayD=findobj(gcf,'tag','sigmayD');
    vals=get(hsigmayD,'string');
    ival=get(hsigmayD,'value');
    sigmayD=str2double(vals{ival});
    %determine the physical spacings
    udat=get(hsigmax,'userdata');
    dx=udat{1};
    dy=udat{2};
    %get the axes handles
%    haxes=get(hsigmay,'userdata');
%     hi=findobj(haxes{2},'type','image');
%     gross=hi.CData;
%     hi=findobj(haxes{3},'type','image');
%     detail=hi.CData;
    %get gross and detail from hspecdat
    hspecdat=findobj(gcf,'tag','specdat');
    udatspecdat=get(hspecdat,'userdata');
    gross=udatspecdat{1};
    detail=udatspecdat{3};
    singcut=udatspecdat{6};
    %apply wavenumber filters
    xx=dx*(1:size(gross,2));
    yy=dy*(1:size(gross,1))';
    [grossa,aG]=wavenumberfilt(gross,xx,yy,sigmaxG,sigmayG);
    [detaila,aD]=wavenumberfilt(detail,xx,yy,sigmaxD,sigmayD);
    %get IOD userdata and update
    hiod=findobj(gcf,'tag','IOD');
    udat=get(hiod,'userdata');
    slice=udat{1};
    slicea=grossa+detaila;
    sliced=slice-slicea;
    udat{2}=slicea;
    udat{3}=sliced;
    set(hiod,'userdata',udat);
    %determine what we are showing in upper and lower axes
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    hiod=findobj(gcf,'tag','IOD');
    hIR=findobj(hiod,'tag','in&out');
    iIR=get(hIR,'value');
    hID=findobj(hiod,'tag','in&diff');
    iID=get(hID,'value');
    hRD=findobj(hiod,'tag','out&diff');
    iRD=get(hRD,'value');
    %udat=get(hiod,'userdata');
    hi1=findobj(haxes{1},'type','image');
    hi2=findobj(haxes{2},'type','image');
    hi3=findobj(haxes{3},'type','image');
    hi4=findobj(haxes{4},'type','image');
    hi5=findobj(haxes{5},'type','image');
    hi6=findobj(haxes{6},'type','image');
    if(iIR==1)
        hi1.CData=slice;
        hi2.CData=gross;
        hi3.CData=detail;
        hi4.CData=slicea;
        hi5.CData=grossa;
        hi6.CData=detaila;
        haxes{1}.Title.String=udat{4};
        haxes{2}.Title.String={'SVD Gross',['singcut= ' int2str(singcut)]};
        haxes{3}.Title.String={'SVD Detail',['singcut= ' int2str(singcut)]};
        haxes{4}.Title.String='Result: after footprint suppression';
        haxes{5}.Title.String={['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};    
    elseif(iID==1)
        hi1.CData=slice;
        hi2.CData=gross;
        hi3.CData=detail;
        hi4.CData=sliced;
        hi5.CData=gross-grossa;
        hi6.CData=detail-detaila;
        haxes{1}.Title.String=udat{4};
        haxes{2}.Title.String={'SVD Gross',['singcut= ' int2str(singcut)]};
        haxes{3}.Title.String={'SVD Detail',['singcut= ' int2str(singcut)]};
        haxes{4}.Title.String='Difference: Input - Result';
        haxes{5}.Title.String={'Difference: Gross - Filtered Gross',['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={'Difference: Detail - Filtered Detail',['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};
    elseif(iRD==1)
        hi1.CData=sliced;
        hi2.CData=gross-grossa;
        hi3.CData=detail-detaila;
        hi4.CData=slicea;
        hi5.CData=grossa;
        hi6.CData=detaila;
        haxes{1}.Title.String='Difference: Input - Result';
        haxes{2}.Title.String={'Difference: Gross - Filtered Gross',['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{3}.Title.String={'Difference: Detail - Filtered Detail',['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};
        haxes{4}.Title.String='Result: after footprint suppression';
        haxes{5}.Title.String={['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};    
    end
    %update hspecdat userdata
    udatspecdat{2}=grossa;
    udatspecdat{4}=detaila;
    udatspecdat{5}=[aG aD];
    udatspecdat{7}=[sigmaxG sigmayG sigmaxD sigmayD];
    udatspecdat{10}=1;
    set(hspecdat,'userdata',udatspecdat);
    
    seisplotsvd_foot('specdat');
    
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
elseif(strcmp(action,'eqzoom'))
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    xl=get(gca,'xlim');
    yl=get(gca,'ylim');
    set([haxes{1} haxes{2} haxes{3} haxes{4} haxes{5} haxes{6}],'xlim',xl,'ylim',yl)
elseif(strcmp(action,'unzoom'))
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    hi=findobj(haxes{1},'type','image');
    x=hi.XData;
    y=hi.YData;
    set([haxes{1} haxes{2} haxes{3} haxes{4} haxes{5} haxes{6}],'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
elseif(strcmp(action,'IOD'))
    %get axes
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    %determine viewing flag
    hiod=findobj(gcf,'tag','IOD');
    hIR=findobj(hiod,'tag','in&out');
    iIR=get(hIR,'value');
    hID=findobj(hiod,'tag','in&diff');
    iID=get(hID,'value');
    hRD=findobj(hiod,'tag','out&diff');
    iRD=get(hRD,'value');
    %unpack specdat userdata
    hspecdat=findobj(gcf,'tag','specdat');
    udatspecdat=get(hspecdat,'userdata');
    gross=udatspecdat{1};
    grossa=udatspecdat{2};
    detail=udatspecdat{3};
    detaila=udatspecdat{4};
    aG=udatspecdat{5}(1);
    aD=udatspecdat{5}(2);
    singcut=udatspecdat{6};
    sigmaxG=udatspecdat{7}(1);
    sigmayG=udatspecdat{7}(2);
    sigmaxD=udatspecdat{7}(3);
    sigmayD=udatspecdat{7}(4);
    udatspecdat{10}=1;
    set(hspecdat,'userdata',udatspecdat);
    %unpack IOD userdata
    udat=get(hiod,'userdata');
    slice=udat{1};
    slicea=udat{2};
    sliced=udat{3};
    %get the image handles
    hi1=findobj(haxes{1},'type','image');
    hi2=findobj(haxes{2},'type','image');
    hi3=findobj(haxes{3},'type','image');
    hi4=findobj(haxes{4},'type','image');
    hi5=findobj(haxes{5},'type','image');
    hi6=findobj(haxes{6},'type','image');
    %update the views
    if(iIR==1)
        hi1.CData=slice;
        hi2.CData=gross;
        hi3.CData=detail;
        hi4.CData=slicea;
        hi5.CData=grossa;
        hi6.CData=detaila;
        haxes{1}.Title.String=udat{4};
        haxes{2}.Title.String={'SVD Gross',['singcut= ' int2str(singcut)]};
        haxes{3}.Title.String={'SVD Detail',['singcut= ' int2str(singcut)]};
        haxes{4}.Title.String='Result: after footprint suppression';
        haxes{5}.Title.String={['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};    
    elseif(iID==1)
        hi1.CData=slice;
        hi2.CData=gross;
        hi3.CData=detail;
        hi4.CData=sliced;
        hi5.CData=gross-grossa;
        hi6.CData=detail-detaila;
        haxes{1}.Title.String=udat{4};
        haxes{2}.Title.String={'SVD Gross',['singcut= ' int2str(singcut)]};
        haxes{3}.Title.String={'SVD Detail',['singcut= ' int2str(singcut)]};
        haxes{4}.Title.String='Difference: Input - Result';
        haxes{5}.Title.String={'Difference: Gross - Filtered Gross',['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={'Difference: Detail - Filtered Detail',['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};
    elseif(iRD==1)
        hi1.CData=sliced;
        hi2.CData=gross-grossa;
        hi3.CData=detail-detaila;
        hi4.CData=slicea;
        hi5.CData=grossa;
        hi6.CData=detaila;
        haxes{1}.Title.String='Difference: Input-Result';
        haxes{2}.Title.String={'Difference: Gross - Filtered Gross',['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{3}.Title.String={'Difference: Detail - Filtered Detail',['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};
        haxes{4}.Title.String='Result: after footprint suppression';
        haxes{5}.Title.String={['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};    
    end
    seisplotsvd_foot('specdat');
elseif(strcmp(action,'specdat'))
    %get specdat userdata
    hspecdat=findobj(gcf,'tag','specdat');
    udatspecdat=get(hspecdat,'userdata');
    displayflag=udatspecdat{10};%1 for data, 0 for spectra, this is what is being displayed
    gross=udatspecdat{1};
    grossa=udatspecdat{2};
    grossd=gross-grossa;
    detail=udatspecdat{3};
    detaila=udatspecdat{4};
    detaild=detail-detaila;
    aG=udatspecdat{5}(1);
    aD=udatspecdat{5}(2);
    singcut=udatspecdat{6};
    sigmaxG=udatspecdat{7}(1);
    sigmayG=udatspecdat{7}(2);
    sigmaxD=udatspecdat{7}(3);
    sigmayD=udatspecdat{7}(4);
    x=udatspecdat{8};
    y=udatspecdat{9};
    %get IOD userdata
    hiod=findobj(gcf,'tag','IOD');
    udat=get(hiod,'userdata');
    slice=udat{1};
    slicea=udat{2};
    sliced=udat{3};
    %determine what is requested
    hdata=findobj(hspecdat,'tag','data');
    idisp=get(hdata,'value');%1 for data 0 for spectra
    if(displayflag==idisp)
        return;
    end
    %get the axes
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    if(displayflag==1 && idisp==0)
        pctfk=10;
        %transform the slices and change coordinates
        [tmp,ky,kx]=fktran(slice,y,x,nan,nan,pctfk);
        slice=abs(tmp);
        hclip=findobj(gcf,'tag','clipxt');
        udat=get(hclip,'userdata');
        clips=udat{1};
        iclip=near(clips,3);
        clip=clips(iclip);
        set(hclip,'value',iclip);
        hclipsvd=findobj(gcf,'tag','clipsvd');
        set(hclipsvd,'value',iclip);
        if(length(udat)==5)
            %this is the first time the spectra button has been pressed
            udat{9}=min(slice(:));
            udat{8}=max(slice(:));
            udat{7}=std(slice(:));
            udat{6}=mean(slice(:));
            set(hclip,'userdata',udat);
            set(hclipsvd,'userdata',udat);
            %save the xdir and ydir settings for the spatial axes
            xdir=get(haxes{1},'xdir');
            ydir=get(haxes{1},'ydir');
            set(hdata,'userdata',{xdir ydir});
        end
        am=udat{6};
        sigma=udat{7};
        tmp=fktran(slicea,y,x,nan,nan,pctfk);
        slicea=abs(tmp);
        tmp=fktran(sliced,y,x,nan,nan,pctfk);
        sliced=abs(tmp);
        tmp=fktran(gross,y,x,nan,nan,pctfk);
        gross=abs(tmp);
        tmp=fktran(grossa,y,x,nan,nan,pctfk);
        grossa=abs(tmp);
        tmp=fktran(grossd,y,x,nan,nan,pctfk);
        grossd=abs(tmp);
        tmp=fktran(detail,y,x,nan,nan,pctfk);
        detail=abs(tmp);
        tmp=fktran(detaila,y,x,nan,nan,pctfk);
        detaila=abs(tmp);
        tmp=fktran(detaild,y,x,nan,nan,pctfk);
        detaild=abs(tmp);
        x=kx;
        y=ky;
        clim=[am-clip*sigma am+clip*sigma];
        xname='xline wavenumber';
        yname='yline wavenumber';
        %update display flag
        udatspecdat{10}=0;
        set(hspecdat,'userdata',udatspecdat);
        xdir='normal';
        ydir='reverse';
    else
        hclip=findobj(gcf,'tag','clipxt');
        udat=get(hclip,'userdata');
        clips=udat{1};
        iclip=near(clips,3);
        clip=clips(iclip);
        am=udat{2};
        sigma=udat{3};
        clim=[am-clip*sigma am+clip*sigma];
        set(hclip,'value',iclip);
        hclipsvd=findobj(gcf,'tag','clipsvd');
        set(hclipsvd,'value',iclip);
        xname='xline number';
        yname='yline number';
        %update display flag
        udatspecdat{10}=1;
        set(hspecdat,'userdata',udatspecdat);
        uudat=get(hdata,'userdata');
        xdir=uudat{1};
        ydir=uudat{2};
    end

    %get the axes settings to restore later
    fs=get(haxes{1},'fontsize');
    xgrid=get(haxes{1},'xgrid');
    ygrid=get(haxes{1},'ygrid');
    gridcolor=get(haxes{1},'gridcolor');
    gridalpha=get(haxes{1},'gridalpha');
    %determine viewing flag
    hiod=findobj(gcf,'tag','IOD');
    hIR=findobj(hiod,'tag','in&out');
    iIR=get(hIR,'value');
    hID=findobj(hiod,'tag','in&diff');
    iID=get(hID,'value');
    hRD=findobj(hiod,'tag','out&diff');
    iRD=get(hRD,'value');
    udat=get(hiod,'userdata');
    if(iIR==1)
        axes(haxes{1});
        imagesc(x,y,slice,clim);
        ylabel(yname);
        set(gca,'xticklabel','');
        axes(haxes{2})
        imagesc(x,y,gross,clim);
        set(gca,'xticklabel','','yticklabel','');
        axes(haxes{3})
        imagesc(x,y,detail,clim);
        set(gca,'xticklabel','','yticklabel','');
        axes(haxes{4})
        imagesc(x,y,slicea,clim);
        ylabel(yname);xlabel(xname);
        axes(haxes{5})
        imagesc(x,y,grossa,clim);
        xlabel(xname);
        set(gca,'yticklabel','');
        axes(haxes{6})
        imagesc(x,y,detaila,clim);
        xlabel(xname);
        set(gca,'yticklabel','');
        haxes{1}.Title.String=udat{4};
        haxes{2}.Title.String={'SVD Gross',['singcut= ' int2str(singcut)]};
        haxes{3}.Title.String={'SVD Detail',['singcut= ' int2str(singcut)]};
        haxes{4}.Title.String='Result: after footprint suppression';
        haxes{5}.Title.String={['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};    
    elseif(iID==1)
        axes(haxes{1});
        imagesc(x,y,slice,clim);
        ylabel(yname);
        set(gca,'xticklabel','');
        axes(haxes{2})
        imagesc(x,y,gross,clim);
        set(gca,'xticklabel','','yticklabel','');
        axes(haxes{3})
        imagesc(x,y,detail,clim);
        set(gca,'xticklabel','','yticklabel','');
        axes(haxes{4})
        imagesc(x,y,sliced,clim);
        ylabel(yname);xlabel(xname);
        axes(haxes{5})
        imagesc(x,y,grossd,clim);
        xlabel(xname);
        set(gca,'yticklabel','');
        axes(haxes{6})
        imagesc(x,y,detaild,clim);
        xlabel(xname);
        set(gca,'yticklabel','');
        haxes{1}.Title.String=udat{4};
        haxes{2}.Title.String={'SVD Gross',['singcut= ' int2str(singcut)]};
        haxes{3}.Title.String={'SVD Detail',['singcut= ' int2str(singcut)]};
        haxes{4}.Title.String='Difference: Input - Result';
        haxes{5}.Title.String={'Difference: Gross - Filtered Gross',['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={'Difference: Detail - Filtered Detail',['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};
    elseif(iRD==1)
        axes(haxes{1});
        imagesc(x,y,sliced,clim);
        ylabel(yname);
        set(gca,'xticklabel','');
        axes(haxes{2})
        imagesc(x,y,grossd,clim);
        set(gca,'xticklabel','','yticklabel','');
        axes(haxes{3})
        imagesc(x,y,detaild,clim);
        set(gca,'xticklabel','','yticklabel','');
        axes(haxes{4})
        imagesc(x,y,slicea,clim);
        ylabel(yname);xlabel(xname);
        axes(haxes{5})
        imagesc(x,y,grossa,clim);
        xlabel(xname);
        set(gca,'yticklabel','');
        axes(haxes{6})
        imagesc(x,y,detaila,clim);
        xlabel(xname);
        set(gca,'yticklabel','');
        haxes{1}.Title.String='Difference: Input-Result';
        haxes{2}.Title.String={'Difference: Gross - Filtered Gross',['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{3}.Title.String={'Difference: Detail - Filtered Detail',['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};
        haxes{4}.Title.String='Result: after footprint suppression';
        haxes{5}.Title.String={['Filtered Gross, scalar=' num2str(aG,3)],['sigmax=' num2str(sigmaxG) ...
        ', sigmay=' num2str(sigmayG)]};
        haxes{6}.Title.String={['Filtered Detail, scalar=' num2str(aD,3)],['sigmax=' num2str(sigmaxD) ...
        ', sigmay=' num2str(sigmayD)]};    
    end
    set([haxes{1} haxes{2} haxes{3} haxes{4} haxes{5} haxes{6}],'fontsize',fs,'xgrid',xgrid,...
        'ygrid',ygrid,'gridcolor',gridcolor,'gridalpha',gridalpha,'xdir',xdir,'ydir',ydir);
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
       transfer='seisplotsvd_foot(''gridcc2'');';
       %the questions
       q1='Physical distance between xlines';
       q2='Physical distance between inlines';
       q3='Acquisition line spacing in xline direction';
       q4='Acquisition line spacing in inline direction';
       q=char(q1,q2,q3,q4);
       askthingsinit(transfer,q);
       return;
    end
    dx=udat(1);
    dy=udat(2);
    dxaq=udat(3);
    dyaq=udat(4);
    %get the six sections
    hthisfig=gcf;
    hiod=findobj(hthisfig,'tag','IOD');
    udat=get(hiod,'userdata');
    slice=udat{1};
    slicea=udat{2};
    dname=udat{4};
    hspec=findobj(hthisfig,'tag','specdat');
    udat=get(hspec,'userdata');
    gross=udat{1};
    grossa=udat{2};
    detail=udat{3};
    detaila=udat{4};
    x=udat{8};
    y=udat{9};
    %get the display mode
    hIR=findobj(hiod,'tag','in&out');
    iIR=get(hIR,'value');
    hID=findobj(hiod,'tag','in&diff');
    iID=get(hID,'value');
    hRD=findobj(hiod,'tag','out&diff');
    iRD=get(hRD,'value');
    hsigmay=findobj(gcf,'tag','sigmayG');
    haxes=get(hsigmay,'userdata');
    if(iIR==1)
        [ccx1,xlags1,ccy1,ylags1,xaqgrid,yaqgrid]=ccfoot(slice,x*dx,y*dy,dxaq,dyaq);
        [ccx2,xlags2,ccy2,ylags2]=ccfoot(gross,x*dx,y*dy,dxaq,dyaq);
        [ccx3,xlags3,ccy3,ylags3]=ccfoot(detail,x*dx,y*dy,dxaq,dyaq);
        [ccx4,xlags4,ccy4,ylags4]=ccfoot(slicea,x*dx,y*dy,dxaq,dyaq);
        [ccx5,xlags5,ccy5,ylags5]=ccfoot(grossa,x*dx,y*dy,dxaq,dyaq);
        [ccx6,xlags6,ccy6,ylags6]=ccfoot(detaila,x*dx,y*dy,dxaq,dyaq);
    end
    if(iID==1)
        [ccx1,xlags1,ccy1,ylags1,xaqgrid,yaqgrid]=ccfoot(slice,x*dx,y*dy,dxaq,dyaq);
        [ccx2,xlags2,ccy2,ylags2]=ccfoot(gross,x*dx,y*dy,dxaq,dyaq);
        [ccx3,xlags3,ccy3,ylags3]=ccfoot(detail,x*dx,y*dy,dxaq,dyaq);
        [ccx4,xlags4,ccy4,ylags4]=ccfoot(slice-slicea,x*dx,y*dy,dxaq,dyaq);
        [ccx5,xlags5,ccy5,ylags5]=ccfoot(gross-grossa,x*dx,y*dy,dxaq,dyaq);
        [ccx6,xlags6,ccy6,ylags6]=ccfoot(detail-detaila,x*dx,y*dy,dxaq,dyaq);
    end
    if(iRD==1)
        [ccx1,xlags1,ccy1,ylags1,xaqgrid,yaqgrid]=ccfoot(slice-slicea,x*dx,y*dy,dxaq,dyaq);
        [ccx2,xlags2,ccy2,ylags2]=ccfoot(gross-grossa,x*dx,y*dy,dxaq,dyaq);
        [ccx3,xlags3,ccy3,ylags3]=ccfoot(detail-detaila,x*dx,y*dy,dxaq,dyaq);
        [ccx4,xlags4,ccy4,ylags4]=ccfoot(slicea,x*dx,y*dy,dxaq,dyaq);
        [ccx5,xlags5,ccy5,ylags5]=ccfoot(grossa,x*dx,y*dy,dxaq,dyaq);
        [ccx6,xlags6,ccy6,ylags6]=ccfoot(detaila,x*dx,y*dy,dxaq,dyaq);
    end
    dname1=haxes{1}.Title.String;
    dname2=haxes{2}.Title.String;
    dname3=haxes{3}.Title.String;
    dname4=haxes{4}.Title.String;
    dname5=haxes{5}.Title.String;
    dname6=haxes{6}.Title.String;
    %new figure and plot
    yup=.05;
    ydown=.05;
    ydiv=.05;
    hccfig=figure;
    hax1=subplot(2,3,1);
    plot(xlags1,ccx1,ylags1,ccy1);
    y0=min([ccx1(:);ccy1(:)]);
    y0=max([0, floor(y0/ydiv)*ydiv-ydown]);
    y1=max([ccx1(:);ccy1(:)]);
    y1=min([1, ceil(y1/ydiv)*ydiv+yup]);
    xlabel('lag');
    ylabel('correlation')
    title(dname1);
    legend('xline correlation','inline correlation','location','south');
    xl=get(gca,'xlim');
    delx=xl(2)/2;
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    
    hax2=subplot(2,3,2);
    plot(xlags2,ccx2,ylags2,ccy2);
    y0a=min([ccx2(:);ccy2(:)]);
    y0a=max([0, floor(y0a/ydiv)*ydiv-ydown]);
    y1a=max([ccx2(:);ccy2(:)]);
    y1a=min([1, ceil(y1a/ydiv)*ydiv+yup]);
    if(y0a<y0);y0=y0a;end
    if(y1a>y1);y1=y1a;end
    xlabel('lag');
    ylabel('correlation')
    title(dname2);
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    
    hax3=subplot(2,3,3);
    plot(xlags3,ccx3,ylags3,ccy3);
    y0a=min([ccx3(:);ccy3(:)]);
    y0a=max([0, floor(y0a/ydiv)*ydiv-ydown]);
    y1a=max([ccx3(:);ccy3(:)]);
    y1a=min([1, ceil(y1a/ydiv)*ydiv+yup]);
    if(y0a<y0);y0=y0a;end
    if(y1a>y1);y1=y1a;end
    xlabel('lag');
    ylabel('correlation')
    title(dname3);
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    
    hax4=subplot(2,3,4);
    plot(xlags4,ccx4,ylags4,ccy4);
    y0a=min([ccx4(:);ccy4(:)]);
    y0a=max([0, floor(y0a/ydiv)*ydiv-ydown]);
    y1a=max([ccx4(:);ccy4(:)]);
    y1a=min([1, ceil(y1a/ydiv)*ydiv+yup]);
    if(y0a<y0);y0=y0a;end
    if(y1a>y1);y1=y1a;end
    xlabel('lag');
    ylabel('correlation')
    title(dname4);
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    
    hax5=subplot(2,3,5);
    plot(xlags5,ccx5,ylags5,ccy5);
    y0a=min([ccx5(:);ccy5(:)]);
    y0a=max([0, floor(y0a/ydiv)*ydiv-ydown]);
    y1a=max([ccx5(:);ccy5(:)]);
    y1a=min([1, ceil(y1a/ydiv)*ydiv+yup]);
    if(y0a<y0);y0=y0a;end
    if(y1a>y1);y1=y1a;end
    xlabel('lag');
    ylabel('correlation')
    title(dname5);
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    
    hax6=subplot(2,3,6);
    plot(xlags6,ccx6,ylags6,ccy6);
    y0a=min([ccx6(:);ccy6(:)]);
    y0a=max([0, floor(y0a/ydiv)*ydiv-ydown]);
    y1a=max([ccx6(:);ccy6(:)]);
    y1a=min([1, ceil(y1a/ydiv)*ydiv+yup]);
    if(y0a<y0);y0=y0a;end
    if(y1a>y1);y1=y1a;end
    xlabel('lag');
    ylabel('correlation')
    title(dname6);
    set(gca,'xtick',xl(1):delx:xl(2));
    grid
    
    set([hax1 hax2 hax3 hax4 hax5 hax6],'ylim',[y0 y1]);
    
    boldlines;
    bigfont(gcf,1.5,1);
    pos=get(hthisfig,'position');
    wid=pos(3)*1;
    ht=pos(4)*.75;
    if(iscell(dname))
        dname= dname{1};
    end
    dname=strrep(dname,'Input: ','');
    set(hccfig,'position',[pos(1)+.5*(pos(3)-wid),pos(2)+.5*(pos(4)-ht),wid,ht],'name',['Grid correlations for ' dname])
    titlefontsize(.8,1,hccfig)
    %showgridsbutton
    uicontrol(hccfig,'style','pushbutton','string','Show grids','units','normalized',...
        'position',[.91,.85,.08,.05],'callback','seisplotsvd_foot(''showgrids'');',...
        'userdata',{xaqgrid yaqgrid x y dname},'tag','showgrids');
    
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
end
end

function [gross,detail,singvalsg,singvalsd]=decompose(singcut)
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

function [grossa,a]=wavenumberfilt(gross,xx,yy,sigmax,sigmay)

%grossa=wavenumber_gaussmask(gross,xx,yy,sigmax,sigmay);
grossa=wavenumber_gaussmask2(gross,sigmax,sigmay);

[~,a]=lsqsubtract(gross(:),grossa(:));
grossa=a*grossa;

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