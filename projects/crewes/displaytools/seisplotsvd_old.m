function datar=seisplotsvd_old_old(seis,t,x,dname)
% seisplotsvd_old: plots a seismic gather and its SVD filtered cousin side-by-side
%
% datar=seisplotsvd_old(seis,t,x,dname)
%
% A new figure is created and divided into two same-sized axes (side-by-side). The seismic
% gather (matrix) is plotted as an image in the left-hand-side and its SVD
% filtered cousin is plotted as an image in the right-hand-side. Controls
% are provided to adjust the singular value threshold
% and to brighten or darken the image plots. The data should be regularly sampled in both t and
% x.
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
%           data{2} ... handle of the svd axes
%           data{3} ... singular values
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
    brighten(.5);
    grid
    dx=x(2)-x(1);
    dt=t(2)-t(1);
    title({dname ,['x-t space dx=' num2str(dx) ', dt=' num2str(dt)]})
    
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
    %make a clip control

    xnow=xnot+xwid;
    wid=.055;ht=.05;sep=.005;
    ynow=ynot+yht-ht;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clipxt','units','normalized',...
        'position',[xnow,ynow,1.2*wid,ht],'callback','seisplotsvd_old(''clipxt'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brightenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_old(''brightenxt'');',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darkenxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotsvd_old(''brightenxt'');',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightnessxt','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    
    set(hax1,'tag','seis');
    
    hax2=subplot('position',[xnot+xwid+xsep ynot xwid yht]);
    [U,S,V]=svd(seis);
    singvals=diag(S);
    thresh=20;%threshold as a percentage of the maximum value
    %make a clip control
    xnow=xnot+2*xwid+xsep;
    ht=.05;
    ynow=ynot+yht-ht;
    %wid=.045;sep=.005;
    hclipsvd=uicontrol(gcf,'style','popupmenu','string','xxx','tag','clipsvd','units','normalized',...
        'position',[xnow,ynow,1.2*wid,ht],'callback','seisplotsvd_old(''clipsvd'');','value',1,...
        'userdata',[],'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped');
    ht=.025;
    ynow=ynow-ht-sep;
    sep=.1*wid;
    uicontrol(gcf,'style','text','string','Threshold:','tag','svdlabel','units','normalized',...
        'position',[xnow,ynow,.6*wid,ht],'tooltipstring',...
        'Expressed as a percent of the maximum singular value');
    uicontrol(gcf,'style','edit','string',num2str(thresh),'tag','thresh','units','normalized',...
        'callback','seisplotsvd_old(''threshold'')','userdata',{U,singvals,V,thresh,dname},...
        'position',[xnow+.6*wid+sep,ynow,.5*wid,ht],'tooltipstring',...
        'Enter a value between 100 and 0');

    [seissvd,numused]=reconstruct;

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seissvd);
    set(hclipsvd,'userdata',{clips,am,sigma,amax,amin,hax2});
    set(hclipsvd,'string',clipstr,'value',iclip);
    
    %clim=[amin am+clip*sigma];
    if(iclip==1)
        clim=[-amax amax];
    else
        clim=[am-clip*sigma am+clip*sigma];
    end
        
    imagesc(x,t,seissvd,clim);colormap(seisclrs);
    brighten(.5);
    grid
    
    title({[dname 'SVD reconstruction'],['Using ' int2str(numused) ' of the ' int2str(length(singvals)) ' singular values']});
    
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
    
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    
    set(hax2,'tag','svd');
    
    set(gcf,'name',dname);
    if(nargout>0)
        datar=cell(1,4);
        datar{1}=hax1;
        datar{2}=hax2;
        datar{3}=singvals;
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
elseif(strcmp(action,'clipsvd'))
    hclip=findobj(gcf,'tag','clipsvd');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
    %amin=udat{5};
    sigma=udat{3};
    hax=udat{6};
    if(iclip==1)
        %clim=[amin amax];
        clim=[-amax amax];
    else
        clip=clips(iclip-1);
        clim=[am-clip*sigma,am+clip*sigma];
        %clim=[amin am+clip*sigma];
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
elseif(strcmp(action,'threshold'))
    hthresh=findobj(gcf,'tag','thresh');
    str=get(hthresh,'string');
    val=str2double(str);
    if(isnan(val) || val<0 || val>100)
       msgbox('Invalid threshold value. Must be a number between 0 and 100'); 
    end
    thresh=val;
    udat=get(hthresh,'userdata');
    udat{4}=thresh;
    set(hthresh,'userdata',udat);
    dname=udat{5};
    singvals=udat{2};
    [seissvd,numused]=reconstruct(thresh);
    haxesvd=findobj(gcf,'tag','svd');
    hi=findobj(haxesvd,'type','image');
    set(hi,'cdata',seissvd);
    title({[dname 'SVD reconstruction'],['Using ' num2str(numused) ' of the ' int2str(length(singvals)) ' singular values']});
end
end

function [seissvd,numused]=reconstruct(thresh)
hthresh=findobj(gcf,'tag','thresh');
udat=get(hthresh,'userdata');
U=udat{1};
singvals=udat{2};
V=udat{3};
if(nargin<1)
    thresh=udat{4}/100;
else
    thresh=thresh/100;
end
ind=find(singvals>thresh*singvals(1));
newsingvals=zeros(size(singvals));
newsingvals(ind)=singvals(ind);
numused=length(ind);
m=size(U,1);
n=size(V,1);
if(m>n)
    S=[diag(newsingvals);zeros(m-n,n)];
elseif(n>m)
    S=[diag(newsingvals) zeros(m,n-m)];
else
    S=diag(newsingvals);
end
seissvd=U*S*V';
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