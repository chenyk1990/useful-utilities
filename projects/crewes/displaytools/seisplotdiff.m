function datar=seisplotdiff(seis1,seis2,t,x,dname1,dname2,xname,yname)
% seisplotdiff: Plots two input seismic gathers and their difference side-by-side in true scale
%
% datar=seisplotdiff(seis1,seis2,t,x,dname1,dname2)
%
% A new figure is created and divided into three same-sized axes (side-by-side). The two input
% gathers are displays as well as their difference in true relative amplitude. A graphical control
% allows the selection of a simple difference or a least-squares subtraction computed by scaling
% seis2 to minimize the difference. 
%
% seis1 ... input seismic matrix #1
% seis2 ... input seismic matrix #2
% t ... time coordinate vector for seis (y or row coordinate). Only used
%       for plotting
% *********** default 1:nrows ***************
% x ... space coordinate vector for seis (x or column coordinate)
% *********** default 1:ncols ***************
% dname1 ... text string giving a name for first dataset
% ************ default dname ='dataset #1' ************
% dname2 ... text string giving a name for second dataset
% ************ default dname ='dataset #2' ************
%
% datar ... Return data which is a length 3 cell array containing
%           data{1} ... handle of the seis1 axes
%           data{2} ... handle of the seis2 axes
%           data{3} ... handle of the difference axes
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


if(~ischar(seis1))
    action='init';
else
    action=seis1;
end

datar=[];%initialize return data to null

if(strcmp(action,'init'))
    
    [nrows,ncols]=size(seis1);
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
    
    if(nargin<5)
        dname1{2}='Amplitudes unchanged';
        dname1{1}='Dataset #1';
    else
        if(iscell(dname1))
            dname1{1}=['#1 ' dname1{1}];
            dname1{2}='Amplitudes unchanged';
        else
            tmp=dname1;
            dname1=cell(1,2);
            dname1{1}=['#1 ' tmp];
            dname1{2}='Amplitudes unchanged';
        end
    end
    if(nargin<6)
        dname2{2}='Amplitudes unchanged';
        dname2{1}='Dataset #2';
    else
        if(iscell(dname2))
            dname2{1}=['#2 ' dname2{1}];
            dname2{2}='Amplitudes unchanged ';
        else
            tmp=dname2;
            dname2=cell(1,2);
            dname2{1}=['#2 ' tmp];
            dname2{2}='Amplitudes unchanged';
        end
    end
    maxmeters=7000;
    if(nargin<7)
        if(max(x)<maxmeters)
            xname='distance (m)';
        else
            xname='distance (ft)';
        end
    end
    if(nargin<8)
        if(max(t)<10)
            yname='time (s)';
        elseif(max(t)<maxmeters)
            yname='distance (m)';
        else
            yname='distance (ft)';
        end
    end
    
    xwid=.25;
    yht=.75;
    xsep=.05;
    xnot=.05;
    ynot=.1;
    factor=1.2;

    figure
    hax1=subplot('position',[xnot+.02 ynot factor*xwid yht]);

    [clips,clipstr,clip,iclip,sigma,am,amax,amin]=getclips(seis1);
    if(iclip==1)
        clim=[-amax amax];
    else
        clim=[am-clip*sigma am+clip*sigma];
    end
        
    imagesc(x,t,seis1,clim);colormap(seisclrs)
%     brighten(.5);
    grid
    xlabel(xname);
    ylabel(yname);
    
    %process dname if it is too long
    toolong=50;
    if(iscell(dname1))
        if(length(dname1{1})>toolong)
            str1=dname1{1};
            str2=dname1{2};
            ind=isspace(str1);
            ind2=find(ind>0);%points to word breaks
            ind3=find(ind2<toolong);
            if(~isempty(ind3))
               str1a=str1(1:ind2(ind3(end)));
               str2a=[str1(ind2(ind3(end))+1:end) ' ' str2];
               dname1{1}=str1a;
               dname1{2}=str2a;
            end
        end
    end
            
    title(dname1,'interpreter','none')
    
    colorbar
    
    

    
    
    set(hax1,'tag','seis1');
    
    hax2=subplot('position',[xnot+factor*xwid+xsep ynot xwid yht]);
    
    imagesc(x,t,seis2,clim);colormap(seisclrs)
%     brighten(.5);
    grid
    
    %process dname if it is too long
    if(iscell(dname2))
        if(length(dname2{1})>toolong)
            str1=dname2{1};
            str2=dname2{2};
            ind=isspace(str1);
            ind2=find(ind>0);%points to word breaks
            ind3=find(ind2<toolong);
            if(~isempty(ind3))
               str1a=str1(1:ind2(ind3(end)));
               str2a=[str1(ind2(ind3(end))+1:end) ' ' str2];
               dname2{1}=str1a;
               dname2{2}=str2a;
            end
        end
    end
            
    title(dname2,'interpreter','none')
    
    xlabel(xname);
    set(hax2,'tag','seis2','yticklabel','');
    
    hax3=subplot('position',[xnot+(1+factor)*xwid+1.5*xsep ynot xwid yht]);
    seisd=seis1-seis2;
    [~,a]=lsqsubtract(seis1(:),seis2(:));
   
    imagesc(x,t,seisd,clim);colormap(seisclrs)
%     brighten(.5);
    grid
    P1=sum(seis1(:).^2);
    PD=sum(seisd(:).^2);
    pctd=round(1000*PD/P1)/10;  
    title( {'Ordinary difference',['Difference power is ' num2str(pctd) '% of dataset #1']},'interpreter','none')
    
    xlabel(xname)
    set(hax3,'tag','seisd','yticklabel','');

    %make a clip control
    xnow=xnot+(factor+2)*xwid+1.5*xsep;
    ht=.05;
    ynow=ynot+yht-ht;
    wid=.04;ht=.05;sep=.005;
    uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdiff(''clip'');','value',iclip,...
        'userdata',{clips,am,sigma,amax,amin,hax1,hax2,hax3,seis1,seis2,a},'tooltipstring',...
        'clip level is the number of standard deviations from the mean at which amplitudes are clipped')
    
    ht=.5*ht;
    ynow=ynow-sep;
    uicontrol(gcf,'style','pushbutton','string','brighten','tag','brighten','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdiff(''brighten'');',...
        'tooltipstring','push once or multiple times to brighten the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','pushbutton','string','darken','tag','darken','units','normalized',...
        'position',[xnow,ynow,wid,ht],'callback','seisplotdiff(''brighten'');',...
        'tooltipstring','push once or multiple times to darken the images');
    ynow=ynow-ht-sep;
    uicontrol(gcf,'style','text','string','lvl 0','tag','brightness','units','normalized',...
        'position',[xnow,ynow,wid,ht],...
        'tooltipstring','image brightness (both images)','userdata',0);
    ynow2=ynow-ht-4*sep;
    hbutgrp=uibuttongroup(gcf,'units','normalized','position',[xnow,ynow2,1.5*wid,1.5*ht],'tag','option');
    uicontrol(hbutgrp,'style','radio','string','Ordinary diff','units','normalized','position',[0,.5,1,.5],...
        'enable','on','tag','ord','backgroundcolor','w','callback',@subtract);
    uicontrol(hbutgrp,'style','radio','string','Least sqs diff','units','normalized','position',[0,0,1,.5],...
        'enable','on','tag','lsq','backgroundcolor','w','callback',@subtract);
    
    %make an info button
%     msg=['The seismic matrix on the left is shown separated into its Gross structure and its Detail. ',...
%         'This separation is done using SVD (singular-value decomposition) and is controlled by the parameter singcut. ',...
%         'There are typically hundreds to thousands of singular values in an image where, sorted from largest to smallest, ',...
%         'the first few control the Gross structure and the remainder control the detail. These are plotted in the ',...
%         'axis on the right. You can adjust the cutoff singular value (singcut) by clicking and dragging the horizontal ',...
%         'red line in the SingVals axis. When you release the line a new separation will be displayed.',...
%         'The singular values defining the Gross structure are the original singular values multiplied by a Gaussian ',...
%         'centered on the largest singular value and whose standard deviation is singcut. The Detail singular values are the original singular values ',...
%         'minus those that define the Gross structure. In this way the original seismic matrix is always equal ',...
%         'to the sum of Gross and Detail.'];
%     uicontrol(gcf,'style','pushbutton','string','Info','tag','info','units','normalized',...
%         'userdata',msg,'position',[xnow,ynow+2*ht,wid,ht],'callback','seisplotdiff(''info'');',...
%         'backgroundcolor','y','fontsize',16);

    
    
    bigfig; %enlarge the figure to get more pixels
    bigfont(gcf,1.6,1); %enlarge the fonts in the figure
    boldlines(gcf,4,2); %make lines and symbols "fatter"
    whitefig;
    titlefontsize(.95)
    
    set(gcf,'name',['Difference of ' dname1{1} ' & ' dname2{1}])
    
    if(nargout>0)
        datar=cell(1,3);
        datar{1}=hax1;
        datar{2}=hax2;
        datar{3}=hax3;
    end
elseif(strcmp(action,'clip'))
    hclip=findobj(gcf,'tag','clip');
    udat=get(hclip,'userdata');
    iclip=get(hclip,'value');    
    clips=udat{1};
    am=udat{2};
    amax=udat{4};
   % amin=udat{5};
    sigma=udat{3};
    hax1=udat{6};
    hax2=udat{7};
    hax3=udat{8};
    if(iclip==1)
        clim=[-amax amax];
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set([hax1 hax2 hax3],'clim',clim);

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

% elseif(strcmp(action,'info'))
%     hinfo=findobj(gcf,'tag','info');
%     udat=get(hinfo,'userdata');
%     if(iscell(udat))
%         msg=udat{1};
%         h=udat{2};
%         if(isgraphics(h))
%             delete(h);
%         end
%     else
%         msg=udat;
%     end
%     h=msgbox(msg,'SVD Separation','help');
%     set(hinfo,'userdata',{msg,h});
end
end

function subtract(~,~)
hclip=findobj(gcf,'tag','clip');
udat=get(hclip,'userdata');
seis1=udat{9};
seis2=udat{10};
a=udat{11};
hax2=udat{7};
hi2=findobj(hax2,'type','image');
hax3=udat{8};
hi3=findobj(hax3,'type','image');
hopt=findobj(gcf,'tag','ord');
opt=get(hopt,'value');
if(opt==1)
    %ordinary subtraction
    seisd=seis1-seis2;
    set(hi2,'cdata',seis2);
    ht=get(hax2,'title');
    dname2=get(ht,'string');
    dname2{2}='Amplitudes unchanged';
    set(ht,'string',dname2);
    set(hi3,'cdata',seisd);
    ht=get(hax3,'title');
    P1=sum(seis1(:).^2);
    PD=sum(seisd(:).^2);
    pctd=round(1000*PD/P1)/10;
    set(ht,'string',{'Ordinary difference',['Difference power is ' num2str(pctd) '% of dataset #1']});
else
    %lsq sub
    seisd=seis1-a*seis2;
    set(hi2,'cdata',a*seis2);
    ht=get(hax2,'title');
    dname2=get(ht,'string');
    dname2{2}=['Amplitudes scaled by ' num2str(a,3)];
    set(ht,'string',dname2);
    set(hi3,'cdata',seisd);
    ht=get(hax3,'title');
    P1=sum(seis1(:).^2);
    PD=sum(seisd(:).^2);
    pctd=round(1000*PD/P1)/10;
    set(ht,'string',{'Least-squares subtraction',['Difference power is ' num2str(pctd) '% of dataset #1']});
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