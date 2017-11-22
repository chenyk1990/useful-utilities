function seisplota(seis,t,x,ampinfo)
% SEISPLOTA ... creat a sesimic image plot in the current axes
%
% seisplota(seis,t,x,ampinfo)
%
% This function provides a similar display to that of plotimage but in the current axes instead
% of a new window and  without the extra features (such as picking). This is suitable when you
% want a quick seismic display in the current axes. Use seisplot for a similar display in a new
% figure window. Seisplot also gives interactive ability to choose the clip level. You might
% also consider seisplotfk and seisplot2. If the default display is too dark, follow this command
% with brighten(.5).  Similarly if it is too light use brighten(-.5).
%
% seis ... seismic matrix (gather).One trace per column
% t ... time coordinate vector. length(t) must equal size(seis,1)
% x ... space coordinate vector. length(x) must equal size(seis,2)
% ampinfo ... max, min, mean, and standard deviation of seismic data used to determine clipping. Specify as
%       a four element vector [max min mean stddev]. This is useful when several different displays need
%       to have the same realtive amplitude.
% ************ default is to compute these from the input data ************
%
% NOTE: The image is plotted with Matlab's imagesc. This function only annotates the axes
% precisely correctly if the x and t vectors are regularly sampled. This is usually the case
% with t but less ofen so with x. For precisly annotated tick marks when x is not regular, the
% only current option is to uses plotseis or plotseismic which both plot wiggle traces not
% images.
%
% 
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

if(strcmp(action,'init'))

if(nargin<4)
    ampinfo=[max(seis(:)) min(seis(:)) mean(seis(:)) std(seis(:))];
end

if(length(t)~=size(seis,1))
    error('time coordinate vector does not match seismic');
end
if(length(x)~=size(seis,2))
    error('space coordinate vector does not match seismic');
end

[clips,clipstr,clip,iclip]=getclips(ampinfo);
if(isempty(clip))
    am=ampinfo(1);
    clim=[-am am];
else
    am=ampinfo(3);
    sigma=ampinfo(4);
    clim=[am-clip*sigma am+clip*sigma];
end

hi=imagesc(x,t,seis,clim);colormap(seisclrs);
% brighten(.5)
grid
hcm=uicontextmenu;
uimenu(hcm,'label','Time-variant spectrum','callback',@showtvspectrum);
set(hi,'uicontextmenu',hcm);
if(max(x)<5000)
    xlabel('distance (m)')
else
    xlabel('distance (ft)')
end
if(max(t)<10)
    ylabel('time (s)')
elseif(max(t)<5000)
    ylabel('depth (m)')
else
    ylabel('depth (ft)')
end


ht=.05;
wid=.05;
%wid=.045;sep=.005;
pos=get(gca,'position');
xnow=pos(1)+pos(3);
ynow=pos(2)+pos(4)-2*ht;
uicontrol(gcf,'style','popupmenu','string',clipstr,'tag','clip','units','normalized',...
    'position',[xnow,ynow,wid,ht],'callback','seisplota(''clip'')','value',iclip,...
    'userdata',{clips,ampinfo,gca},'visible','on','tooltipstring',...
    'clip level is the number of standard deviations from the mean at which amplitudes are clipped')

elseif(strcmp(action,'clip'))
    hclip=gcbo;
    iclip=get(hclip,'value');
    udat=get(hclip,'userdata');
    clips=udat{1};
    ampinfo=udat{2};
    haxe=udat{3};
    amax=ampinfo(1);
    am=ampinfo(3);
    sigma=ampinfo(4);
    if(iclip==1)
        clim=[-amax amax];
    else
        clip=clips(iclip);
        clim=[am-clip*sigma,am+clip*sigma];
    end
    set(haxe,'clim',clim);
end

end

function [clips,clipstr,clip,iclip]=getclips(ampinfo)
%
% 
% clips ... determined clip levels
% clipstr ... cell array of strings for each clip level for use in popup menu
% clip ... starting clip level
% iclip ... index into clips where clip is found
%

sigma=ampinfo(4);
amin=ampinfo(2);
amax=ampinfo(1);
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
