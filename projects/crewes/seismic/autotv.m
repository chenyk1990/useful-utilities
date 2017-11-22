function [autos,tlags]=autotv(seis,t,t1s,t2s,tmlags,pflag,dname,x,xname)
% AUTOTV: compute and optionally display time variant autocorrelations
%
% [autos,tlags]=autotv(seis,t,t1s,t2s,tmlags,pflag,dname,x,xname)
%
% For each specified window, the autocorrelation functions are computed and
% displayed for each trace in the seismic matrix. If plotting is requested,
% then it is best to have fewer than about 5 windows (3 is optimal). In each
% window, the data is multiplied by an mwindow (boxcar with raised cosine
% tapers) to minimize edge effects. This is a useful tool to study the
% multiple contamination of s seismic dataset.
%
% seis ... input seismic section or gather
% t ... time coordinate for seis
% NOTE: length(t) must equal size(seis,1)
% t1s ... start times of the autocorrelation windows
% t2s ... end times of the windows
% tmlags ... max lag in each window
% Note: t1s, t2s, and tmlags must all be vectors of the same size. However
% tmlags of length 1 is allowed meaning each window has the same lags.
% pflag ... plotting flag, 1 to plot, 0 to not plot
% ******* default = 1 ************
% dname ... string for a title on the plot
% ******* default = '' ****************
% x ... x coordinate
% ******* default = trace number *********
% xname ... string naming the x coordinate 
% ******* default = 'trace number' *********
% 
% G.F. Margrave, Devon, 2016
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

nwins=length(t1s);
if(nwins~=length(t2s))
    error('t1s and t2s must be same size vectors');
end

if(length(tmlags)==1)
    tmlags=tmlags*ones(size(t1s));
end

if(nwins~=length(tmlags))
    error('tmlags must be the same size as t1s');
end

if(nargin<6)
    pflag=1;
end
if(nargin<7)
    dname='';
end
if(nargin<8)
    x=1:size(seis,2);
end
if(nargin<9)
    xname='trace number';
end
    

autos=cell(1,nwins);
tlags=autos;
dt=t(2)-t(1);
nx=size(seis,2);

for k=1:nwins
    tmlag=round(tmlags(k)/dt)*dt;
    tau=-tmlag:dt:tmlag;
    n=tmlag/dt;
    ind=near(t,t1s(k),t2s(k));
    mw=mwindow(length(ind));
    otto=zeros(length(tau),nx);
    for kk=1:nx
        otto(:,kk)=auto2(mw.*seis(ind,kk),n);
    end
    autos{k}=otto;
    tlags{k}=tau;
end

if(pflag)
    figure
    axht=.8/nwins;
    axwid=.8;
    xnot=.1;
    sep=.02;
    ynow=.95;
    for k=1:nwins
        ynow=ynow-axht;
        axes('position',[xnot,ynow,axwid,axht-sep])
        imagesc(x,tlags{k},autos{k},[-.8,1]);
        text(x(1),.5*tmlags(k),['Autos from ' num2str(t1s(k)) 's to ' num2str(t2s(k)) 's'])
        ylabel('lag time (s)');
        if(k<nwins)
            set(gca,'xticklabel',[])
        else
            xlabel(xname);
        end
        if(k==1)
            ht=title(dname);
            set(ht,'interpreter','none');
        end
        grid
        colorbar
    end
    if(~iscell(dname))
        set(gcf,'name',['autotv: ' dname])
    else
        set(gcf,'name',['autotv: ' dname{1}])
    end
    prepfiga
end
        
