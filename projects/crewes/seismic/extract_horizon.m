function ahor=extract_horizon(seis3D,t,thor,method,twidth)
%
% seis3D ... 3D seismic volume
% t ... time coordinate for seis3D
% NOTE: length(t) must equal size(seis3D,1)
% thor ... 2D matrix of horizon times
% NOTE: size(thor) must equal size(squeeze(seis3D(1,:,:))). It is a
%   confusing point that this means that the x coordinate (xline) of thor
%   must be the row coordinate not the column coordainte.
% method ... string. must be one of 'amp','aveamp','rmsamp'
% twidth ... length 2 vector giving the width of the extraction window
%       about the times defined in thor. That is, for 'aveamp' and 'rmsamp'
%       the averages will be take over the window thor-twidth(1) to
%       thor+twidth(2). Twidth may be specified as a single value if the up
%       and down widths are the same. Has no effect for 'amp'.
% NOTE: twidth must be positive values in units of seconds
%
% thor is a matrix of tsize equal to the spatial dimensions of seis3D. thor
% may contain nans if the hroizon is not defined in places.
% 

if(nargin<5)
    twidth=[0 0];
end

[nt,nx,ny]=size(seis3D);
if(length(t)~=nt)
    error('t and seis3D have incompatible sizes');
end
[nr,nc]=size(thor);
if(nr~=nx)
    error('horizon has incorrect xline dimension');
end
if(nc~=ny)
    error('horizon has incorrect inline dimension');
end
if(~strcmp(method,'amp')&&~strcmp(method,'aveamp')&&~strcmp(method,'rmsamp'))
   error('unrecognized method') 
end
if(length(twidth)==1)
    twidth=[twidth twidth];
end
if(length(twidth)>2)
    error('twidth must contain only two values');
end
if(any(twidth)<0)
    error('twidth must have positive entries');
end
tmin=min(t);
tmax=max(t);
if(any(thor-twidth(1))<tmin)
    error('horizon window falls outsize seismic time range');
end
if(any(thor+twidth(2))>tmax)
    error('horizon window falls outsize seismic time range');
end

ahor=zeros(nx,ny,'like',seis3D);
dt=abs(t(2)-t(1));
%loop over xlines
for k=1:nx
    switch method
        case 'amp'
            times=thor(k,:);
            ntimes=round((times-t(1))/dt)+1;
            for kk=1:ny
                ahor(k,kk)=seis3D(ntimes(kk),k,kk);
            end
        case 'aveamp'
            t1=thor(k,:)-twidth(1);
            n1=round((t1-t(1))/dt)+1;
            t2=thor(k,:)+twidth(1);
            n2=round((t2-t(1))/dt)+1;
            for kk=1:ny
                ahor(k,kk)=mean(seis3D(n1(kk):n2(kk),k,kk));
            end
        case 'rmsamp'
            t1=thor(k,:)-twidth(1);
            n1=round((t1-t(1))/dt)+1;
            t2=thor(k,:)+twidth(1);
            n2=round((t2-t(1))/dt)+1;
            for kk=1:ny
                nn=length(n1(kk):n2(kk));
                ahor(k,kk)=sqrt(sum(seis3D(n1(kk):n2(kk),k,kk).^2)/nn);
            end
    end
end
