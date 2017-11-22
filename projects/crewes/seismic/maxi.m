function [ccmax,taumax]=maxi(cc,tau,taufrac)
% MAXI: find the maximum of a time series interpolated to the nearest 1/10 th sample
%
% [ccmax,taumax]=maxi(cc,tau,taufrac)
%
% MAXI first finds the maximum value on a sample. It then takes 10 samples to each side of this
% maximum and does a band-limited interpolation (using interpbl) to 1/10th of the sample size and
% finds the maximum of that interpolated segment.
%
% cc ... input time series (usually a crosscorrelation but can be anything)
% tau ... time coordinate vector for cc (same size as cc)
% taufrac ... fractional sample size of the interpolated result
% ************ default = 10 ***********
% The default means that the result is interpolated to 1/10 th of the input sample size
% flag2 ... if 1, then the actual maximum is determined. If 2, then the max(abs()) is determined
% ************ default =1 ***************
%
% ccmax ... interpolated maximum of cc
% taumax ... time at which the interpolated maximum occurs
%
% G.F. Margrave, Devon Energy, 2017
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

if(nargin<4)
    flag2=1;
end
if(nargin<3)
    taufrac=10;
end

if(flag2==1)
[tmp,itmp]=max(cc); %#ok<*ASGLU>
else
[tmp,itmp]=max(abs(cc));   
end
nt=length(tau);
i1=max([1 itmp-10]);
i2=min([nt itmp+10]);
dt=abs(tau(2)-tau(1));
tint=tau(i1):dt/taufrac:tau(i2);
cci=interpbl(tau(i1:i2),cc(i1:i2),tint);

if(flag2==1)
    [ccmax,imax]=max(cci);
else
    [ccmax,imax]=max(abs(cci));
end
taumax=tint(imax);