function ss=stretcht(s,t,delt,tdel)
% STRETCHT: stretch a time series from t to t+delt
%
% ss=stretcht(s,t,delt,tdel)
%
% s... input trace
% t... time coordinate for s
% delt ... vector of shifts, same size as t. A positive value moves samples
%      to greater times, a negative value to earlier times.
% tdel ... times at which delt is specified
% ************ defaults to t *********
%
% ss... stretched trace, ss(t) corresponds to s(t-delt)
%
% by: G.F. Margrave, Devon Canada, 2016
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
    tdel=t;
end
if(length(tdel)~=length(delt))
    error('tdel and delt must have the same length');
end
if(length(t)~=length(s))
    error('s and t must have the same length');
end
s=s(:);
t=t(:);
delt=delt(:);
tdel=tdel(:);
if(length(tdel)~=length(t))
    %need to interpolate the delays
    if(tdel(1)>t(1))
        tdel=[t(1);tdel];
        delt=[delt(1);delt];
    end
    if(tdel(end)<t(end))
        tdel=[tdel;t(end)];
        delt=[delt;delt(end)];
    end
   delt=interp1(tdel,delt,t); 
end
        
ss=interpbl(t,s,t-delt);