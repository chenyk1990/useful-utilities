function pep=penpred(s,r,w,tw)
% PENPRED: portion of energy predicted
% pep=penpred(s,r,w,tw)
% 
% PEP = portion of energy predicted
%
% Method: Given a seismic trace, a reflectivity, and a wavelet, the wavelet and reflectivity are
% convolved together and the result is compared to the seimic trace. The PEP is defined as 
% PEP= 1-Energy(s-conv(r,w))/Energy(s).
% Where Energy(x)=sum(x.^2). Notice that if s is exactly matched by conv(r,w) then PEP is 1. Also
% notice that the formula assumed that s and conv(r,w) are balanced in overall amplitude. This
% function takes care of that by determining the best(least-squares) scalar to be applied to
% conv(r,w) to minimize Energy(s-conv(r,w)). (See lsqsubtract for more on this).
% 
% s ... seismic trace
% r ... reflectivity
% NOTE: s and r must be the same length
% w ... wavelet
% tw ... time coordinate for wavelet.(Needed to deterine time zero of the
%       wavelet)
% NOTE: If w and tw are not provided, then r is used directly as if it is already convolved.
%
%
% G.F. Margrave, Devon Energy, 2016
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

if(length(s)~=length(r))
    error(' s and r must be the same length');
end

if(nargin>2)
    if(length(w)~=length(tw))
        error(' w and tw must be the same length');
    end
    izero=near(tw,0);
    s2=convz(r,w,izero(1),length(r),0);
else
    s2=r;
end

%s2b=balans(s2,s);
%s2b=s2;
[~,a]=lsqsubtract(s,s2);
s2b=a*s2;

pep=1-sum((s-s2b).^2)/sum(s.^2);





