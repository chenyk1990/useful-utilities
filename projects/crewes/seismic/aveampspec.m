function [A,f]=aveampspec(seis,t,t1,t2)
% AVEAMPSPEC ... compute the average amplitude spectrum of a seismic gather
%
% [A,f]=aveampspec(seis,t,t1,t2)
% 
% A default raised-cosine window (see mwindow) is used.
%
% seis ... the seismic gather (matrix)
% t ... time coordinate for seis
% NOTE: length(t) must equal size(seis,1)
% t1 ... start time for spectrum. Must be a scalar
% t2 end time for spectrum. Must be a scalar
%
% G.F. Margrave, Devon Canada, July 2017
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

if(t2<t1)
    error('t2 must be greater than t1');
end

if(length(t)~=size(seis,1))
    error('t and seis are incompatible');
end

nt=length(t);
ntpad=2*2^nextpow2(nt);

ind=near(t,t1,t2);
mw=mwindow(length(ind));
nx=size(seis,2);
for k=1:nx
    [S,f]=fftrl(seis(ind,k).*mw,t(ind),10,ntpad);
    if(k==1)
        A=abs(S);
    else
        A=A+abs(S);
    end
    
end

A=A/nx;