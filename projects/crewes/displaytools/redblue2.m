function a=redblue2(m)
%
% a=alpine(m)
%
% Returns the CREWES alpine colormap with m colors
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

%

if( nargin < 1) 
	m=size(get(gcf,'colormap'),1); 
end
imid=floor(m/2);
m=2*imid+1;
% r=.667*[ones(1,iend) linspace(1,0,m-2*iend) zeros(1,iend)]';
r0=.3;r1=.8;
b0=.3;b1=.8;
g0=0;g1=.8;
r=[linspace(r0,r1,m-imid) zeros(1,imid)]';
%g=.333*(.5*cos(-pi*(x-imid)/(1-imid))+.5)';
%g=(.5*cos(-pi*(x-imid)/(1-imid))+.5)';
ge=linspace(g1,g0,imid+1);
g=[linspace(g0,g1,m-imid) ge(2:end)]';
%b = .667*[zeros(1,iend) linspace(0,1,m-2*iend) ones(1,iend)]';
b=[zeros(1,imid) linspace(b1,b0,m-imid)]';
a = [r g.^2 b];
a = brighten(a,.5);