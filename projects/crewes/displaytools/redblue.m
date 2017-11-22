function a=redblue(m)
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
n=3;
if( nargin < 1) 
	m=size(get(gcf,'colormap'),1); 
end
imid=m/2;
ind =linspace(1,m,m);
%r=((cos(pi*(ind-1)/(m-1))+1)/2)';
r=.8*mwhalf(m,55)+.2;
g=r(imid)*(.3+.7*((cos(-pi*(ind-imid)/(imid))+1)/2)');
b = r(imid)*(.5+.5*((cos(-pi*(ind-imid)/(imid))+1)/2)');
ind=1:imid;
tmp=g(ind);
g(ind)=b(ind);
b(ind)=tmp;
a = [r.^n g.^n b.^n];