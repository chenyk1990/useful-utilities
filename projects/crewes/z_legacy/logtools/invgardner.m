function v=invgardner(rho,a,m)
% INVGARDNER - deduce velocity from density given Gardner parameters
% 
% v=invgardner(rho,a,m)
%
% INVGARDNER estimates a velocity vector given a density vector
% and values for the empirical parameters a and m as defined the Gardner relation
%	rho=a*v.^m (a times v to the mth power)
% Thus v = (rho/a).^(1/m)
%
% 	rho= density vector
%	a= scalar multplier
%		*********** default is .31 **********
%	m= scalar exponent
%		*********** default is .25 **********
%	v= a vector of velocities computed from v = (rho/a).^(1/m)
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

if(nargin<3)
	m=.25;
end
if(nargin<2)
	a=.31;
end

v= (rho/a).^(1/m);