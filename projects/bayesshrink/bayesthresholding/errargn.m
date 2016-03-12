function err = errargn(ndfct,nbargin,argin,nbargout,argout)
%ERRARGN Check function arguments number.
%   ERR = ERRARGN('function',NUMARGIN,ARGIN,NUMARGOUT,ARGOUT)
%   is equal to 1 if either the number of input
%   ARGIN or output (ARGOUT) arguments of the specified
%   function do not belong to the vector of allowed values
%   (NUMARGIN and NUMARGOUT, respectively). 
%   Otherwise ERR = 0.
%
%   If ERR = 1, ERRARGN displays an error message in the
%   command window. The header of this error message contains
%   the string 'function'.
%
%   See also ERRARGT.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 01-May-96.
%   Last Revision: 29-Jun-1999.
%   Copyright 1995-2001 The MathWorks, Inc.
% $Revision: 1.9 $

% Special case:
% -------------
%  If ARGIN is not a numeric array, the number of input arguments
%  is not controlled. The same holds for ARGOUT.
%  example:
%    err = errargn('function',[0:2],'var',[1:4],2);
%    returns err = 0;

err = (isnumeric(argin)  & isempty(find(argin==nbargin))) | ...
      (isnumeric(argout) & isempty(find(argout==nbargout)));
if err , errargt(ndfct,'invalid number of arguments','msg'); end
