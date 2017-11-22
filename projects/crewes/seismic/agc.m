function [trout,envsm,ma]=agc(trin,sampint,op_length,trip,profileflag)
% AGC: This is just a wrapper around AEC (automatic envelope correction)
%
% Syntax
%  trout=aec(trin);
%  trout=aec(trin,sampint,op_length,trip)
%  trout=aec(trin,sampint,op_length,trip,profileflag)
%
% Description 
%   AEC performs an automatic amplitude adjustment.
%
% Method
%   1) Compute Hilbert envelope of the input trace TRIN
%   2) Convolve envelope with triangular smoother of half-length
%      OP_LENGTH
%   3) Divide input trace by smoothed envelope
%   4) Balance the output to a maximum of 1.0
%
% Inputs
%   trin      = input trace or gather of traces.
%   t or sampint   = sample interval for trin 
%               For backwards compatibility, if sampint is supplied as a time coordinate 
%               vector, the difference of the first two elements is used as
%               the sample interval.
%   **********  Default is 0.001 seconds.  (1 millisecond) ********
%   op_length = half-length of triangular smoother in seconds
%   **********  default is 1/8th of the trace length *******
%               ***** must be less than half the trace length *****
%   trip      = front end time before which the smoothed envelope is
%               set to a constant
%   **********  default is op_length/10 *******
%  profileflag= 0 or 1, this applies only if the input is a profile of
%               traces. In that case, if profileflag is 0, then each trace
%               is divided by its own envelope, and this is called single
%               trace more. If profileflag is 1, then the hilbert
%               envelopes of all the traces are stacked and the result is
%               smoothed to get a single smooth envelope. Each trace is
%               then divided by this envelope and normalized.
%   ***********  default is 0 ***********
% Outputs
%   trout     = output trace or gather of traces
%   envsm     = smoothed hilbert envelope of trace or traces
%   ma        = maximum absolut value of the corretect trace befor
%               normalization to 1.
% NOTES:
%   1) To remove trace normalization to 1, trout_unnorm=trout*ma; (single
%           trace) or trout_unnorm(:,k)=trout(:,k)*ma(a); (multitrace)
%   2) To recover the input trace: trin_recovered=trout*ma.*envsm; (single
%           trace) or trin_recovered(:,k)=trout(:,k)*ma.*envsm(:,k) (multitrace)
%
% by G.F. Margrave, May 1991-2016 (Henry Bland and Kevin Hall too)
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
   
% set defaults
if nargin < 2 || isempty(sampint)
    sampint = 0.001;   % sample interval in seconds
else
    % Backwards compatibility with time coordinate vector
    if length(sampint) > 1
        sampint = sampint(2) - sampint(1);
    end
end
if nargin < 3 || isempty(op_length)
    op_length = sampint*length(trin)/8;
end
if nargin < 4 || isempty(trip);
    trip=op_length/10.;
end
if(nargin<5)
    profileflag=0;
end

[trout,envsm,ma]=aec(trin,sampint,op_length,trip,profileflag);