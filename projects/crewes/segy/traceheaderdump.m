function [dump,words,inotempty]=traceheaderdump(trchdr)
% TRACEHEADERDUMP: does a complete dump of the SEGY trace headers from a SegyFile object
%
% [dump,words,inotempty]=traceheaderdump(trchdr);
%      or
% traceheaderdump(trchdr);
%
% TRACEHEADERDUMP examines all fields in the traceheader structure and prints on the screen the names
% of those fields that are not empty and their SEGY word number. There are currently 91 possible words
% in the SEGY trace header according to SEGY Rev 1. Because these are the trace headers, each word has
% a possibly unique value for each trace; however, it is often true that these values are either
% constant or zero. Usually, trace header words that are entirely zero are not useful. This function
% makes it easy to find the non-zero entries.
%
% trchdr ... the trace header structure read by a SegyFile object.
% dump ... nwords by ntraces array of header values. nwords is the number
%   of defined header words that can be stored in a trace header according to SEGY Rev 1.  Currently
%   nwords=91.
% words ... cell array of all 91 defined header words
% inotempty ... vector of indices into words of the non-empty parts of the trace headers
% 
% To see a list of header words that are non-empty and hence can be useful, do:
%
% sf=SegyFile([path,fname],'r'); %create a SegyFile object, open for read.
% [trcdat,trchdr]=sf.Trc.read; %read traces and trace headers
% [dump,words,inotempty]=traceheaderdump(trchdr); %dump the trace headers
% words(inotempty) %note round brackets
%
% The last line here is not really necessary since traceheaderdump nows prints this information on
% the screen.
%
% You may wish to use traceheaderdump_g to see graphically what is in each nonzero header. 
%
% To retrieve the kth header value, if you ran traceheaderdump with returns, use:
% value = dump(k,:);
%
% Or if you ran it without returns use:
% value=double(getfield(trchdr,words{(k)}));
% note carefully the round and curly brackets in this statement. k must be
% an integer between 1 and length(inotempty). If you omit the double command then the header value
% will likely be of class int32.
% 
% by: G.F. Margrave, Devon Canada, 2017
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
words=fieldnames(trchdr);

ntraces=length(getfield(trchdr,words{1})); %#ok<*GFLD>
nwords=length(words);

dump=zeros(nwords,ntraces);

for k=1:nwords 
 dump(k,:)=getfield(trchdr,words{k});
end

test=sum(abs(dump),2);
inotempty=find(test~=0);
disp('Non-empty headers are as follows')
disp('Header word       ... Header word number')
ncols=17;
for k=1:length(inotempty)
   nblanks=ncols-length(words{inotempty(k)});
   disp([words{inotempty(k)} char(32*ones(1,nblanks)) ' ... ' char(32*ones(1,5)) int2str(inotempty(k))]) 
end

if(nargout==0)
clear dump words inotempty
end