function sf=writesegy(sgyfile,trcdat,segyrev,sampint,segfmt,txtfmt, ...
                   bytord,txthdr,binhdr,exthdr,trchdr,bindef,trcdef,gui)
% function  writesegy(sgyfile,trcdat,segyrev,sampint,segfmt,txtfmt, ...
%                  bytord,txthdr,binhdr,exthdr,trchdr,bindef,trcdef,gui)
%
% Write Matlab structure, array and scalar data to a new big-endian SEG-Y file.
%
% Required inputs:
%   sgyfile   - filename (should end in .sgy or .segy)
%   trcdat    - 2-D matrix of samples indexed by (sample, trace)
%
% Optional inputs, [] is accepted:
%   segyrev   - segy revision number (0, 1 or 2)
%   sampint   - sample interval in seconds (default 1e-3)
%   segfmt    - output segfmt (default=5)
%               1=IBM floating point, 2=4-byte integer, 3=2-byte integer,
%               5=IEEE floating point,8=1-byte integer
%   txtfmt    - text file format (default='ascii')
%               'ascii', 'ebcdic'
%   bytord    - output byte order (default='b')
%               'b' = big-endian (SEG-Y standard), 'l' = little-endian
%   txthdr    - 2-D char array containing an ASCII SEG-Y textual file header. This
%               will be written to disk as EBCDIC if txtfmt='ebcdic'
%   binhdr    - binary header struct
%   exthdr    - 2-D char array containing ASCII SEG-Y extended textual file
%               header(s)
%   trchdr    - trace header struct (optional)
%   bindef    - 4 column binary header definition cell array such as provided by
%               @BinaryHeader/new; See uiSegyDefinition().
%               NOTE! writesegy requires the same custom bindef used with 
%               readsegy unless you modify the binhdr struct to match!
%   trcdef    - 5 column trace header definition cell array such as provided by
%               @BinaryHeader/new; See uiSegyDefinition()
%               NOTE! writesegy requires the same custom trcdef used with 
%               readsegy unless you modify the binhdr struct to match!
%   gui       - 0 (no progress bar or prompts; guess and go...), 1 (text
%               progress bar and prompts, [] (default; gui progress bar and 
%               prompts), h=figure handle (same as [], but an attempt is
%               made to center gui popups on the figure represented by h)
%
% Examples:
%   writesegy('test_ieee.sgy',dataout,1,sampint)
%   writesegy(sgyfile, dataout, segyrev, sampint, segfmt, txtfmt, ...
%                     byteorder, txthdr, binhdr, exthdr, trchdr)
%   writesegy('test_ieee.sgy',dataout,[],sampint,[],[],[],txthdr,binhdr,[],trchdr)
%
% NOTE: writesegy is a wrapper for @SegyFile
%
% Authors: Kevin Hall, 2017
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

% Make certain we have between three and eleven input arguments
%narginchk(2,14)
%                        1       2        3       4       5       6
% function  writesegy(sgyfile, trcdat, sampint, segfmt, txtfmt, bytord,

narginchk(2,14)

% Check inputs
if nargin <3 || isempty(segyrev)
    segyrev=1; 
end
if nargin < 4 || isempty(sampint)
    sampint = 1e-3; %let's just pretend it's 1 ms
end
if nargin <5 || isempty(segfmt)
    segfmt = 5; %4-byte IEEE float
end
if nargin < 6 || isempty(txtfmt)
    txtfmt = 'ascii';
end
if nargin < 7 || isempty(bytord)
    bytord = 'b';
end
if nargin < 8
    txthdr = [];
end
if nargin < 9
    binhdr = [];
end
if nargin < 10
    exthdr = [];
end
if nargin < 11
    trchdr = [];
end
if nargin < 12
    bindef = [];
end
if nargin < 13
    trcdef = [];
end
if nargin < 14
    gui=[];
end
    
[nsamps,numtrc] = size(trcdat);

sf = SegyFile(sgyfile,'w',segyrev,sampint,nsamps,...
    segfmt,txtfmt,bytord,bindef,trcdef,gui);
                    
%Binary Header sanity checks
if segyrev==0 && segfmt>4
    if isempty(gui) || isa(gui, 'handle')
        mm_warndlg(['writesegy: Trace data Format Code ' num2str(segfmt) ...
            ' is not supported in SEG-Y Rev 0. ' ...
            'Continuing...'],'Warning!',gui);
    elseif gui
        warning(['writesegy: Trace data Format Code ' num2str(segfmt) ...
            ' is not supported in SEG-Y Rev 0. ' ...
            'Continuing...']);
    end
elseif segyrev==1 && ((segfmt>5 && segfmt<8) || segfmt >8)
    if isempty(gui) || isa(gui, 'handle')
        mm_warndlg(['writesegy: Trace data Format Code ' num2str(segfmt) ...
            ' is not supported in SEG-Y Rev 1. ' ...
            'Continuing...'],'Warning!',gui);
    elseif gui
        warning(['writesegy: Trace data Format Code ' num2str(segfmt) ...
            ' is not supported in SEG-Y Rev 0. ' ...
            'Continuing...']);
    end
elseif segyrev==2
    error('This code does not yet support writing SEG-Y Rev 2 Files');
end

if isempty(binhdr)
    binhdr = sf.BinHdr.new(nsamps,sampint,segfmt); %rev1 by default
else
    binhdr.(sf.BinHdr.byte2word(16)) = sampint*1e6; %microseconds
    binhdr.(sf.BinHdr.byte2word(18)) = sampint*1e6; %microseconds
    binhdr.(sf.BinHdr.byte2word(20)) = nsamps;  
    binhdr.(sf.BinHdr.byte2word(22)) = nsamps;    
    binhdr.(sf.BinHdr.byte2word(24)) = segfmt;
end

if segyrev==0 && segfmt==5 %update segy revision major #. Warn user??
    binhdr.(sf.BinHdr.byte2word(300)) = 1;  
end

if isempty(txthdr)
    txthdr = sf.TxtHdr.new();
end

if ~isempty(exthdr)
    error('Extended Textual File Headers are not yet supported')
end

if ~isempty(trcdef)
    sf.Trc.HdrDef=trcdef; %custom trace header def
end

if isempty(trchdr)
    trchdr = sf.Trc.new(numtrc,nsamps,sampint,'headers');
else
    trchdr.(sf.Trc.byte2word(114))(:) = nsamps;
    trchdr.(sf.Trc.byte2word(116))(:) = sampint*1e6;
end

%     if updatetrchdr  %update to segyrev 1
%         ntrchdr = sf.Trc.new(numtrc,numsamp,sampint,'headers'); %rev1 by default
%         f = fieldnames(ntrchdr);        
%         for ii=1:length(f)
%             if isfield(trchdr,f(ii)) && isfield(ntrchdr,f(ii))
%                 ntrchdr.(f{ii}) = trchdr.(f{ii});
%             end
%         end
%         trchdr=ntrchdr;
%     end


%% Do the dirty deed
sf.TxtHdr.write(txthdr);
sf.BinHdr.write(binhdr);
sf.ExtTxtHdr.write(exthdr);
sf.Trc.write(trcdat,trchdr);

%end function writesegy