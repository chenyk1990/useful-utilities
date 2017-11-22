function file_snippet(infile,outfile,nbytes)
%function file_snippet(infile,outfile,nbytes)
%
% Reads nbytes from infile and writes them to outfile
%

infid = fopen(infile,'r');
if infid < 1
    error('fopen for file read failed')
end

if exist(outfile,'file')
    error([outfile ' already exists, refusing to overwrite'])
end
outfid = fopen(outfile,'w');
if outfid < 1
    error('fopen for file write failed')
end

[data,count] = fread(infid, nbytes, 'uint8=>uint8');
if ~isequal(count,nbytes)
    error('fread failed')
end

count = fwrite(outfid, data, 'uint8');
if ~isequal(count,nbytes)
    error('fwrite failed')    
end

fclose(infid);
fclose(outfid);

end %end function






