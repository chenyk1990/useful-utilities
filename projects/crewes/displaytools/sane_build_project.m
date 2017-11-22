function sane_build_project(segyfiles,segypath,segyrev,projfilename,projpath,projname,timeshifts,loadedflags,displayflags)
% sane_build ... builds a SANE project
%
% sane_build(segyfiles,segypath,segyrev,projfilename,projpath,projname,timeshifts,loadedflags,displayflags)
%
% segyfiles ... cell array of segy file names. Must end in .sgy. If a single file, then it can be
%               just a string.
% segypath ... cell array of paths for the segy files. If the files all have the same path then this
%               can be a single string.
% segyrev ... SEGY revision number for each file. Must be a velctor of length(segyfiles) or a single
%               value if they are all the same. An entry of [] causes the revision number in the
%               binary header to be used. An entry of 1 causes all fiels to be read with segy
%               revision 1 regadless of what is in the binary header. Same idea for 2.
% projfilename ... file name (string) for the finished project. A .mat will be appended automatically.
% projpath ... path (string) for the finished project
% projname ... name (string) of the project. Any existing file with this name will be deleted.
%          *********** default is projfilename without the .mat **********
% timeshifts ... array of datum shifts (in seconds or depth units) one for each dataset
%          *********** default is all zeros ***************
% loadedflags ... vector of same length as segyfiles containing either 1 or 0. 1 indicates the dataset
%               is to be loaded into memory upon opening the project, 0 indicates the file is to
%               remain on disk.
%          *********** default is all zeros *********
% displayflags ... vector of same length as segyfiles containing either 1 or 0. 1 indicates the
%               dataset is to be displayed in plotimage3D upon opening the project. 
%          *********** default is all zeros *********
% NOTE: if displayflags(k) is 1 then loadedflags(k) will be forced to 1.
%
%

if(iscell(segyfiles))
    nfiles=length(segyfiles);
else
    nfiles=1;
    if(~ischar(segyfiles))
        error('segyfiles must be either a string or a cell array');
    end
    segyfiles={segyfiles};
end

if(iscell(segypath))
    npaths=length(segypath);
else
    npaths=1;
    if(~ischar(segypath))
        error('segypath must be either a string or a cell array');
    end
    segypath={segypath};
end
if(nfiles>1 && npaths==1)
    paths2=cell(1,nfiles);
    for k=1:nfiles
        paths2{k}=segypath{1};
    end
    segypath=paths2;
end
if(~ischar(projfilename))
    error('projfilename must be a string');
else
    ind=strfind(projfilename,'.mat');
    if(isempty(ind))
        projfilename=[projfilename '.mat'];
    else
        if(ind(1)~=length(projfilename)-3)
            error('projfilename must end in .mat')
        end
    end
end
if(~ischar(projpath))
    error('projpath must be a string');
end

if(length(segyrev)==1)
    if(isempty(segyrev))
        segyrev=nan*ones(1,nfiles);
    else
        segyrev=segyrev*ones(1,nfiles);
    end
end
if(length(segyrev)~=nfiles)
    error('there must be a segyrev for each input');
end
segyrev=round(segyrev);
if(any(segyrev<0) || any(segyrev>2))
    error('Invalid value for segyrev');
end

%check for existance on input files
existflags=zeros(1,nfiles);
datanames=cell(1,nfiles);
for k=1:nfiles
    if(exist([segypath{k} segyfiles{k}],'file'));
        existflags(k)=1;
    end
    ind=strfind(segyfiles{k},'.sgy');
    if(~isempty(ind))
        datanames{k}=segyfiles{k}(1:ind(1)-1);
    else
        datanames{k}=segyfiles{k};
    end
end
if(any(existflags==0))
    ind=find(existflags==0);
    for k=ind
        disp([segypath{k} segyfiles{k} ' does not exist']);
    end
    error('One or more input files does not exist');
end
if(nargin<6)
    ind=strfind(projfilename,'.mat');
    projname=projfilename(1:ind(1)-1);
end
if(nargin<7)
    timeshifts=zeros(1,nfiles);
end
if(nargin<8)
    loadedflags=zeros(1,nfiles);
end
if(nargin<9)
    displayflags=zeros(1,nfiles);
end
if(length(timeshifts)~=nfiles)
    error('there must be one time shift for each input file');
end
if(length(loadedflags)~=nfiles)
    error('there must be one loadedflag for each input file');
end
if(length(displayflags)~=nfiles)
    error('there must be one displayflag for each input file');
end
for k=1:nfiles
    if(displayflags(k)==1)
        loadedflags(k)=1;
    end
end

if(exist([projpath projfilename],'file'))
    delete([projpath projfilename]);
    disp(['Existing file ' projpath projfilename ' deleted'])
end

% ok now make a project structure
proj.name=projname;
proj.filenames=segyfiles;
proj.projfilename=projfilename;
proj.projpath=projpath;
proj.paths=segypath;
proj.datanames=datanames;
proj.isloaded=loadedflags;
proj.isdisplayed=displayflags;
proj.xcoord=cell(1,nfiles);
proj.ycoord=cell(1,nfiles);
proj.tcoord=cell(1,nfiles);
proj.tshift=timeshifts;
proj.datasets=cell(1,nfiles);
proj.xcdp=cell(1,nfiles);
proj.ycdp=cell(1,nfiles);
proj.dx=ones(1,nfiles);
proj.dy=ones(1,nfiles);
proj.depth=zeros(1,nfiles);
proj.texthdr=cell(1,nfiles);
proj.texthdrfmt=cell(1,nfiles);
proj.segfmt=cell(1,nfiles);
proj.byteorder=cell(1,nfiles);
proj.binhdr=cell(1,nfiles);
proj.exthdr=cell(1,nfiles);
proj.tracehdr=cell(1,nfiles);
proj.bindef=cell(1,nfiles);
proj.trcdef=cell(1,nfiles);
proj.segyrev=segyrev;
proj.kxline=cell(1,nfiles);
proj.gui=cell(1,nfiles);
proj.rspath=[];
proj.wspath=[];
proj.rmpath=[];
proj.wmpath=[];
proj.pifigures=cell(1,nfiles);
proj.isdeleted=zeros(1,nfiles);
proj.deletedondisk=zeros(1,nfiles);
proj.saveneeded=zeros(1,nfiles);

%read in the datasets
datasets=cell(1,nfiles);
t0=clock;
for k=1:nfiles
    disp(['reading file #' int2str(k) ' as sgyrev ' int2str(segyrev(k))]); 
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
    [seis,sgrv,dt,proj.segfmt{k},proj.texthdrfmt{k},proj.byteorder{k},proj.texthdr{k},proj.binhdr{k},...
        proj.exthdr{k},proj.tracehdr{k},proj.bindef{k},proj.trcdef{k}] =readsegy([segypath{k} segyfiles{k}],[],segyrev(k),[],[],...
        [],[],[],[],[],1); %#ok<ASGLU>
    
    t=dt*(0:size(seis,1)-1)'+timeshifts(k);
    proj.tcoord{k}=t;
    dt=abs(t(2)-t(1));
    if(dt>.02)
        proj.depth(k)=1;
    end
    %we need InlineNum and XlineNum to be populated
%     y=double(proj.tracehdr{k}.InlineNum);
%     x=double(proj.tracehdr{k}.XlineNum);
%     if(sum(abs(x))==0)
%         disp(['dataset ' segypath{k} segyfiles{k}])
%         error('Xline header information is empty. Cannot proceed');
%     end
%     if(sum(abs(y))==0)
%         disp(['dataset ' segypath{k} segyfiles{k}])
%         error('Inline header information is empty. Cannot proceed');
%     end
%     cdpx=proj.tracehdr{k}.CdpX;
%     cdpy=proj.tracehdr{k}.CdpY;
    if(isfield(proj.tracehdr{k},'InlineNum')) %this is a field for segyrev1 but not rev0
        %we need InlineNum and XlineNum to be populated
        y=double(proj.tracehdr{k}.InlineNum);
        x=double(proj.tracehdr{k}.XlineNum);
        if(sum(abs(x))==0)
            %try the Kingdom header locations
            if(isfield(proj.tracehdr{k},'TrcNumCDP'))
                x=double(proj.tracehdr{k}.TrcNumCDP);
                y=double(proj.tracehdr{k}.SourcePoint);
            else
                x=0;
            end
            if(sum(abs(x)==0))
                msgbox('Xline/Inline header information is empty. Cannot proceed');
                waitsignaloff;
                return
            end
        end
    else
        %try the Kingdom header locations
        if(isfield(proj.tracehdr{k},'TrcNumCDP'))
            x=double(proj.tracehdr{k}.TrcNumCDP);
            y=double(proj.tracehdr{k}.SourcePoint);
        else
            x=0;
        end
        if(sum(abs(x)==0))
            msgbox('Xline/Inline header information is empty. Cannot proceed');
            waitsignaloff;
            return
        end
    end
    if(isfield(proj.tracehdr{k},'CdpX'))
        cdpx=proj.tracehdr{k}.CdpX;
        cdpy=proj.tracehdr{k}.CdpY;
    else
        cdpx=zeros(size(x));
        cdpy=zeros(size(y));
    end
    [datasets{k},proj.xcoord{k},proj.ycoord{k},xcdp,ycdp,proj.kxline{k}]=make3Dvol(seis,x,y,cdpx,cdpy);
    
    proj.dx(k)=max([abs(xcdp(2)-xcdp(1)), 1]);
    proj.dy(k)=max([abs(ycdp(2)-ycdp(1)), 1]);
    proj.xcdp{k}=xcdp;
    proj.ycdp{k}=ycdp;
    tnow=clock;
    timeused=etime(tnow,t0);
    disp(['Finished dataset #' int2str(k) ' of ' int2str(nfiles)]);
    disp(['Time used = ' num2str(timeused/60) ' minutes'])
end

mObj=matfile([projpath projfilename],'writable',true);
mObj.proj=proj;
mObj.datasets=datasets;
tnow=clock;
timeused=etime(tnow,t0);
disp('Project saved to disk');
disp(['Total project build time = ' num2str(timeused/60) ' minutes'])
