function sanev(action)

if(nargin<1)
    action='init';
end

if(strcmp(action,'init'))
    hsanfig=figure;
    ssize=get(0,'screensize');
    width=floor(ssize(3)*.25);
    height=floor(ssize(4)*.25);
    xnot=(ssize(3)-width)*.5;
    ynot=(ssize(4)-height)*.5;
    set(hsanfig,'position',[xnot,ynot,width,height]);
    set(hsanfig,'menubar','none','toolbar','none','numbertitle','off','name','SANEV main window',...
        'nextplot','new');
    
    hfile=uimenu(hsanfig,'label','File');
    uimenu(hfile,'label','Load existing project');
    hread=uimenu(hfile,'Label','Read','callback','sanev(''readproj'');');
    uimenu(hread,'label','*.mat file','callback','sanev(''readmat'');');
    uimenu(hread,'label','*.sgy file','callback','sanev(''readsegy'');');
    uimenu(hfile,'label','Save project','callback','sanev(''saveproj'');');
    hwrite=uimenu(hfile,'label','Write');
    uimenu(hwrite,'label','*.mat file','callback','sanev(''writemat'');');
    uimenu(hwrite,'label','*.sgy file','callback','sanev(''writesegy'');');
    
elseif(strcmp(action,'readsegy'))
    %read in a segy dataset
    
elseif(strcmp(action,'readmat'))
    %read in a *.mat file
    [fname,path]=uigetfile('*.mat','Choose the .mat file to import');
    if(fname==0)
        return
    end
    m=matfile([path fname]);
    varnames=fieldnames(m);
    varnames(1)=[];
    varsizes=cell(size(varnames));
    threed=zeros(size(varsizes));
    %find any 3D matrices. only one is allowed
    for k=1:length(varnames)
        varsizes{k}=size(m,varnames{k});
        if(length(varsizes{k})==3)
            threed(k)=1;%points to 3D matrices
        end
    end
    if(sum(threed)>1)
        msgbox('Dataset contains more than 1 3D matrix. Unable to proceed.','Sorry!');
        return
    end
    %look for t,iline and xline
    i3d=find(threed==1);
    sz3d=size(m,varnames{i3d});
    nt=sz3d(1);%time is always the first dimension
    nx=sz3d(2);%xline is always the second dimension
    ny=sz3d(3);%inline is always the third dimension
    it=zeros(size(threed));
    ix=it;
    iy=it;
    %find things that are the size of nt, nx, and ny
    itchoice=0;
    ixchoice=0;
    iychoice=0;
    inamechoice=0;
    for k=1:length(varnames)
        if(k~=i3d)
            szk=varsizes{k};
            if((min(szk)==1)&&(max(szk)>1))
                %ok its a vector
                n=max(szk);
                if(n==nt)
                    %mark as possible time coordinate
                    it(k)=1;
                    if(strcmp(varnames{k},'t'))
                        itchoice=k;
                    end
                end
                if(n==nx)
                    %mark as possible xline
                    ix(k)=1;
                    if(strcmp(varnames{k},'xline'))
                        ixchoice=k;
                    end
                end
                if(n==ny)
                    %mark as possible iline
                    iy(k)=1;
                    if(strcmp(varnames{k},'iline')||strcmp(varnames{k},'inline')||strcmp(varnames{k},'yline'))
                        iychoice=k;
                    end
                end
            end
            if(contains(varnames{k},'dname'))
                inamechoice=k;
            end
        end
    end
    %ok, the best case is that it, ix, and iin all sum to 1 meaning there is only 1 possible
    %coordinate vector for each. If any one sums to zero, then we cannot continue. If any one sums
    %to greater than 1 then we have ambiguity that must be resolved.
    failmsg='';
    if(sum(it)==0)
        failmsg={failmsg 'Dataset contains no time coordinate vector. '};
    end
    if(sum(ix)==0)
        failmsg={failmsg 'Dataset contains no xline coordinate vector. '};
    end
    if(sum(iy)==0)
        failmsg={failmsg 'Dataset contains no inline coordinate vector. '};
    end
    if(~isempty(failmsg))
        msgbox(failmsg,'Sorry, dataset is not compaible with SANEV.')
        return
    end
    if(itchoice==0)
        ind=find(it==1);
        itchoice=it(ind(1));
    end
    if(ixchoice==0)
        ind=find(ix==1);
        ixchoice=ix(ind(1));
    end
    if(iychoice==0)
        ind=find(iy==1);
        iychoice=iy(ind(1));
    end
    ambig=[0 0 0];
    if(sum(it)>1)
        ambig(1)=1;
    end
    if(sum(ix)>1)
        ambig(2)=1;
    end
    if(sum(iy)>1)
        ambig(3)=1;
    end
    if(sum(ambig)>1)
        % put up dialog to resolve ambiguity
        [itchoice,ixchoice,iychoice]=ambigdialog(ambig,it,ix,iy,varnames,itchoice,ixchoice,iychoice);
    end
    seis=getfield(m,varnames{i3d}); %#ok<*GFLD>
    t=getfield(m,varnames{itchoice});
    xline=getfield(m,varnames{ixchoice});
    iline=getfield(m,varnames{iychoice});
    dname='';
    if(inamechoice>0)
        dname=getfield(m,varnames{inamechoice});
    end
    if(inamechoice==0 || ~ischar(dname))
        dname=fname;
    end
    plotimage3D(seis,t,xline,iline,dname)
end

end

function [itchoice,ixchoice,iychoice]=ambigdialog(ambig,it,ix,iin,varnames,itchoice,ixchoice,iychoice)
hsanfig=gcf;
pos=get(hsanfig,'position');
if(pos(3)<450);pos(3)=450;end
%hdial=dialog;
hdial=figure('windowstyle','modal');
set(hdial,'position',pos,'menubar','none','toolbar','none','numbertitle','off',...
    'name','SANEV mat file input ambiguity dialog','nextplot','new');
indt= it==1;
indx= ix==1;
indin= iin==1;
columnformat={};
data={};
columnnames={};
if(ambig(1)==1)
   columnformat=[columnformat {varnames(indt)'}];
   data=[data varnames{itchoice}];
   columnnames=[columnnames {'time coordinate'}];
end
if(ambig(2)==1)
   columnformat=[columnformat {varnames(indx)'}];
%    if(strcmp(varnames{indt(itchoice)},varnames{indx(ixchoice)}))
%     data=[data varnames{indx(2)}];
%     ixchoice=2;
%    else
    data=[data varnames{ixchoice}];
%    end
   columnnames=[columnnames {'xline coordinate'}];
end
if(ambig(3)==1)
   columnformat=[columnformat {varnames(indin)'}];
%    if(strcmp(varnames{indin(iychoice)},varnames{indx(ixchoice)}))
%       iychoice=ixchoice+1;
%    end
   data=[data varnames{iychoice}];
   columnnames=[columnnames {'inline coordinate'}];
end

ynow=.4;ht=.3;
htab=uitable(hdial,'data',data,'columnformat',columnformat,'columnname',columnnames,...
    'rowname','choose:','columneditable',true,'units','normalized','position',[.1 ynow .8 ht],...
    'tag','table','userdata',varnames);
htab.Position(3)=htab.Extent(3);
htab.Position(4)=htab.Extent(4);
ynow=.6;
msg={'Unable to determine coordinate vectors based on size only.'...
    'Please choose a unique name for each coordinate'};
uicontrol(hdial,'style','text','units','normalized','position',[.1 ynow .8 ht],'string',msg,...
    'tag','msg');
uicontrol(hdial,'style','pushbutton','string','Done','units','normalized',...
    'position',[.1 .1 .3 .1],'callback',@checkambig);

uiwait(hdial)

%function [itchoice,ixchoice,iychoice]=checkambig(~,~)
function checkambig(~,~)
htable=findobj(gcf,'tag','table');
choices=htable.Data;
if(strcmp(htable.Data{1},htable.Data{2})||strcmp(htable.Data{1},htable.Data{3})...
        ||strcmp(htable.Data{2},htable.Data{3}))
    hmsg=findobj(gcf,'tag','msg');
    set(hmsg,'string','You must choose unique names for each!!!','foregroundcolor',[1 0 0],...
        'fontsize',10,'fontweight','bold');
    itchoice=0;
    ixchoice=0;
    iychoice=0;
    return
end
varnames=get(htable,'userdata');
for k=1:length(varnames)
   if(strcmp(varnames{k},choices{1}))
       itchoice=k;
   end
   if(strcmp(varnames{k},choices{2}))
       ixchoice=k;
   end
   if(strcmp(varnames{k},choices{3}))
       iychoice=k;
   end
end
close(gcf)
end

end