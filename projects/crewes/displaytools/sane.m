function sane(action,arg2)
% SANE: Seismic ANalysis Environment (3D)
%
% Just type SANE and the GUI will appear
%
% SANE establishes an analysis environment for 3D seismic datasets. There is currently no facility
% for 2D. SANE is meant to be used in an energy company environment where 3D seismic volumes must be
% compared, QC'd, and perhaps enhanced. To make effective use of SANE you should run this on a
% workstation with lots of RAM. SANE was developed on a computer with 128GB of RAM. It is suggested
% that your workstation should have at least 3 times the RAM of the size of the SEGY volumes that
% you wish to analyze. So, for eaxmple, if you have a 3D dataset that is 10GB as a SEGY file on
% disk, then you should have at least 30GB of RAM available. SANE allows you to read in one or more
% 3D volumes into a project. Its a good idea if the volumes are all somehow related or similar. For
% example maybe they are different processing of the same data. You load these into SANE using the
% "Read SEGY" option and then save the project. Reading SEGY is quite slow but once you save the
% project (as a Matlab binary) further reads are much faster. SANE saves your data internally in
% single precision because that reduces memory and that is how SEGY files are anyway. If your 3D
% dataset has a very irregular patch size, then forming the data into a 3D volume will require
% padding it with lots of zero traces and this can significantly increase memory usage. SANE allows
% you to control which datasets in the project are in memory and which are displayed at any one
% time. Thus it is quite possible to have many more datasets in a project than you can possibly
% display at any one time. Data display is done by sending the data to plotimage3D and you might
% want to check the help for that function. So each dataset you display is in a separate plotimage3D
% window and the windows can be "grouped" to cause them all to show the same view. Plotimage3D can
% show 2D slices of either inline, xline (crossline), or timeslice. In each view there are a number
% of analysis tools available and these are accessed by a right-click of the mouse directly on the
% image plot. Such a right-click brings up a "context menu" of available analysis tools. These tools
% just operate directly on the 2D image that is being displayed. SANE also offers a gradually
% expanding list of "tasks" that are accesible from the "Compute" menu that operate on an entire 3D
% volume and usually produces a same-size 3D volume. Examples are filtering and deconvolution.
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

% To Do:
%
% 
% Plotimage3D needs to be able to send signals to sane. For example when grouping changes
%   - Done: SANE sets the 'tag' of the PI3D figure to 'fromsane' and the userdata to {idata hsane}
%   where idata is the dataset number and hsane is the handle of the SANE window.
%   - PI3D then calls SANE for: group and ungroup: sane('pi3d:group',{idata hsane}). This causes
%   SANE to examine the globals and reset the group buttons as needed
%
% Need ability to edit time shifts and to apply coordinate scalars
%
% Write segy and mat file options. They need a way of choosing the dataset to write.
%
% It would be nice to be able to write out a subset.
%
% Need a way to edit the project structure fields like dx, dy, tshift, filenames, 
%
% Need a way to edit the SEGY text header
% 
% 

% How to edit a SANE Project dataset. SANE project datasets are saved as .mat files with just two
% variables visible at the top level: proj and datasets. proj is a structure with lots of fields and
% datasets is a cell array with the 3D datasets stored in it. The easiest way to edit a dataset from
% the command line is to use the "matfile" facility. If you are unfamiliar with matfile then you
% should read the online help. Suppose saneproject.mat is the name of a SANE project file that needs
% to be edited. Then open it like this
% >> m=matfile('saneproject.mat','writable',true);
% This command does not read the objects in the file, it just opens them. To read the project
% structure into memory do
% >> proj=m.proj
% where I've left off the semicolon to list the structure fields. If you choose to edit this
% structure, beware that there are lots of implied logical connections between the fields and making
% a change in one field can require a corresponding change in another in order that SANE will
% understand. Also, never change the field names. Notice that the field 'filenames' is a cell array
% with a certain length, there is one filename for each 3D survey in the project. Most, but not all,
% of the fields in proj must also have this same length. So, if you make an edit that changes a
% field length, then you must change all of the other fields in the same way.  Also notice that most
% of the fields are cell arrays but some are ordinary arrays. Be sure to preserve this. There is a
% field in proj called 'datasets' that is a cell array to put the seismic datasets in. When you read
% proj from disk like is done here, this will always be empty and the seismic volumes are all in the
% datasets cell array on disk. The field 'isloaded' is an ordinary array of 1 or 0, 1 for each
% dataset. If 1, then the dataset is read from disk into proj.datasets when the project is first
% opened. The isloaded field reflects the load status of things when the project was last saved.
% SANE reads this upon opening the project and puts up a dialog allowing you to decide what is to be
% loaded and what is to be displayed (there is also an isdisplayed field). When a dataset is
% deleted, the various fields are not shortened in length, rather, the corresponding entry in
% datasets is set to null and the fields 'isdeleted' and 'deletedondisk' are set to indicate
% deletion. To read a dataset from the datasets cell array on disk you must know the number of the
% dataset. This is just the integer number of the dataset in the cell array. In determining this, be
% sure to allow for any deleted datasets by checking the isdeleted field. To read in dataset number
% 3 the two-step syntax is
% >> cdat=m.datasets(1,3);
% >> seis=cdat{1};
% The first line here uses regular round brackets even though we are reading from a cell array. This
% is a 'feature' of the matfile behaviour. Note also that you must use two index addressing like
% (1,3) and not simply (3). This line reads the dataset into a cell array of length 1 and the second
% line is required to unpack the dataset into a conventional 3D matrix. This matrix is stored in the
% order expected by plotimage3D which is with time as dimension 1, xline as dimension 2, and inline
% as dimension3. If you wish to write a new dataset into location 3, the syntax is
% >> m.datasets(1,3)={seis};
% Where seis is a 3D volume of the proper size. If you changed the dimensions of this volume in the
% course of altering it, then you must also update the various coordinate fields that have the same
% dimensions. You may also encounter problems if you try to output a SEGY volume from SANE after
% changing the dimensions. This is because SANE remembers the SEGY headers from the original read
% and tries to reuse them. Deletion of a dataset is similar
% >> m.datasets(1,3)={[]};
% which must be followed with
% >> proj.isdeleted(3)=1;
% >> proj.deletedondisk(3)=1;
% >> m.proj=proj;
% Deletion is the one exception to changing the size of a dataset that does not require
% corresponding changes in the coordinate arrays.

% plan for tasks. There will be two types of tasks, those accesible via the context menu of the
% current view in plotimage3D and those accessible by the tasks menu in Sane. I discuss the latter
% here and the former in plotimage. Sane tasks will be applied to entire datasets with the
% subsequent option of saving the result in the project, writing to segy, or discarding. Since many
% tasks in Sane will be similar to those in plotimage3d, there needs to be a mechanism to share
% parmsets. Sane tasks will execute via a callback to an action in sane. Each callback needs to do:
% (1) identify the input dataset, (2) identify the parameters, (3) run the task, and (4) determine
% the disposition of the output.

% To implement a new SANE task, do the following
% 1) Add a new entry to the Compute menu. The tag of this menu will be the name of the task.
% 2) Create a parmset function. This is an internal-to-SANE function that does two things (i) it
% defines the parameters needed to run the task and (ii) it checks a parmset edited by the user for
% validity. See parmsetfilter and parmsetdecon (in this file) for examples.
% 3) Create a new case in the switch in the internal function getparmset (in this file). This new
% case must call the parmset function created in step 2.
% 4) Create a new case in the switch in 'starttask' action for the new task. This function calls the
% internal function sanetask (in this file) that puts up a GUI showing the current parameters in the
% parmset and their current values. The user then changes them and pushes the done button which
% calls the action 'dotask'. Internal function sanetask is meant to automatically adapt to parmsets
% of different lengths and types and (hopefully) will not require modification. 
% 5) In the 'dotask' action there are two switch's than need new cases. The first is the switch
% statment that calls the appropriate parmset function to check the current parmset for validity.
% The second is the switch that actually does the computation. This switch will generally need more
% work and thought than the others.
%
% NOTE: Initially, there are no parmsets saved in a project and so the parmset presented to the
% user is the default one. However, once a task is run, the parmset is saved in the proj structure.
% The next time that same task is run, then the starting parmset is the saved one.
%
% NOTE: At present all of the tasks have two things in common: (i) They are either trace-by-trace
% operation or slice-by-slice. As such it is easy to put them in a loop and save the results in the
% input matrix. This means that if a computation is partially complete, then the input matrix is
% part input and part output. For this reason, if a task is interrupted when partially done, the
% input dataset is unloaded from memory. Rerunning the task will therefore require a reload (which
% is automatic). (ii) They all have the same 4 options for dealing with the output dataset. The
% options are established in the internal function sanetask and are {'Save SEGY','Save SEGY and
% display','Replace input in project','Save in project as new'} . If these two behaviors are not
% appropriate, then more work will be required to implement the new task. The four output options are
% implemented in the action 'dotask' .


% SANE project structure fields
% name ... name of the project
% projfilename ... file name of the project
% projpath ... path of the project.
% filenames ... cell array of filenames for datasets
% paths ... cell array of paths for datasets
% datanames ... cell array of dataset names
% isloaded ... array of flags for loaded or not, one per dataset
% isdisplayed ... array of flags for displayed or not, one per dataset
% xcoord ... cell array of x (xline) coordinate vectors, will be empty if dataset not loaded.
% ycoord ... cell array of y (inline) coordinate vectors, will be empty if dataset not loaded.
% tcoord ... cell array of t (time or depth) coordinate vectors, will be empth if not loaded.
% datasets ... cell array of datasets as 3D matrices, will be empty if not loaded
% xcdp ... cell array of xcdp numbers, empty if not loaded
% ycdp ... cell array of ycdp numbers, empty if not loaded
% dx ... array of physical grid in x direction, will be empty if not loaded
% dy ... array of physical grid in x direction, will be empty if not loaded
% depth ... array of flags, 1 for depth, 0 for time
% texthdr ... cell array of text headers
% texthdrfmt ... cell array of text header formats
% segfmt ... cell array of seg formats 
% byteorder ... cell array of byte orders
% binhdr ... cell array of binary headers
% exthdr ... cell array of extended headers
% tracehdr ... cell array of trace headers
% bindef ... cell array of bindef, nan indicates default behaviour and should be passed to readsegy and writesegy as [];
% trcdef ... cell array of trcdef, nan indicates default behaviour and should be passed to readsegy and writesegy as [];
% segyrev ... cell array of segyrev. nan indicates default behaviour and should be passed to readsegy and writesegy as []; 
% kxline ... cell array of kxlineall values as returned from make3Dvol
% gui ... array of handles of the data panels showing name and status
% rspath ... last used path for reading segy
% wspath ... last used path for writing segy
% rmpath ... last used path for reading matlab
% wmpath ... last used path for writing matlab
% pifigures ... Cell array of plotimage3D figure handles. One for each dataset. Will be empty if not displayed.
% isdeleted ... array of flags signalling deleted or no
% isdeletedondisk ... array of flags signalling deleted on disk or not
% saveneeded ... array of flags indicating a dataset needs to be saved
% parmsets ... cell array of parmsets which are also cell arrays. A parmset holds parameters for functions like
%           filters, decon, etc. Each parmset is a indefinite length cell array of name value triplets.
%           However the first entry is always a string giving the name of the parmset. Thus a
%           parmset always has length 3*nparms+1. A name value triple consists of (1) parameter
%           name, (2) parameter value, (3) tooltip string. The latter being a hint or instruction.
%           The parameter value is either a string or a cell array. If the parameter is actually a
%           number then it is read from the string with str2double. If the parameter is a choice,
%           then it is encoded as a cell array like this: param={'choice1' 'choice2' 'choice3' 1}.
%           The last entry is numeric and in this example meanse that the default is choice1.



%userdata assignments
%hfile ... the project structure
%hmpan ... (the master panel) {hpanels geom thispanel}
%hpan ... the idata (data index) which is the number of the dataset for that data panel
%hreadsegy ... path for the most recently read segy
%hreadmat ... path for the most recently read .mat
%any plotimage figure ... {idata hsane}, idata is the number of the dataset, hsane is the sane figure
%any figure menu in "view" ... the handle of the figure the menu refers to
% 

global PLOTIMAGE3DTHISFIG PLOTIMAGE3DFIGS HWAIT CONTINUE


if(nargin<1)
    action='init';
end

if(strcmp(action,'init'))
%     test=findsanefig;
%     if(~isempty(test))
%         msgbox('You already have SANE running, only one at a time please');
%         return
%     end
    hsane=figure;
    ssize=get(0,'screensize');
    figwidth=1000;
    figheight=floor(ssize(4)*.4);
    xnot=(ssize(3)-figwidth)*.5;
    ynot=(ssize(4)-figheight)*.5;
    set(hsane,'position',[xnot,ynot,figwidth,figheight],'tag','sane');
    set(hsane,'menubar','none','toolbar','none','numbertitle','off','name','SANE New Project',...
        'nextplot','new','closerequestfcn','sane(''close'')');
    
    hfile=uimenu(hsane,'label','File','tag','file');
    uimenu(hfile,'label','Load existing project','callback','sane(''loadproject'');','tag','loadproject');
    uimenu(hfile,'label','Save project','callback','sane(''saveproject'');');
    uimenu(hfile,'label','Save project as ...','callback','sane(''saveprojectas'');');
    uimenu(hfile,'label','New project','callback','sane(''newproject'');');
    hread=uimenu(hfile,'Label','Read datasets','tag','read');
    uimenu(hread,'label','*.sgy file','callback','sane(''readsegy'');','tag','readsegy');
    uimenu(hread,'label','Multiple sgy files','callback','sane(''readmanysegy'');','tag','readmanysegy');
    uimenu(hread,'label','*.mat file','callback','sane(''readmat'');','tag','readmat');
    hwrite=uimenu(hfile,'label','Write datasets');
    uimenu(hwrite,'label','*.sgy file','callback','sane(''writesegy'');');
    uimenu(hwrite,'label','*.mat file','callback','sane(''writemat'');');
    uimenu(hfile,'label','Quit','callback','sane(''close'')');
    
    
    uimenu(hsane,'label','View','tag','view');
    
    hcompute=uimenu(hsane,'label','Compute','tag','compute');
    uimenu(hcompute,'label','Bandpass filter','callback','sane(''starttask'')','tag','filter');
    uimenu(hcompute,'label','Spiking decon','callback','sane(''starttask'')','tag','spikingdecon');
    uimenu(hcompute,'label','Wavenumber lowpass filtering','callback','sane(''starttask'')','tag','wavenumber','enable','on');
    if(exist('tvfdom','file')==2)
        uimenu(hcompute,'label','Dominant frequency volumes','callback','sane(''starttask'')','tag','fdom','enable','on');
    end
    uimenu(hcompute,'label','Phase maps','callback','sane(''starttask'')','tag','phasemap','enable','off');
    
    proj=makeprojectstructure;
    set(hfile,'userdata',proj);
    
    x0=.05;width=.2;height=.05;
    ysep=.02;
    xsep=.02;
    xnow=x0;
    ynow=1-height-ysep;
    fs=10;
    %a message area
    uicontrol(gcf,'style','text','string','Load an existing project or read a SEGY file','units','normalized',...
        'position',[xnow,ynow,.8,height],'tag','message','fontsize',fs,'fontweight','bold');
    %project name display
    ynow=ynow-height-ysep;
    uicontrol(gcf,'style','text','string','Project Name:','tag','name_label','units','normalized',...
        'position',[xnow,ynow-.25*height,width,height],'fontsize',fs,'horizontalalignment','right');
    xnow=xnow+xsep+width;
    uicontrol(gcf,'style','edit','string',proj.name,'tag','project_name','units','normalized',...
        'position',[xnow,ynow,2*width,height],'fontsize',fs,'callback','sane(''projectnamechange'');');
    panelwidth=1-2*x0;
    panelheight=1.2*height;
    xnow=x0;
    ynow=ynow-height-ysep;
    %the master panel
    hmpan=uipanel(gcf,'tag','master_panel','units','normalized','position',...
        [xnow ynow panelwidth panelheight]);
    xn=0;yn=0;wid=.5;ht=.8;ht2=1.1;
    ng=.94*ones(1,3);
    dg=.7*ones(1,3);
    uicontrol(hmpan,'style','text','string','Dataset','tag','dataset_label','units','normalized',...
        'position',[xn,yn,wid,ht],'horizontalalignment','center','fontsize',fs,'backgroundcolor',ng);
    uicontrol(hmpan,'style','text','string','','units','normalized','position',...
        [xn+wid+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
    xn=xn+wid+xsep;
    wid2=(1-wid-3*xsep)/4.5;
    uicontrol(hmpan,'style','text','string','Info','tag','info_label','units','normalized',...
        'position',[xn,yn,.75*wid2,ht],'horizontalalignment','center','fontsize',fs,'backgroundcolor',ng);
    uicontrol(hmpan,'style','text','string','','units','normalized','position',...
        [xn+.75*wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
    xn=xn+.75*wid2+xsep;
    uicontrol(hmpan,'style','text','string','In memory','tag','memory_label','units','normalized',...
        'position',[xn,yn,wid2,ht],'horizontalalignment','center','fontsize',fs,'backgroundcolor',ng);
    uicontrol(hmpan,'style','text','string','','units','normalized',...
        'position',[xn+wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
    xn=xn+wid2+xsep;
    uicontrol(hmpan,'style','text','string','Displayed','tag','display_label','units','normalized',...
        'position',[xn,yn,wid2,ht],'horizontalalignment','center','fontsize',fs,'backgroundcolor',ng);
    uicontrol(hmpan,'style','text','string','','units','normalized',...
        'position',[xn+wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
    xn=xn+wid2+xsep;
    uicontrol(hmpan,'style','text','string','Delete','tag','delete_label','units','normalized',...
        'position',[xn,yn,.5*wid2,ht],'horizontalalignment','center','fontsize',fs,'backgroundcolor',ng);
    uicontrol(hmpan,'style','text','string','','units','normalized',...
        'position',[xn+.5*wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
    xn=xn+.5*wid2+xsep;
    uicontrol(hmpan,'style','text','string','Group','tag','group_label','units','normalized',...
        'position',[xn,yn,.5*wid2,ht],'horizontalalignment','center','fontsize',fs,'backgroundcolor',ng);
    %userdata of hmpan will be a cell array. The first entry is an array of panel handles, one for each dataset in the project
    %the second entry is geometry information: panelwidth panelheight wid ht xsep ysep
    ysep=.01;
    set(hmpan,'userdata',{[],[panelwidth panelheight xnow ynow wid ht xsep ysep ynow]});
        
    
elseif(strcmp(action,'readsegy'))
    %read in a segy dataset
    hsane=findsanefig;
    hfile=findobj(hsane,'label','File');
    hreadsegy=findobj(hsane,'tag','readsegy');
    hmsg=findobj(hsane,'tag','message');
    startpath=get(hreadsegy,'userdata');
    if(isempty(startpath))
        spath='*.sgy';
    else
        spath=[startpath '*.sgy'];
    end
    [fname,path]=uigetfile(spath,'Choose the .sgy file to import');
    if(fname==0)
        return
    end
    ind=strfind(fname,'.sgy');
    if(isempty(ind))
        ind=strfind(fname,'.segy');
        if(isempty(ind))
            msgbox('Chosen file is not a either .sgy or .segy, cannot proceed');
            return;
        end
    end
    dname=fname(1:ind(1)-1);
    waitsignalon
    segyrev=1;
    set(hmsg,'string',['Reading SEGY dataset ' path fname]);
    [seis,segyrev,dt,segfmt,texthdrfmt,byteorder,texthdr,binhdr,...
        exthdr,tracehdr,bindef,trcdef] =readsegy([path fname],[],segyrev,[],[],...
        [],[],[],[],[],hsane);

    t=double(dt)*(0:size(seis,1)-1)';
    %determine if time or depth
    dt=abs(t(2)-t(1));
    depthflag=0;
    if(dt>.02)
        depthflag=1;
    end
    if(isfield(tracehdr,'InlineNum')) %this is a field for segyrev1 but not rev0
        %we need InlineNum and XlineNum to be populated
        y=double(tracehdr.InlineNum);
        x=double(tracehdr.XlineNum);
        if(sum(abs(x))==0)
            %try the Kingdom header locations
            if(isfield(tracehdr,'TrcNumCDP'))
                x=double(tracehdr.TrcNumCDP);
                y=double(tracehdr.SourcePoint);
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
        if(isfield(tracehdr,'TrcNumCDP'))
            x=double(tracehdr.TrcNumCDP);
            y=double(tracehdr.SourcePoint);
        else
            x=0;
        end
        if(sum(abs(x)==0))
            msgbox('Xline/Inline header information is empty. Cannot proceed');
            waitsignaloff;
            return
        end
    end
    if(isfield(tracehdr,'CdpX'))
        cdpx=tracehdr.CdpX;
        cdpy=tracehdr.CdpY;
    end
    [seis3D,xline,iline,xcdp,ycdp,kxline]=make3Dvol(seis,x,y,cdpx,cdpy);
    
    dx=max([abs(xcdp(2)-xcdp(1)), 1]);
    dy=max([abs(ycdp(2)-ycdp(1)), 1]);
    
    %ask for a few things
    tshift=0;
    if(depthflag==1)
        q4='Datum shift (depth units)';
    else
        q4='Datum shift (seconds)';
    end
    q={'Specify dataset name:','Physical distance between crosslines:','Physical distance between inlines:',q4};
    a={dname num2str(dy) num2str(dx) num2str(tshift)};
    set(hmsg,'string','SANE is waiting for your response. Look for the dialog box');
    a=askthingsle('name','Please double check these values','questions',q,'answers',a);
    if(isempty(a))
        msgbox('SEGY input aborted');
        set(hmsg,'string','SEGY input was cancelled');
        waitsignaloff;
        return;
    end
    dname=a{1};
    dy=str2double(a{2});
    dx=str2double(a{3});
    tshift=str2double(a{4});
    if(isnan(tshift));tshift=0;end
    %insist on a positive number for dx and dy
    if(isnan(dy)); dy=0; end
    if(isnan(dx)); dx=0; end
    if(dx<0); dx=0; end
    if(dy<0); dy=0; end
    while(dx*dy==0)
        q={'Specify dataset name:','Physical distance between crosslines:','Physical distance between inlines:',q4};
        a={dname num2str(dy) num2str(dx) num2str(tshift)};
        a=askthingsle('name','Inline and crossline distances must be positive numbers!!','questions',q,'answers',a);
        if(isempty(a))
            msgbox('SEGY input aborted');
            set(hmsg,'string','SEGY input was cancelled');
            waitsignaloff;
            return;
        end
        dname=a{1};
        dy=str2double(a{2});
        dx=str2double(a{3});
        if(isnan(dy)); dy=0; end
        if(isnan(dx)); dx=0; end
        if(dx<0); dx=0; end
        if(dy<0); dy=0; end
    end
    if(depthflag==0 && tshift>10)
        tshift=tshift/1000;%assume they mean milliseconds
    end
    
    hpan=newdatapanel(dname,1,1);
    
    %update the project structure
    
    proj=get(hfile,'userdata');
    nfiles=length(proj.filenames)+1;
    proj.filenames{nfiles}=fname;
    proj.paths{nfiles}=path;
    proj.datanames{nfiles}=dname;
    proj.isloaded(nfiles)=1;
    proj.isdisplayed(nfiles)=1;
    proj.xcoord{nfiles}=xline;
    proj.ycoord{nfiles}=iline;
    proj.tcoord{nfiles}=t+tshift;
    proj.tshift(nfiles)=tshift;
    proj.xcdp{nfiles}=xcdp;
    proj.ycdp{nfiles}=ycdp;
    proj.dx(nfiles)=dx;
    proj.dy(nfiles)=dy;
    proj.depth(nfiles)=depthflag;
    proj.texthdr{nfiles}=texthdr;
    proj.texthdrfmt{nfiles}=texthdrfmt;
    proj.segfmt{nfiles}=segfmt;
    proj.byteorder{nfiles}=byteorder;
    proj.tracehdr{nfiles}=tracehdr;
    proj.binhdr{nfiles}=binhdr;
    proj.exthdr{nfiles}=exthdr;
    proj.bindef{nfiles}=bindef;
    proj.trcdef{nfiles}=trcdef;
    proj.segyrev{nfiles}=segyrev;
    proj.kxline{nfiles}=kxline;
    proj.datasets{nfiles}=seis3D;
    proj.gui{nfiles}=hpan;
    proj.isdeleted(nfiles)=0;
    proj.deletedondisk(nfiles)=0;
    proj.saveneeded(nfiles)=1;
    

    %call plotimage3D
    plotimage3D(seis3D,t,xline,iline,dname,'seisclrs',dx,dy);
    set(gcf,'tag','fromsane','userdata',{nfiles hsane});
    set(gcf,'closeRequestFcn','sane(''closepifig'');')
    hview=findobj(hsane,'tag','view');
    uimenu(hview,'label',dname,'callback','sane(''popupfig'');','userdata',gcf);
    proj.pifigures{nfiles}=gcf;
    %save the path
    set(hreadsegy,'userdata',path);
    
    %save the project structure
    set(hfile,'userdata',proj);
    figure(hsane)
    waitsignaloff
    set(hmsg,'string',['File ' fname ' imported and displayed. Data will be written to disk when you save the project.'])
    
elseif(strcmp(action,'readmanysegy'))
    %Here we read in a bunch of SEGY's at once. 
    %If the existing project is not new, then we first ask if we want to add to it. Otherwise we suggest save.  
    hsane=findsanefig;
    hmsg=findobj(hsane,'tag','message');
    hfile=findobj(hsane,'tag','file');
    proj=hfile.UserData;
    newproject=true;
    if(~isempty(proj.filenames))
       Q='Merge new data into existing project?';
       Q0='What about the current project?';
       A1='Yes';
       A2='No, start new project';
       A3='Cancel, I forgot to save';
       answer=questdlg(Q,Q0,A1,A2,A3,A1);
       if(strcmp(answer,A3))
           set(hmsg,'string','Multilple SEGY load cancelled, Save your project first');
           return;
       elseif(strcmp(answer,A2))
           sane('newproject');
       else
           newproject=false;
       end
    end
    %put up dialog
    set(hmsg,'string','Beginning Multiple SEGY load');
    multiplesegyload(newproject);
    return;
elseif(strcmp(action,'readmanysegy2'))
    hdial=gcf;
    hsane=findsanefig;
    hmsg=findobj(hsane,'tag','message');
    set(hmsg,'string','Beginning multiple SEGY read');
    hfile=findobj(hsane,'tag','file');
    proj=hfile.UserData;
    ht=findobj(hdial,'tag','table');
    data=ht.Data;
    datapaths=ht.UserData;
    if(isempty(data))
        msgbox('You need to choose some datasets to import','Oooops!');
        return;
    end
    hprojfile=findobj(hdial,'tag','projsavefile');
    projfile=get(hprojfile,'string');
    
    if(strcmpi('undefined',projfile))
        msgbox('You need define the project save file','Oooops!');
        return;
    end
    proj.projfilename=projfile;
    proj.projpath=get(hprojfile,'userdata');
    matobj=matfile([proj.projpath proj.projfilename],'writable',true);%open the mat file
    nd=length(proj.datanames);%number of datasets currently in project
    ndnew=size(data,1);%number of new datasets
    proj=expandprojectstructure(proj,ndnew);
    set(hfile,'userdata',proj);
    delete(hdial);
    t0=clock;
    figure(hsane);
    waitsignalon
    for k=nd+1:nd+ndnew
        filename=data{k-nd,1};
        datapath=datapaths{k-nd};
        proj.datanames{k}=data{k-nd,2};
        dispopt=data{k-nd,3};
        if(dispopt)
            proj.isdisplayed(k)=1;
            proj.isloaded(k)=1;
        end
        tshift=str2double(data{k-nd,4});
        if(isnan(tshift))
            tshift=0;
        elseif(tshift>10)
            tshift=tshift/1000;
        end
        proj.tshift(k)=tshift;
        
        set(hmsg,'string',['reading file ' int2str(k-nd) ' of ' int2str(ndnew) ', ' filename ...
            ' from path ' datapath ', see main Matlab window for progress']);
        drawnow
        [seis,sgrv,dt,proj.segfmt{k},proj.texthdrfmt{k},proj.byteorder{k},proj.texthdr{k},proj.binhdr{k},...
            proj.exthdr{k},proj.tracehdr{k},proj.bindef{k},proj.trcdef{k}] =readsegy([datapath filename],[],1,[],[],...
            [],[],[],[],[],1); %#ok<ASGLU>
        
        t=dt*(0:size(seis,1)-1)'+tshift;
        proj.tcoord{k}=t;
        dt=abs(t(2)-t(1));
        if(dt>.02)
            proj.depth(k)=1;
        end
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
        [seis3D,proj.xcoord{k},proj.ycoord{k},xcdp,ycdp,proj.kxline{k}]=make3Dvol(seis,x,y,cdpx,cdpy);
        
        proj.dx(k)=max([abs(xcdp(2)-xcdp(1)), 1]);
        proj.dy(k)=max([abs(ycdp(2)-ycdp(1)), 1]);
        proj.xcdp{k}=xcdp;
        proj.ycdp{k}=ycdp;
        proj.gui{k}=newdatapanel(proj.datanames{k},proj.isloaded(k),proj.isdisplayed(k));%make a data panel
        matobj.datasets(1,k)={seis3D};
        
        if(proj.isloaded(k)==1)
            proj.datasets{k}=seis3D;
            plotimage3D(seis3D,t,proj.xcoord{k},proj.ycoord{k},proj.datanames{k},'seisclrs',proj.dx(k),proj.dy(k));
            set(gcf,'tag','fromsane','userdata',{k hsane});
            set(gcf,'closeRequestFcn','sane(''closepifig'');')
            hview=findobj(hsane,'tag','view');
            uimenu(hview,'label',proj.datanames{k},'callback','sane(''popupfig'');','userdata',gcf);
            proj.pifigures{k}=gcf;
            figure(hsane)
        end
        
        
        
        tnow=clock;
        timeused=etime(tnow,t0);
        
        set(hmsg,'string',['Finished reading file ' int2str(k-nd) ' of ' int2str(ndnew) ', '...
            filename ', Time used = ' num2str(timeused/60) ' minutes']);
    end
    set(hfile,'userdata',proj);
    for k=nd+1:nd+ndnew
        proj.datasets{k}=[];%don't save a dataset in the project structure
        proj.pifigures{k}={};%don't want to save a graphics handle
        proj.gui{k}={};
    end
    matobj.proj=proj;%save project structure
    
    
    waitsignaloff
    
elseif(strcmp(action,'selectnewdataset'))
%     hsane=findsanefig;
    hdial=gcf;
%     hfile=findobj(hsane,'label','File');
%     hreadsegy=findobj(hsane,'tag','readsegy');
%     hmsg=findobj(hsane,'tag','message');
%     hreadmany=findobj(hsane,'tag','readmanysegy');
    hnew=findobj(hdial,'tag','new');
    startpath=get(hnew,'userdata');
    if(isempty(startpath))
        spath='*.sgy';
    else
        spath=[startpath '*.sgy'];
    end
    [fname,path]=uigetfile(spath,'Choose the .sgy file to import');
    if(fname==0)
        return
    end
    ind=strfind(fname,'.sgy');
    if(isempty(ind))
        ind=strfind(fname,'.segy');
        if(isempty(ind))
            msgbox('Chosen file is not a either .sgy or .segy, cannot proceed');
            return;
        end
    end
    dname=fname(1:ind(1)-1);
    ht=findobj(hdial,'tag','table');
    data=ht.Data;
    nd=size(data,1);
    data=[data;{fname,dname,false,0.0}];
    ht.Data=data;
    paths=ht.UserData;
    if(~iscell(paths))
        paths={};
    end
    paths{nd+1}=path;
    ht.UserData=paths;%cell array of all paths
    set(hnew,'userdata',path);%the is always the path of the last chosen dataset
elseif(strcmp(action,'defineprojectsavefile'))
%     hsane=findsanefig;
    hdial=gcf;
    %hfile=findobj(hsane,'label','File');
    %hreadsegy=findobj(hsane,'tag','readsegy');
%     hmsg=findobj(hsane,'tag','message');
    hdialmsg=findobj(hdial,'tag','dialmsg');
%     hreadmany=findobj(hsane,'tag','readmanysegy');
    hnew=findobj(hdial,'tag','new');
    ht=findobj(hdial,'tag','table');
    paths=get(ht,'userdata');
    if(isempty(paths))
        startpath=get(hnew,'userdata');
    else
        startpath=paths{1};
    end
    if(isempty(startpath))
        spath='*.mat';
    else
        spath=[startpath '*.mat'];
    end
    [fname,path]=uiputfile(spath,'Specify the Project .mat file');
    if(fname==0)
        set(hdialmsg,'string','UNDEFINED');
        return
    end
    ind=strfind(fname,'.sgy');
    if(~isempty(ind))
        if(isempty(ind))
            msgbox('Chosen file mut be a .mat file not a .sgy. Try again');
            set(hdialmsg,'String','Project save file must be a .mat file');
            return;
        end
    end
%     ind=strfind(fname,'.mat');
%     if(isempty(ind)) %#ok<STREMP>
    if(~contains(fname,'.mat'))
        fname=[fname '.mat'];
    end
    if(exist([path fname],'file'))
        response=questdlg('The specified prohect file already exists. Overwrite?','Project file question.',...
            'Yes','No','Yes');
        if(strcmp(response','No'))
           set(hdialmsg,'string','Choose a different Project file');
           return;
        end
    end
    hprojfile=findobj(hdial,'tag','projsavefile');
    set(hprojfile,'string',fname,'userdata',path);
elseif(strcmp(action,'cancelmultipleload'))
    delete(gcf)
elseif(strcmp(action,'reloaddataset'))
    hsane=findsanefig;
    hmsg=findobj(hsane,'tag','message');
    %this is called from a datapanel to load a dataset not in memory. The data panel is identified
    %by the third entry of the userdata of the master panel
    hmpan=findobj(hsane,'tag','master_panel');
    udat=hmpan.UserData;
    idata=udat{3};%this is the dataset we are loading
    %hpan=udat{1}(udat{3});
    hfile=findobj(hsane,'tag','file');
    proj=hfile.UserData;
    
    set(hmsg,'string',['Recalling dataset ' proj.datanames{idata} ' from disk'])
    waitsignalon
    matobj=matfile([proj.projpath proj.projfilename]);
    cseis3D=matobj.datasets(1,idata);%this reads from disk
    seis3D=cseis3D{1};
    t=proj.tcoord{idata};
    xline=proj.xcoord{idata};
    iline=proj.ycoord{idata};
    dname=proj.datanames{idata};
    
    %update the project structure
    proj.datasets{idata}=seis3D;
    proj.isloaded(idata)=1;
    proj.saveneeded(idata)=0;
    %call plotimage3D
    if(proj.isdisplayed(idata)==0)%don't display it a second time
        plotimage3D(seis3D,t,xline,iline,dname,'seisclrs',proj.dx(idata),proj.dy(idata));
        set(gcf,'tag','fromsane','userdata',{idata hsane});
        set(gcf,'closeRequestFcn','sane(''closepifig'');')
        hview=findobj(hsane,'tag','view');
        uimenu(hview,'label',dname,'callback','sane(''popupfig'');','userdata',gcf);
        proj.pifigures{idata}=gcf;
        proj.isdisplayed(idata)=1;
    end
    memorybuttonon(idata);
    waitsignaloff;
    
    set(hfile,'userdata',proj);
    figure(hsane);
    
    set(hmsg,'string',['Dataset ' dname ' reloaded']);
    
elseif(strcmp(action,'readmat'))
    hsane=findsanefig;
    hreadmat=gcbo;
    %basic idea is that the .mat file can contain one and only one 3D matrix. It should also contain
    %coordinate vectors for each of the 3 dimensions. So, an abort will occur if there are more than
    %1 3D matrices (or none), and also if there are no possible coordinate vectors for the 3
    %dimensions. If there are coordinate vectors but it is not clear which is which, then an
    %ambiguity dialog is put up to resolve this.
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
            if(strfind(varnames{k},'dname'))
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
        msgbox(failmsg,'Sorry, dataset is not compaible with SANE.')
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
    %ok now get stuff from the matfile
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
    
    %determine if time or depth
    dt=abs(t(2)-t(1));
    depthflag=0;
    if(dt>.02)
        depthflag=1;
    end
    
    %ask for a few things
    tshift=0;
    if(depthflag==1)
        q4='Datum shift (depth units)';
    else
        q4='Datum shift (seconds)';
    end
    dx=1;
    dy=1;
    q={'Specify dataset name:','Physical distance between crosslines:','Physical distance between inlines:',q4};
    a={dname num2str(dy) num2str(dx) num2str(tshift)};
    a=askthingsle('name','Please double check these values','questions',q,'answers',a);
    if(isempty(a))
        msgbox('SEGY input aborted');
        return;
    end
    dname=a{1};
    dy=str2double(a{2});
    dx=str2double(a{3});
    tshift=str2double(a{4});
    if(isnan(tshift));tshift=0;end
    %insist on a positive number for dx and dy
    if(isnan(dy)); dy=0; end
    if(isnan(dx)); dx=0; end
    if(dx<0); dx=0; end
    if(dy<0); dy=0; end
    while(dx*dy==0)
        q={'Specify dataset name:','Physical distance between crosslines:','Physical distance between inlines:',q4};
        a={dname num2str(dy) num2str(dx) num2str(tshift)};
        a=askthingsle('name','Inline and crossline distances must be positive numbers!!','questions',q,'answers',a);
        if(isempty(a))
            msgbox('.mat input aborted');
            return;
        end
        dname=a{1};
        dy=str2double(a{2});
        dx=str2double(a{3});
        if(isnan(dy)); dy=0; end
        if(isnan(dx)); dx=0; end
        if(dx<0); dx=0; end
        if(dy<0); dy=0; end
    end
    if(depthflag==0 && tshift>2)
        tshift=tshift/1000;%assume they mean milliseconds
    end
    
    hpan=newdatapanel(dname,1,1);
    
    %update the project structure
    hfile=findobj(hsane,'label','File');
    proj=get(hfile,'userdata');
    nfiles=length(proj.filenames)+1;
    proj.filenames{nfiles}=fname;
    proj.paths{nfiles}=path;
    proj.datanames{nfiles}=dname;
    proj.isloaded(nfiles)=1;
    proj.isdisplayed(nfiles)=1;
    proj.xcoord{nfiles}=xline;
    proj.ycoord{nfiles}=iline;
    proj.tcoord{nfiles}=t+tshift;
    proj.tshift(nfiles)=tshift;
    proj.xcdp{nfiles}=dx*(1:length(xline));
    proj.ycdp{nfiles}=dy*(1:length(iline));
    proj.dx(nfiles)=dx;
    proj.dy(nfiles)=dy;
    proj.depth(nfiles)=depthflag;
    proj.texthdr{nfiles}='Data from .mat file. No text header available';
    proj.texthdrfmt{nfiles}=nan;
    proj.segfmt{nfiles}=nan;
    proj.byteorder{nfiles}=nan;
    proj.tracehdr{nfiles}=nan;
    proj.binhdr{nfiles}=nan;
    proj.exthdr{nfiles}=nan;
    proj.kxline{nfiles}=nan;
    proj.datasets{nfiles}=seis;
    proj.gui{nfiles}=hpan;
    proj.isdeleted(nfiles)=0;
    proj.deletedondisk(nfiles)=0;
    proj.saveneeded(nfiles)=1;
    
    %save the path
    set(hreadmat,'userdata',path);
    
    plotimage3D(seis,t,xline,iline,dname,'seisclrs',dx,dy);
    set(gcf,'tag','fromsane','userdata',{nfiles hsane});
    set(gcf,'closeRequestFcn','sane(''closepifig'');')
    hview=findobj(hsane,'tag','view');
    uimenu(hview,'label',dname,'callback','sane(''popupfig'');','userdata',gcf);
    proj.pifigures{nfiles}=gcf;
    set(hfile,'userdata',proj);
    figure(hsane)
    
elseif(strcmp(action,'datainfo'))
    hinfo=gco;
    hsane=findsanefig;
    hpan=get(hinfo,'parent');
    idata=get(hpan,'userdata');
    hfile=findobj(hsane,'label','File');
    proj=get(hfile,'userdata');
    %msgbox(hproj.texthdr{dataindex});
    pos=get(hsane,'position');
    x0=pos(1);
    y0=max([1 pos(2)-1.2*pos(4)]);
    hinfofig=figure('position',[x0 y0 .7*pos(3) 1.4*pos(4)],'name',...
        ['Info for dataset ' proj.datanames{idata}],'menubar','none','toolbar','none','numbertitle','off');
    uicontrol(hinfofig,'style','text','units','normalized','position',[.1 .1 .8 .7],...
        'string',proj.texthdr{idata},'HorizontalAlignment','left');
    nt=length(proj.tcoord{idata});
    nx=length(proj.xcoord{idata});
    ny=length(proj.ycoord{idata});
    dt=abs(proj.tcoord{idata}(2)-proj.tcoord{idata}(1));
    dx=proj.dx(idata);
    dy=proj.dy(idata);
    tshift=proj.tshift(idata);
    ntraces=nx*ny;
    tmax=proj.tcoord{idata}(end)-proj.tcoord{idata}(1);
    mb=round(ntraces*nt*4/10^6);
    datasummary=cell(1,4);
    datasummary{1,1}=['Dataset consists of ' int2str(ntraces) ' traces, ' num2str(tmax) ' seconds long, datum shift= ' num2str(tshift)];
    datasummary{1,2}=['Number of inlines= ' int2str(ny) ', number of crosslines= ' int2str(nx) ', number of time samples=' int2str(nt)];
    datasummary{1,3}=['Time sample size= ' num2str(dt) ' seconds, inline separation= ' num2str(dy) ', crossline separation= ' num2str(dx)];
    datasummary{1,4}=['Dataset size (without headers) = ' int2str(mb) ' megabytes'];
    datasummary{1,5}='SEGY text header follows:';
    uicontrol(hinfofig,'style','text','units','normalized','position',[.1 .85 .8 .1],...
        'string',datasummary,'HorizontalAlignment','left','fontsize',10);
elseif(strcmp(action,'projectnamechange'))
    hfile=findobj(gcf,'tag','file');
    proj=get(hfile,'userdata');
    hprojname=gcbo;
    proj.name=get(hprojname,'string');
    set(gcf,'name',['SANE, Project: ' proj.name])
    set(hfile,'userdata',proj);
elseif(strcmp(action,'saveproject')||strcmp(action,'saveproj'))
    % The project file is a mat file with two variables: proj and datasets. Proj is a structure and
    % datasets is a cell array. Proj has a field called datasets but this is always saved to disk as
    % null and the datasets are sames separately in the cell array. Also in the proj structure are
    % fields isloaded and isdisplayed. When a project is loaded, only proj is read at first and a
    % dialog is presented showing which datasets were previously loaded and/or displayed. The user
    % can then choose to load and display as before or make any changes desired. Datasets that are
    % currently in memory are included in project as proj.datasets. All datasets, in memory or not,
    % are found in the datasets cell array. When a dataset is moved out of memory, it is saved into
    % the datasets array where it can be retrieved when desired. Syntax for retrieving a dataset
    % matobj=matfile([path filename]);
    % cdataset=matobj.datasets(1,thisdataset);%variable thisdataset is the index of the dataset
    % dataset=cdataset{1};
    % Syntax for saveing a dataset
    % matobj=matfile([path filename],'writable',true);
    % matobj.datasets(1,thisdataset)={dataset};
    % There does not appear to be a need to formally close a mat file after writing to it.
    % There does not seem to be a way to load a portion of a dataset. I have to load it all at once.
    %
    
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    hmsg=findobj(hsane,'tag','message');
    proj=get(hfile,'userdata');
    filename=proj.projfilename;
    path=proj.projpath;
    %newprojectfile=0;
    if(isempty(filename) || isempty(path))
        [filename, path] = uiputfile('*.mat', 'Specify project filename');
        if(filename==0)
            msgbox('Project save cancelled');
            return;
        end
        %ind=strfind(filename,'.mat');
        %if(isempty(ind))
        if(~contains(filename,'.mat'))
            filename=[filename '.mat'];
        end
        proj.projfilename=filename;
        proj.projpath=path;
        set(hfile,'userdata',proj);
        %newprojectfile=1;
        if(exist([path filename],'file'))
            delete([path filename]);
        end
    end
    set(hmsg,'string','Saving project ...')
    waitsignalon
    %plan: strip out the datasets from the project. Make sure we have all of the datasets saved in
    %the datasets cell array on disk. Active datasets will be automatically moved into the proj
    %structure on loading.
    datasets=proj.datasets;
    pifigs=proj.pifigures;
    gui=proj.gui;
    proj.datasets={};%save it as an empty cell
    proj.pifigures={};%don't save graphics handles
    proj.gui={};
    disp('opening mat file')
    matobj=matfile([path filename],'Writable',true);
    
   
    isave=proj.saveneeded;
    %save any new datasets
    if(sum(isave)>0)
        disp('saving datasets')
        set(hmsg,'string','Saving new datasets');
        for k=1:length(isave)
            if(isave(k)==1)
                if(isempty(datasets(1,k)))
                    error('attempt to save empty dataset');
                end
                matobj.datasets(1,k)=datasets(1,k);%writes to disk
            end
        end
    end
    %check for newly deleted datasets
    ndatasets=length(datasets);
    for k=1:ndatasets
        if(proj.isdeleted(k)==1 && proj.deletedondisk(k)==0)
            matobj.datasets(1,k)={[]};
            proj.deletedondisk(k)=1;
        end
    end
    disp('writing project structure')
    set(hmsg,'string','Writing project structure');
    proj.saveneeded=zeros(1,ndatasets);
    matobj.proj=proj;%this writes the project structure
    proj.datasets=datasets;
    proj.pifigures=pifigs;
    proj.gui=gui;
    set(hfile,'userdata',proj)
    waitsignaloff
    set(hsane,'name',['SANE, Project file: ' filename]);
    set(hmsg,'string',['Project ' filename ' saved'])
    
elseif(strcmp(action,'saveprojectas'))
    % In this case we squeeze out any deleted datasets
    
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    hmsg=findobj(hsane,'tag','message');
%     filename=proj.projfilename;
%     path=proj.projpath;
    [filename, path] = uiputfile('*.mat', 'Specify project filename');
    if(filename==0)
        msgbox('Project save cancelled');
        return;
    end
    if(exist([path filename],'file'))
        delete([path filename]);
    end
    %     ind=strfind(filename,'.mat');
    %     if(isempty(ind))
    if(~contains(filename,'.mat'))
        filename=[filename '.mat'];
    end
    if(exist([path filename],'file'))
        delete([path filename]);
    end
    set(hmsg,'string',['Saving project into ' path filename])
    %open the old file
    mOld=matfile([proj.projpath proj.projfilename]);
    
%     proj.projfilename=filename;
%     proj.projpath=path;
    set(hfile,'userdata',proj);
    waitsignalon
    nfiles=length(proj.datanames);
    for k=1:nfiles
       if(proj.isdeleted(k)~=1)
           if(isempty(proj.datasets{k}))
               proj.datasets(1,k)=mOld.datasets(1,k);
           end
       end
    end
    %ok, proj,datasets is fully loaded with those that we are keeping
    ind=proj.isdeleted~=1;
    datasets=proj.datasets(ind);%these are the keepers
    newproj=makeprojectstructure;
    newproj.projfilename=filename;
    newproj.projpath=path;
    newproj.filenames=proj.filenames(ind);
    newproj.paths=proj.paths(ind);
    newproj.datanames=proj.datanames(ind);
    newproj.isloaded=proj.isloaded(ind);
    newproj.isdisplayed=proj.isdisplayed(ind);
    newproj.xcoord=proj.xcoord(ind);
    newproj.ycoord=proj.ycoord(ind);
    newproj.tcoord=proj.tcoord(ind);
    newproj.tshift=proj.tshift(ind);
    newproj.xcdp=proj.xcdp(ind);
    newproj.ycdp=proj.ycdp(ind);
    newproj.dx=proj.dx(ind);
    newproj.dy=proj.dy(ind);
    newproj.depth=proj.depth(ind);
    newproj.texthdr=proj.texthdr(ind);
    newproj.texthdrfmt=proj.texthdrfmt(ind);
    newproj.segfmt=proj.segfmt(ind);
    newproj.byteorder=proj.byteorder(ind);
    newproj.binhdr=proj.binhdr(ind);
    newproj.exthdr=proj.exthdr(ind);
    newproj.tracehdr=proj.tracehdr(ind);
    newproj.kxline=proj.kxline(ind);
    newproj.isdeleted=proj.isdeleted(ind);
    newproj.deletedondisk=proj.deletedondisk(ind);
    newproj.saveneeded=zeros(size(ind));
    
    matobj=matfile([path filename],'Writable',true);
    matobj.proj=newproj;
    matobj.datasets=datasets(ind);
 
    newproj.pifigures=proj.pifigures(ind);
    newproj.datasets=datasets;
    
    %delete any gui panels
    for k=1:length(proj.gui)
        if(isgraphics(proj.gui{k}))
            delete(proj.gui{k});
        end
    end
    hmpan=findobj(hsane,'tag','master_panel');
    udat=get(hmpan,'userdata');
    geom=udat{2};
    geom(4)=geom(9);%resets the initial y coordinates of the panels
    udat{2}=geom;
    set(hmpan,'userdata',udat);
    set(hsane,'name',['SANE, Project file: ' filename]);
    set(hmsg,'string',['Project ' filename ' saved'])
    proj=newproj;
    %hmpan=findobj(hsane,'tag','master_panel');
    %udat=get(hmpan,'userdata');
    udat{1}=[];
    set(hmpan,'userdata',udat);
    
    ndatasets=length(proj.datanames);
    proj.gui=cell(1,ndatasets);
    
    % put up datapanels for each dataset, read and display if needed
    %hpanels=cell(1,ndatasets);
    %PIfigs=cell(1,ndatasets);
    for k=1:ndatasets
        if(proj.isdeleted(k)~=1)
            proj.gui{k}=newdatapanel(proj.datanames{k},proj.isloaded(k),proj.isdisplayed(k));
        end
    end
    waitsignaloff
    set(hfile,'userdata',proj);
    hpn=findobj(hsane,'tag','project_name');
    set(hpn,'string',proj.name)
elseif(strcmp(action,'loadproject'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    hmsg=findobj(hsane,'tag','message');
    [fname,path]=uigetfile('*.mat','Choose the SANE Project file (.mat) to load');
    if(fname==0)
        return
    end
    m=matfile([path fname]);
    varnames=fieldnames(m);
    ivar=[];
    for k=1:length(varnames)
        ind=strcmp(varnames{k},'proj');
        if(ind~=0)
            ivar=k;
        end
    end
    if(isempty(ivar))
        msgbox('Chosen file is not a SANE Project file, Nothing has been loaded');
        return
    end
    set(hmsg,'string',['Loading project ' path fname ' ... '])
    %check for existing project and delete any data panels
    projold=get(hfile,'userdata');
    if(~isempty(projold))
       hpanels=projold.gui;
       for k=1:length(hpanels)
           if(isgraphics(hpanels{k}))
               delete(hpanels{k});
           end
       end
       for k=1:length(projold.pifigures)
           if(isgraphics(projold.pifigures{k}))
               delete(projold.pifigures{k});
           end
       end
       hview=findobj(hsane,'tag','view');
       hk=get(hview,'children');
       delete(hk);
       hmpan=findobj(gcf,'tag','master_panel');
       udat=get(hmpan,'userdata');
       udat{1}={};
       geom=udat{2};
       geom(4)=geom(9);%resets the initial y coordinates of the panels
       udat{2}=geom;
       set(hmpan,'userdata',udat);
    end
    waitsignalon 
    proj=getfield(m,varnames{ivar});
    ndatasets=length(proj.datanames);
    proj.projfilename=fname;
    proj.projpath=path;
    proj.pifigures=cell(1,ndatasets);
    proj.gui=cell(1,ndatasets);
    proj.datasets=cell(1,ndatasets);
    set(hfile,'userdata',proj);
    loadprojectdialog
    return;
elseif(strcmp(action,'cancelprojectload'))
    hdial=gcf;
    hsane=get(hdial,'userdata');
    delete(hdial);
    hmsg=findobj(hsane,'tag','message');
    set(hmsg,'string','Project load cancelled');
    waitsignaloff
    return
elseif(strcmp(action,'loadprojdial'))
    hdial=gcf;
    hsane=get(hdial,'userdata');
    hbutt=gcbo;
    subaction=get(hbutt,'tag');
    hh=findobj(hdial,'tag','loaded');
    hloaded=get(hh,'userdata');
    hh=findobj(hdial,'tag','display');
    hdisplayed=get(hh,'userdata');
    switch subaction
        case 'allyes'
            set([hloaded hdisplayed],'value',2);
            
        case 'allno'
            set([hloaded hdisplayed],'value',1);
            
        case 'continue'
            hfile=findobj(hsane,'tag','file');
            proj=get(hfile,'userdata');
            ndata=length(hloaded);
            for k=1:ndata
                if(proj.isdeleted(k)~=1)
                    proj.isloaded(k)=get(hloaded(k),'value')-1;
                    proj.isdisplayed(k)=get(hdisplayed(k),'value')-1;
                end
            end
            set(hfile,'userdata',proj);
            delete(hdial);
            figure(hsane);
            sane('loadproject2');  
    end
    
elseif(strcmp(action,'loadproject2'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    hmsg=findobj(hsane,'tag','message');
    t0=clock;
    proj=get(hfile,'userdata');
    ndatasets=length(proj.datanames);
    m=matfile([proj.projpath proj.projfilename]);
    % put up datapanels for each dataset, read and display if needed
    %hpanels=cell(1,ndatasets);
    %PIfigs=cell(1,ndatasets);
    for k=1:ndatasets
        if(proj.isdeleted(k)~=1)
            proj.gui{k}=newdatapanel(proj.datanames{k},proj.isloaded(k),proj.isdisplayed(k));
            if(proj.isloaded(k)==1 || proj.isdisplayed(k)==1)%see if we read the dataset from disk
                cseis=m.datasets(1,k);%this reads it
                if(proj.isloaded(k)==1)
                    proj.datasets(1,k)=cseis;
                end
            end
            if(proj.isdisplayed(k)==1)%see if we display the dataset
                plotimage3D(cseis{1},proj.tcoord{k},proj.xcoord{k},proj.ycoord{k},...
                    proj.datanames{k},'seisclrs',proj.dx(k),proj.dy(k));
                set(gcf,'tag','fromsane','userdata',{k hsane});
                set(gcf,'closeRequestFcn','sane(''closepifig'');');
                hview=findobj(hsane,'tag','view');
                uimenu(hview,'label',proj.datanames{k},'callback','sane(''popupfig'');','userdata',gcf);
                proj.pifigures{k}=gcf;
                figure(hsane)
            end
        end
    end
    waitsignaloff
    set(hfile,'userdata',proj);
    hpn=findobj(hsane,'tag','project_name');
    set(hpn,'string',proj.name)
    
    tnow=clock;
    timeused=etime(tnow,t0)/60;
    if(timeused>1)
        timeused=round(100*timeused)/100;
        set(hmsg,'string',['Project ' proj.name ' loaded in ' num2str(timeused) ' min'])
    else
        timeused=round(60*10*timeused)/10;
        set(hmsg,'string',['Project ' proj.name ' loaded in ' num2str(timeused) ' sec'])
    end
    set(hsane,'name',['SANE, Project: ' proj.name ])

elseif(strcmp(action,'newproject'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    %check for existing project and delete any data panels
    projold=get(hfile,'userdata');
    if(~isempty(projold))
        hpanels=projold.gui;
        for k=1:length(hpanels)
            if(isgraphics(hpanels{k}))
                delete(hpanels{k});
            end
        end
        for k=1:length(projold.pifigures)
            if(isgraphics(projold.pifigures{k}))
                delete(projold.pifigures{k});
            end
        end
        hview=findobj(hsane,'tag','view');
        hk=get(hview,'children');
        delete(hk);
        hmpan=findobj(gcf,'tag','master_panel');
        udat=get(hmpan,'userdata');
        geom=udat{2};
        geom(4)=geom(9);%resets the initial y coordinates of the panels
        udat{1}={};
        udat{2}=geom;
        set(hmpan,'userdata',udat);
    end
    proj=makeprojectstructure;
    set(hfile,'userdata',proj);
    hpn=findobj(hsane,'tag','project_name');
    set(hpn,'string',proj.name)
    hmsg=findobj(hsane,'tag','message');
    set(hmsg,'string','Now read some data into your project')
    set(hsane,'name',['SANE, Project: ' proj.name ])
elseif(strcmp(action,'datamemory'))
    %This gets called if the "in memory" radio buttons are toggled
    %We want to be able to control whether a dataset is in memory or not. When in memory, then it is
    %present in the "datasets" field of the proj structure. If it is displayed, then it is also
    %present in the userdata of a plotimage3D window. Once a dataset is displayed, we may want to
    %clear it from memory, otherwise there are effectively two copies of it in memory. We would
    %really only want to save it in memory if we planned on applying an operation to it. If a
    %dataset has just been loaded to SEGY and not yet saved in the project, then we need to write it
    %to disk first before clearing it from memory.
    hsane=findsanefig;
    hbut=gcbo;%will be either 'Y' or 'N'
    hbg=get(hbut,'parent');
    hpan=get(hbg,'parent');
    hmsg=findobj(hsane,'tag','message');
    idata=hpan.UserData;%the dataset number
    choice=get(hbut,'string');%this will be 'Y' or 'N'
    hfile=findobj(hsane,'tag','file');
    proj=hfile.UserData;
    if(strcmp(choice','N'))
        %dataset is being cleared from memory
        %first check to ensure that the project has been save so that a project file exists on disk
        if(isempty(proj.projfilename)||isempty(proj.projpath))
            msgbox('Please save the project to disk before clearing data from memory');
            return;
        end
        %we also close any display
%         hfig=proj.pifigures(idata);
%         if(isgraphics(hfig))
%             close(hfig);
%         end
%         proj.pifigures{idata}=[];
%         proj.isdisplayed(idata)=0;
        proj.isloaded(idata)=0;
        %now we need to be sure that the data exists on disk before clearing it from memory
        if(proj.saveneeded(idata)==1)
            %open the project file
            matobj=matfile([proj.projpath proj.projfilename],'writable',true);
            [meh,ndatasetsondisk]=size(matobj,'datasets'); %#ok<ASGLU>%this is the number datasets that exist on disk
            ndatasets=length(proj.datasets);%this is how many datasets there are in total
            if(ndatasetsondisk<ndatasets)
                nnew=ndatasetsondisk+1:ndatasets;
                matobj.datasets(1,nnew)=cell(1,length(nnew));
            end
            if(isempty(proj.datasets(1,idata)))
                error('attempt to save empty dataset');
            end
            matobj.datasets(1,idata)=proj.datasets(1,idata);
            proj.saveneeded(idata)=0;
            proj.datasets{idata}=[];
            pifigs=proj.pifigures;
            proj.pifigures=[];
            matobj.proj=proj;%need to save the project for consistency on disk
            proj.pifigures=pifigs;
        end
        proj.datasets{idata}=[];
%         %set the isdisplayed button to no
%         hno=findobj(hpan,'tag','displayno');
%         set(hno,'value',0);
        set(hfile,'userdata',proj);
        set(hmsg,'string',['Datset ' proj.datanames{idata} ' cleared from memory but may still be displayed']);
    else
        %dataset is being loaded into memory and displayed
        hmpan=findobj(gcf,'tag','master_panel');
        udat=hmpan.UserData;
        udat{3}=idata;%this flags to reload which dataset we are reading
        hmpan.UserData=udat;
        sane('reloaddataset');%this updates proj and displays the dataset
        %set the isdisplayed button to yes
        hyes=findobj(hpan,'tag','displayyes');
        set(hyes,'value',1);
    end
elseif(strcmp(action,'datadisplay'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    hmsg=findobj(hsane,'tag','message');
    proj=get(hfile,'userdata');
    %determine choice
    hbut=gcbo;
    butname=get(hbut,'string');
    val=get(hbut,'value');
    if((strcmp(butname,'Y')&&val==1)||(strcmp(butname,'N')&&val==0))
        choice='display';
    else
        choice='dontdisplay';
    end
    %deternime dataset number
    hpan=get(get(hbut,'parent'),'parent');
    idata=get(hpan,'userdata');
    
    switch choice
        case 'display'
            waitsignalon
            %check if dataset is loaded
            if(~proj.isloaded(idata)||isempty(proj.datasets{idata}))
                hmpan=findobj(gcf,'tag','master_panel');
                udat=hmpan.UserData;
                udat{3}=idata;%this flags to reload which dataset we are reading
                hmpan.UserData=udat;
                sane('reloaddataset');%this loads and displays and updates proj
                %set the isdisplayed button to yes
                hyes=findobj(hpan,'tag','displayyes');
                set(hyes,'value',1);
                hyes=findobj(hpan,'tag','memoryyes');
                set(hyes,'value',1);
            else
                %now display
                plotimage3D(proj.datasets{idata},proj.tcoord{idata},proj.xcoord{idata},...
                    proj.ycoord{idata},proj.datanames{idata},'seisclrs',proj.dx(idata),proj.dy(idata))
                set(gcf,'tag','fromsane','userdata',{idata hsane});
                set(gcf,'closeRequestFcn','sane(''closepifig'');')
                hview=findobj(hsane,'tag','view');
                uimenu(hview,'label',proj.datanames{idata},'callback','sane(''popupfig'');','userdata',gcf);
                proj.isdisplayed(idata)=1;
                proj.pifigures{idata}=gcf;
                set(hfile,'userdata',proj);
                set(hmsg,'string',['Dataset ' ,proj.datanames{idata} ' displayed'])
            end
            waitsignaloff
        case 'dontdisplay'
            hpifig=proj.pifigures{idata};
            if(isgraphics(hpifig))
                figure(hpifig);
                sane('closepifig');
            else
                proj.isdisplayed(idata)=0;
                proj.pifigures{idata}=[];
                set(hfile,'userdata',proj)
            end
            
    end
    figure(hsane)
elseif(strcmp(action,'closepifig'))
    hthisfig=gcbf;
    if(strcmp(get(hthisfig,'tag'),'sane'))
        %ok, the call came from SANE
        hpan=get(get(gcbo,'parent'),'parent');%should be the panel of the dataset whose figure is closing
        idat=get(hpan,'userdata');
        hfile=findobj(hthisfig,'tag','file');
        proj=get(hfile,'userdata');
        hthisfig=proj.pifigures{idat};
    end
    test=get(hthisfig,'name');
    %ind=strfind(test,'plotimage3D');
    tag=get(hthisfig,'tag');
    udat=get(hthisfig,'userdata');
    if(length(udat)>1)
        hsane=udat{2};
    else
        hsane=findsanefig;
    end
    hview=findobj(hsane,'tag','view');
    if(contains(test,'plotimage3D')&&strcmp(tag,'fromsane')&&iscell(udat)&&length(udat)==2)
        hmenus=get(hview,'children');
        for k=1:length(hmenus)
           fig=get(hmenus(k),'userdata');
           if(fig==hthisfig)
              delete(hmenus(k));
           end
        end
%         PLOTIMAGE3DTHISFIG=proj.pifigures{idata};
%         flag=get(hg,'value');
%         if(flag==1)
%             plotimage3D('groupex');
%         else
%             plotimage3D('ungroupex')
%         end
        idata=udat{1};
        plotimage3D('closesane');
        if(isgraphics(hthisfig))
            return; %this happens if they choose not to close
        end
        %delete(hthisfig);
        hmsg=findobj(hsane,'tag','message');
        
        hfile=findobj(hsane,'tag','file');
        proj=get(hfile,'userdata');
        proj.pifigures{idata}=[];
        hpanels=proj.gui;
        hpan=hpanels{idata};
        hno=findobj(hpan,'tag','displayno');
        set(hno,'value',1);
        proj.isdisplayed(idata)=0;
        %ungroup if needed
        hg=findobj(hpan,'tag','group');
        val=get(hg,'value');
        if(val==1)
            set(hg,'value',0);
        end
        set(hfile,'userdata',proj);
        bn=questdlg('Do you want to keep the data in memory?','Memory question','Yes','No','No');
        if(strcmp(bn,'No'))
            %remove from memory
            waitsignalon
            %need to check if it has been saved before we delete it
            if(~isempty(proj.projfilename))
                %this means the project has been save at least onece. However, the dataset might not
                %have been. So, we check for dataset saved.
                mObj=matfile([proj.projpath proj.projfilename],'writable',true);
                [meh,ndatasetsondisk]=size(mObj,'datasets'); %#ok<ASGLU>
                ndatasets=length(proj.datasets);
                if(ndatasets>ndatasetsondisk)
                    jdata=ndatasetsondisk+1:ndatasets;
                    nnew=length(jdata);
                    mObj.datasets(1,jdata)=cell(1,nnew);
                end
                if(isempty(mObj.datasets(1,idata)))
                    set(hmsg,'string','Saving dataset to disk');
                    if(isempty(proj.datasets(1,idata)))
                        error('attempt to save empty dataset');
                    end
                    mObj.datasets(1,idata)=proj.datasets(1,idata);
                end
                waitsignaloff
            else
                sane('saveproject');
                proj=get(hfile,'userdata');
            end
            proj.datasets{idata}=[];
            hno=findobj(hpan,'tag','memoryno');
            set(hno,'value',1);
            proj.isloaded(idata)=0;
            set(hmsg,'string','Dataset removed from memory')
            set(hfile,'userdata',proj);
        else
            set(hmsg,'string','Display closed but data retained in memory');
        end
        
        figure(hsane)
    end
elseif(strcmp(action,'datanamechange'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    hmsg=findobj(hsane,'tag','message');
    hnamebox=gcbo;
    idata=get(get(hnamebox,'parent'),'userdata');
    oldname=proj.datanames{idata};
    newname=get(hnamebox,'string');
    proj.datanames{idata}=newname;
    set(hfile,'userdata',proj);
    %look for open pi3d windows and change those
    hview=findobj(hsane,'tag','view');
    hkids=get(hview,'children');
    for k=1:length(hkids)
        if(strcmp(get(hkids(k),'label'),oldname))%see if the name matches
            hpifig=get(hkids(k),'userdata');
            if(isgraphics(hpifig))
                udat=get(hpifig,'userdata');
                jdata=udat{1};%there might be two datasets with the same name so jdata must match idata
                if(jdata==idata)
                    %this is it
                    plotimage3D('datanamechange',hpifig,newname);
                    set(hkids(k),'label',newname);
                end
            end
        end
    end
    set(hmsg,'string',[' dataset name changed to ' newname])
elseif(strcmp(action,'close'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    hmsg=findobj(hsane,'tag','message');
    bn=questdlg('Do you want to save the project first?','Close SANE','Yes','No','Cancel','Yes');
    switch bn
        case 'Cancel'
            set(hmsg,'string','Close cancelled')
            return;
        case 'Yes'
            sane('saveproject')
            delete(hsane);
            for k=1:length(proj.pifigures)
               if(isgraphics(proj.pifigures{k}))
                   delete(proj.pifigures{k});
               end
            end
        case 'No'
            delete(hsane);
            for k=1:length(proj.pifigures)
               if(isgraphics(proj.pifigures{k}))
                   delete(proj.pifigures{k});
               end
            end
    end
elseif(strcmp(action,'popupfig'))
    hmenu=gcbo;
    fig=get(hmenu,'userdata');
    if(isgraphics(fig))
        figure(fig);
    end
elseif(strcmp(action,'datadelete'))
    hsane=findsanefig;
    hbutt=gcbo;
    idata=get(get(hbutt,'parent'),'userdata');
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    hmsg=findobj(hsane,'tag','message');
    %confirm
    choice=questdlg(['Please confirm the deletion of dataset ' proj.datanames{idata}],'Data deletion','Yes','No','Cancel','Yes');
    switch choice
        case 'No'
            set(hmsg,'string','Data deletion cancelled');
        case 'Cancel'
            set(hmsg,'string','Data deletion cancelled');
        case 'Yes'
            proj.isdeleted(idata)=1;%set deletion flag
            proj.deletedondisk(idata)=0;
            if(isgraphics(proj.pifigures{idata}))
                delete(proj.pifigures{idata});
            end
            proj.filenames{idata}=[];
            proj.paths{idata}=[];
            deadname=proj.datanames{idata};
            proj.datanames{idata}=[];
            proj.isloaded(idata)=0;
            proj.isdisplayed(idata)=0;
            proj.xcoord{idata}=[];
            proj.ycoord{idata}=[];
            proj.tcoord{idata}=[];
            proj.tshift(idata)=0;
            proj.datasets{idata}=[];
            proj.xcdp{idata}=[];
            proj.ycdp{idata}=[];
            proj.dx(idata)=0;
            proj.dy(idata)=0;
            proj.depth(idata)=0;
            proj.texthdr{idata}=[];
            proj.texthdrfmt{idata}=[];
            proj.segfmt{idata}=[];
            proj.byteorder{idata}=[];
            proj.binhdr{idata}=[];
            proj.exthdr{idata}=[];
            proj.tracehdr{idata}=[];
            proj.kxline{idata}=[];
            if(isgraphics(proj.gui{idata}))
                pos=get(proj.gui{idata},'position');
                delete(proj.gui{idata});
                uicontrol(hsane,'style','text','string',{['dataset ' deadname ' has been deleted.'],...
                    'This space will disappear when you save and reload the Project. Deletion on disk does not happen until you save the project'},...
                    'units','normalized','position',pos);
            end
            set(hmsg,'string',['dataset ' deadname ' has been deleted.'])
            set(hfile,'userdata',proj);
    end
elseif(strcmp(action,'group'))
    hsane=findsanefig;
    hg=gcbo;
    idata=get(get(hg,'parent'),'userdata');
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    if(proj.isdisplayed(idata)==0)
        msgbox({'Grouping only works if the dataset is displayed.' 'Display it first, then group.'});
        return
    end
    if(isgraphics(proj.pifigures{idata}))
        PLOTIMAGE3DTHISFIG=proj.pifigures{idata};
        flag=get(hg,'value');
        if(flag==1)
            plotimage3D('groupex');
        else
            plotimage3D('ungroupex')
        end
    else
        error('Logic error in SANE when trying to group/ungroup figures')
    end
elseif(strcmp(action,'pi3d:group'))
    %this is a message from plotimage3D saying either a group or an ungroup has happened
%     hpi3d=gcf;%if this happens then a pi3d figure is current
%     %now we verify that the pi3d figure is from sane
%     udat=get(hpi3D,'userdata');
%     tag=get(hpi3d,'tag');
%     if(strcmp(tag,'fromsane')&&length(udat)==2)
%         hsane=udat{2};
%     else
%         return; %if this happens then the logic has failed
%     end
    hsane=arg2{2};
    
    hgroup=PLOTIMAGE3DFIGS;%these are the grouped figures
    
    hmpan=findobj(hsane,'tag','master_panel');%the master panel is the key to SANE data panels
    udat2=get(hmpan,'userdata');
    hpanels=udat2{1};%the sane data panels
    idatas=zeros(size(hgroup));
    %loop over the grouped figures and find their sane data numbers
    for k=1:length(hgroup)
        udat3=get(hgroup(k),'userdata');
        if(strcmp(get(hgroup(k),'tag'),'fromsane'))
            idatas(k)=udat3{1};
        end    
    end
    %now loop over panels and compare their data numbers to those in the group
    for k=1:length(hpanels)
        hg=findobj(hpanels{k},'tag','group');
        idata=get(hpanels{k},'userdata');
        ind=find(idata==idatas, 1);
        if(~isempty(ind))
            set(hg,'value',1);
        else
            set(hg,'value',0);
        end
    end
elseif(strcmp(action,'writesegy'))
    hsane=findsanefig;
    hmsg=findobj(hsane,'tag','message');
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    iout=listdlg('Promptstring','Choose the dataset for output','liststring',proj.datanames,...
        'selectionmode','single','listsize',[500,300]);
    if(isempty(iout))
        return;
    end
    [fname,path]=uiputfile('*.sgy','Select the output file');
    if(isequal(fname,0) || isequal(path,0))
        msgbox('Output cancelled');
        return;
    end
    if(exist([path fname],'file'))
        delete([path fname]);
    end
    if(isempty(proj.datasets{iout}))
        %recall the dataset from disk
        mObj=matfile([pro.projpath proj.projfilename]);
        cseis=mObj.datasets(1,iout);
    else
        cseis=proj.datasets(1,iout);
    end
    waitsignalon
    set(hmsg,'string','Forming output array');
    seis=unmake3Dvol(cseis{1},proj.xcoord{iout},proj.ycoord{iout},proj.xcdp{iout},proj.ycdp{iout},...
        'kxlineall',proj.kxline{iout});
    dt=abs(proj.tcoord{iout}(2)-proj.tcoord{iout}(1));
    set(hmsg,'string','Beginning SEGY output');
    writesegy([path fname],seis,proj.segyrev(iout),dt,proj.segfmt{iout},proj.texthdrfmt{iout},...
        proj.byteorder{iout},proj.texthdr{iout},proj.binhdr{iout},proj.exthdr{iout},...
        proj.tracehdr{iout},proj.bindef{iout},proj.trcdef{iout},hsane);
    
    waitsignaloff
    set(hmsg,'string',['Dataset ' [path fname] ' written']);
elseif(strcmp(action,'writemat'))
    msgbox('Sorry, feature not yet implemented')
    return;
elseif(strcmp(action,'starttask'))
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    proj=hfile.UserData;
    if(isempty(proj.datasets))
        msgbox('You need to load some data before you can do this!','Oh oh ...');
        return
    end
    %determine the task
    task=get(gcbo,'tag');
    parmset=getparmset(task);
    switch task
        case 'filter'
            sanetask(proj.datanames,parmset,task);
        case 'phasemap'
            return;
        case 'spikingdecon'
            sanetask(proj.datanames,parmset,task);
            return;
        case 'fdom'
            sanetask(proj.datanames,parmset,task);
            return;
        case 'wavenumber'
            sanetask(proj.datanames,parmset,task);
    end
elseif(strcmp(action,'dotask'))
    htaskfig=gcf;
    htask=findobj(htaskfig,'tag','task');
    udat=get(htask,'userdata');
    task=udat{1};
    parmset=udat{2};%these are the parameters before user modification
    nparms=(length(parmset)-1)/3;
    %determine the dataset name, get the project
    hdat=findobj(htaskfig,'tag','datasets');
    idat=get(hdat,'value');
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    hmsg=findobj(hsane,'tag','message');
    proj=hfile.UserData;
    %pull the updated parameters out of the gui
    for k=1:nparms
       hobj=findobj(htaskfig,'tag',parmset{3*(k-1)+2});
       val=get(hobj,'string');
       parm=parmset{3*(k-1)+3};
       if(iscell(val))
           ival=get(hobj,'value');
           parm{end}=ival;
       else
           parm=val;
       end
       parmset{3*(k-1)+3}=parm;
    end
    %check the parms for validity
    switch task
        case 'filter'
            parmset=parmsetfilter(parmset);
            %if parmset comes back as a string, then we have failure
        case 'phasemap'
            
            
        case 'spikingdecon'
            parmset=parmsetdecon(parmset,proj.tcoord{idat});
            
        case 'fdom'
            parmset=parmsetfdom(parmset,proj.tcoord{idat});
            
        case 'wavenumber'
           parmset=parmsetwavenumber(parmset); 
            
    end
    if(ischar(parmset))
        msgbox(parmset,'Oh oh, there are problems...');
        return;
    end
    %save the updated parmset
    setparmset(parmset);
    %get the dataset
    if(~proj.isloaded(idat))
        %load the dataset
        figure(hsane)
        waitsignalon
        set(hmsg,'string',['Loading ' proj.datanames{idat} ' from disk']);
        hmpan=findobj(hsane,'tag','master_panel');
        udat=hmpan.UserData;
        udat{3}=idat;%this flags to reload which dataset we are reading
        hmpan.UserData=udat;
        sane('reloaddataset');%this updates proj and displays the dataset
        proj=hfile.UserData;
        waitsignaloff
    end
    %determine output dataset's fate
    hout=findobj(htaskfig,'tag','outputs');
    outopts=get(hout,'string');
    fate=outopts{get(hout,'value')};
    seis=proj.datasets{idat};
    t=proj.tcoord{idat};
    x=proj.xcoord{idat};
    y=proj.ycoord{idat};
    dx=proj.dx(idat);
    dy=proj.dy(idat);
    dname=proj.datanames{idat};
    %close the task window
    close(htaskfig);
    %start the task
    hcompute=findobj(hsane,'label','Compute');
    switch task
        case 'filter'
            set(hcompute,'userdata',idat); %this allows a "cancel" to determine the input dataset
            set(hmsg,'string','Bandpass filtering in progress')
            fmin=getparm(parmset,'fmin');
            dfmin=getparm(parmset,'dfmin');
            fmax=getparm(parmset,'fmax');
            dfmax=getparm(parmset,'dfmax');
            phase=getparm(parmset,'phase');
            if(strcmp(phase,'zero'))
                phase=0;
            else
                phase=1;
            end
            nx=length(x);ny=length(y);
            itrace=0;
            %The waitbar implements a cancel operation that works through the globals HWAIT and
            %CONTINUE. When the cancel button is hit (on the waitbar) the callback 'sane(''canceltaskafter'')'
            %sets the value of CONTINUE to false which causes the loop to stop. The callback also
            %removes the input dataset from memory and deletes the waitbar. Thus to restart the task
            %(it always must start from the beginning) then the dataset must be reloaded.
            HWAIT=waitbar(0,'Please wait for bandpass filtering to complete','CreateCancelBtn','sane(''canceltaskafter'')');
            ntraces=nx*ny;
            t0=clock;
            CONTINUE=true;
            for k=1:nx
                for j=1:ny
                    tmp=seis(:,k,j);
                    if(sum(abs(tmp))>0)
                        seis(:,k,j)=filtf(tmp,t,[fmin dfmin],[fmax dfmax],phase);
                    end
                    itrace=itrace+1;
                    if(~CONTINUE)
                        break;
                    end
                end
                if(~CONTINUE)
                        break;
                end
                t1=clock;
                timeused=etime(t1,t0);
                timeleft=(timeused/itrace)*(ntraces-itrace)/60;%in minutes
                timeleft=round(timeleft*100)/100;
                waitbar(itrace/ntraces,HWAIT,['Estimated time remaining ' num2str(timeleft) ' minutes']);
            end
            t1=clock;
            timeused=etime(t1,t0);
            set(hmsg,'string',['Completed task ' task ' for ' dname ' in ' num2str(timeused/60) ' minutes'])
            delete(HWAIT)
            set(hcompute,'userdata',[]);
        case 'spikingdecon'
            set(hcompute,'userdata',idat); %this allows a "cancel" to determine the input dataset
            set(hmsg,'string','Spiking decon in progress')
            oplen=getparm(parmset,'oplen');
            stab=getparm(parmset,'stab');
            topgate=getparm(parmset,'topgate');
            botgate=getparm(parmset,'botgate');
            fmin=getparm(parmset,'fmin');
            dfmin=getparm(parmset,'dfmin');
            fmax=getparm(parmset,'fmax');
            dfmax=getparm(parmset,'dfmax');
            phase=getparm(parmset,'phase');
            dt=t(2)-t(1);
            if(strcmp(phase,'zero'))
                phase=0;
            else
                phase=1;
            end
            nx=length(x);ny=length(y);
            itrace=0;
            HWAIT=waitbar(0,'Please wait for decon to complete','CreateCancelBtn','sane(''canceltaskafter'')');
            ntraces=nx*ny;
            t0=clock;
            idesign=near(t,topgate,botgate);
            nop=round(oplen/dt);
            CONTINUE=true;
            for k=1:nx
                for j=1:ny
                    tmp=seis(:,k,j);
                    tmpd=seis(idesign,k,j);
                    if(sum(abs(tmpd))>0)
                        tmpdecon=deconw(tmp,tmpd,nop,stab);
                        seis(:,k,j)=filtf(tmpdecon,t,[fmin dfmin],[fmax dfmax],phase);
                    end
                    itrace=itrace+1;
                    if(~CONTINUE)
                        break;
                    end
                end
                if(~CONTINUE)
                        break;
                end
                t1=clock;
                timeused=etime(t1,t0);
                timeleft=(timeused/itrace)*(ntraces-itrace)/60;%in minutes
                timeleft=round(timeleft*100)/100;
                waitbar(itrace/ntraces,HWAIT,['Estimated time remaining ' num2str(timeleft) ' minutes']);
            end
            t1=clock;
            timeused=etime(t1,t0);
            set(hmsg,'string',['Completed task ' task ' for ' dname ' in ' num2str(timeused/60) ' minutes'])
            delete(HWAIT)
            set(hcompute,'userdata',[]);
        case 'wavenumber'
            set(hcompute,'userdata',idat); %this allows a "cancel" to determine the input dataset
            set(hmsg,'string','Wavenumber filtering in progress')
            sigmax=getparm(parmset,'sigmax');
            sigmay=getparm(parmset,'sigmay');
            
            nt=length(t);
            HWAIT=waitbar(0,'Please wait for wavenumber filtering to complete','CreateCancelBtn','sane(''canceltaskafter'')');
            t0=clock;
            ievery=10;
            CONTINUE=true;
            for k=1:nt

                slice=squeeze(seis(k,:,:));
                slice2=wavenumber_gaussmask2(slice,sigmax,sigmay);
                slice2=slice2*norm(slice)/norm(slice2);
                seis(k,:,:)=shiftdim(slice2,-1);
                if(~CONTINUE)
                        break;
                end
                if(rem(k,ievery)==0)
                    tnow=clock;
                    timeused=etime(tnow,t0);
                    timeperslice=timeused/k;
                    timeleft=timeperslice*(nt-k)/60;
                    timeleft=round(timeleft*100)/100;
                    waitbar(k/nt,HWAIT,['Estimated time remaining ' num2str(timeleft) ' minutes']);
                end

            end
            t1=clock;
            timeused=etime(t1,t0);
            set(hmsg,'string',['Completed task ' task ' for ' dname ' in ' num2str(timeused/60) ' minutes'])
            delete(HWAIT)
            set(hcompute,'userdata',[]);
        case 'phasemap'
            
            
        case 'fdom'
            ny=length(y);
            twin=getparm(parmset,'twin');
            ninc=getparm(parmset,'ninc');
            fmax=getparm(parmset,'Fmax');
            tfmax=getparm(parmset,'tfmax');
            tinc=ninc*(t(2)-t(1));
            interpflag=1;
            p=2;
            fc=1;
            
            
            if(isnan(tfmax))
                fmt0=fmax;
            else
                fmt0=[fmax tfmax];
            end
            
            set(hcompute,'userdata',idat); %this allows a "cancel" to determine the input dataset
            set(hmsg,'string','Dominant frequency volume computation in progress')
            
            HWAIT=waitbar(0,'Please wait for fdom computation to complete','CreateCancelBtn','sane(''canceltaskafter'')');
            t0=clock;
            CONTINUE=true;
            ievery=1;
            for k=1:ny
                %process each iline as a panel in tvfdom3
                spanel=squeeze(seis(:,:,k));
                test=sum(abs(spanel));
                ilive=find(test~=0);
                if(~isempty(ilive))
                    fd=tvfdom3(spanel(:,ilive),t,twin,tinc,fmt0,interpflag,p,fc);
                    ind=find(fd<0);
                    if(~isempty(ind))
                        fd(ind)=0;
                    end
                end
                seis(:,ilive,k)=single(fd);
                if(rem(k,ievery)==0)
                    time_used=etime(clock,t0);
                    time_per_line=time_used/k;
                    timeleft=(ny-k-1)*time_per_line/60;
                    timeleft=round(100*timeleft)/100;
                    waitbar(k/ny,HWAIT,['Estimated time remaining ' num2str(timeleft) ' minutes']);
                end
            end
            t1=clock;
            timeused=etime(t1,t0)/60;
            timeused=round(timeused*100)/100;
            set(hmsg,'string',['Completed task ' task ' for ' dname ' in ' num2str(timeused) ' minutes'])
            delete(HWAIT)
            set(hcompute,'userdata',[]);
    end
    if(~CONTINUE)
        waitsignaloff
        set(hmsg,'string','Computation interrupted by user, input dataset unloaded');
        return
    end
    %deal with the output
    switch fate
        case 'Save SEGY'
            [fname,path]=uiputfile('*.sgy','Select the output file');
            if(isequal(fname,0) || isequal(path,0))
                msgbox('Output cancelled');
                return;
            end
            if(exist([path fname],'file'))
                delete([path fname]);
            end
            waitsignalon
            set(hmsg,'string','Forming output array');
            seis=unmake3Dvol(seis,proj.xcoord{idat},proj.ycoord{idat},proj.xcdp{idat},proj.ycdp{idat},...
                'kxlineall',proj.kxline{idat});
            dt=abs(proj.tcoord{idat}(2)-proj.tcoord{idat}(1));
            set(hmsg,'string','Beginning SEGY output');
            writesegy([path fname],seis,proj.segyrev(idat),dt,proj.segfmt{idat},proj.texthdrfmt{idat},...
                proj.byteorder{idat},proj.texthdr{idat},proj.binhdr{idat},proj.exthdr{idat},...
                proj.tracehdr{idat},proj.bindef{idat},proj.trcdef{idat},hsane);
            
            waitsignaloff
            set(hmsg,'string',['Dataset ' [path fname] ' written']);
        case 'Save SEGY and display'
            %display
            plotimage3D(seis,t,x,y,[dname ' ' task],dx,dy)
            %output
            [fname,path]=uiputfile('*.sgy','Select the output file');
            if(isequal(fname,0) || isequal(path,0))
                msgbox('Output cancelled');
                return;
            end
            if(exist([path fname],'file'))
                delete([path fname]);
            end
            waitsignalon
            set(hmsg,'string','Forming output array');
            seis=unmake3Dvol(seis,proj.xcoord{idat},proj.ycoord{idat},proj.xcdp{idat},proj.ycdp{idat},...
                'kxlineall',proj.kxline{idat});
            dt=abs(proj.tcoord{idat}(2)-proj.tcoord{idat}(1));
            set(hmsg,'string','Beginning SEGY output');
            writesegy([path fname],seis,proj.segyrev(idat),dt,proj.segfmt{idat},proj.texthdrfmt{idat},...
                proj.byteorder{idat},proj.texthdr{idat},proj.binhdr{idat},proj.exthdr{idat},...
                proj.tracehdr{idat},proj.bindef{idat},proj.trcdef{idat},hsane);
            
            waitsignaloff
            set(hmsg,'string',['Dataset ' [path fname] ' written']);
        case 'Replace input in project'
            proj.datasets{idat}=seis;
            %see if data is displayed
            if(proj.isdisplayed{idat}==1)
                delete(proj.pifigures{idat})
                plotimage3D(seis,t,x,y,dname,'seisclrs',proj.dx(idat),proj.dy(idat));
                set(gcf,'tag','fromsane','userdata',{idat hsane});
                set(gcf,'closeRequestFcn','sane(''closepifig'');')
                hview=findobj(hsane,'tag','view');
                uimenu(hview,'label',dname,'callback','sane(''popupfig'');','userdata',gcf);
                proj.pifigures{idat}=gcf;
                proj.isdisplayed(idat)=1;
            end
            proj.saveneeded(idat)=1;
            set(hfile,'userdata',proj)
            set(hmsg,'string',[proj.datanames{idat} ' replaced. Data will be written to disk when you save the project.'])
        case 'Save in project as new'
            waitsignalon
            %update project structure
            ndatasets=length(proj.datanames)+1;
            proj.filenames{ndatasets}=proj.filenames{idat};
            proj.paths{ndatasets}=proj.paths{idat};
            proj.datanames{ndatasets}=[proj.datanames{idat} '_' task];
            proj.isloaded(ndatasets)=1;
            proj.isdisplayed(ndatasets)=1;
            proj.xcoord{ndatasets}=proj.xcoord{idat};
            proj.ycoord{ndatasets}=proj.ycoord{idat};
            proj.tcoord{ndatasets}=proj.tcoord{idat};
            proj.datasets{ndatasets}=seis;
            proj.tshift(ndatasets)=proj.tshift(idat);
            proj.xcdp{ndatasets}=proj.xcdp{idat};
            proj.ycdp{ndatasets}=proj.ycdp{idat};
            proj.dx(ndatasets)=proj.dx(idat);
            proj.dy(ndatasets)=proj.dy(idat);
            proj.depth(ndatasets)=proj.depth(idat);
            proj.texthdr{ndatasets}=proj.texthdr{idat};
            proj.texthdrfmt{ndatasets}=proj.texthdrfmt{idat};
            proj.segfmt{ndatasets}=proj.segfmt{idat};
            proj.byteorder{ndatasets}=proj.byteorder{idat};
            proj.binhdr{ndatasets}=proj.binhdr{idat};
            proj.exthdr{ndatasets}=proj.exthdr{idat};
            proj.tracehdr{ndatasets}=proj.tracehdr{idat};
            proj.bindef{ndatasets}=proj.bindef{idat};
            proj.trcdef{ndatasets}=proj.trcdef{idat};
            proj.kxline{ndatasets}=proj.kxline{idat};
            
            proj.segyrev(ndatasets)=proj.segyrev(idat);
            proj.pifigures{ndatasets}=[];
            proj.isdeleted(ndatasets)=0;
            proj.deletedondisk(ndatasets)=0;
            proj.saveneeded(ndatasets)=1;
            
            hpan=newdatapanel(proj.datanames{ndatasets},1,1);
            proj.gui{ndatasets}=hpan;
            
            %call plotimage3D
            plotimage3D(seis,t,x,y,proj.datanames{ndatasets},'seisclrs',dx,dy);
            set(gcf,'tag','fromsane','userdata',{ndatasets hsane});
            set(gcf,'closeRequestFcn','sane(''closepifig'');')
            hview=findobj(hsane,'tag','view');
            uimenu(hview,'label',dname,'callback','sane(''popupfig'');','userdata',gcf);
            proj.pifigures{ndatasets}=gcf;
            
            %save the project
            set(hfile,'userdata',proj);
            figure(hsane)
            waitsignaloff
            set(hmsg,'string',[proj.datanames{ndatasets} ' saved and displayed. Data will be written to disk when you save the project.'])
    end
    
    
    
elseif(strcmp(action,'canceltask'))
    %this is called if the task is cancelled before it actually begins to compute
    delete(gcf);
    hsane=findsanefig;
    figure(hsane);
    hmsg=findobj(hsane,'tag','message');
    set(hmsg,'string','Computation task cancelled');
    return
    
    
elseif(strcmp(action,'canceltaskafter'))
    %this is call if the task has already started to compute. This is a special problem because the
    %tasks are designed to replace the input data array either trace-by-trace or slice-by-slice. So
    %cancelling the task once the computation has begun means the data volume must be discarded.
    drawnow
    hsane=findsanefig;
    hmsg=findobj(hsane,'tag','message');
    set(hmsg,'string','Computation task cancelled');
    CONTINUE=false;
    delete(HWAIT);
    hcompute=findobj(hsane,'label','Compute');
    idat=get(hcompute,'userdata');
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    proj.isloaded(idat)=0;%this will cause a reload if the task is restarted
    memorybuttonoff(idat);
    proj.datasets{idat}=[];
    set(hfile,'userdata',proj);
    return
end

end

function [itchoice,ixchoice,iychoice]=ambigdialog(ambig,it,ix,iin,varnames,itchoice,ixchoice,iychoice)
%this is called when importing a .mat file that may have many variables in it.
hsane=findsanefig;
pos=get(hsane,'position');
if(pos(3)<450);pos(3)=450;end
%hdial=dialog;
hdial=figure('windowstyle','modal');
set(hdial,'position',pos,'menubar','none','toolbar','none','numbertitle','off',...
    'name','SANE mat file input ambiguity dialog','nextplot','new');
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
    'position',[.1 .1 .3 .1],'callback',@checkambig,'userdata',ambig,'tag','done');

uiwait(hdial)

%function [itchoice,ixchoice,iychoice]=checkambig(~,~)
function checkambig(~,~)
htable=findobj(gcf,'tag','table');
choices=htable.Data;
if(length(choices)==3)
    if(strcmp(choices{1},choices{2})||strcmp(choices{1},choices{3})...
            ||strcmp(choices{2},choices{3}))
        hmsg=findobj(gcf,'tag','msg');
        set(hmsg,'string','You must choose unique names for each!!!','foregroundcolor',[1 0 0],...
            'fontsize',10,'fontweight','bold');
        itchoice=0;
        ixchoice=0;
        iychoice=0;
        return
    end
elseif(length(choices)==2)
    if(strcmp(choices{1},choices{2}))
        hmsg=findobj(gcf,'tag','msg');
        set(hmsg,'string','You must choose unique names for each!!!','foregroundcolor',[1 0 0],...
            'fontsize',10,'fontweight','bold');
        itchoice=0;
        ixchoice=0;
        iychoice=0;
        return
    end
end
varnames=get(htable,'userdata');
hbut=findobj(gcf,'tag','done');
ambig=get(hbut,'userdata');
colnames=htable.ColumnName;
if(sum(ambig)==3)
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
elseif(sum(ambig)==2)
    for k=1:length(varnames)
        if(strcmp(varnames{k},choices{1}))
            if(colnames{1}(1)=='t')
                itchoice=k;
            elseif(colnames{1}(1)=='x')
                ixchoice=k;
            end
        end
        if(strcmp(varnames{k},choices{2}))
            if(colnames{2}(1)=='x')
                ixchoice=k;
            elseif(colnames{2}(1)=='i')
                iychoice=k;
            end
        end
    end
else
    for k=1:length(varnames)
        if(strcmp(varnames{k},choices{1}))
            if(colnames{1}(1)=='t')
                itchoice=k;
            elseif(colnames{1}(1)=='x')
                ixchoice=k;
            elseif(colnames{1}(1)=='i')
                iychoice=k;
            end
        end
    end
end
close(gcf)
end

end

function proj=makeprojectstructure

proj.name='New Project';
proj.filenames={};
proj.projfilename=[];
proj.projpath=[];
proj.paths={};
proj.datanames={};
proj.isloaded=[];
proj.isdisplayed=[];
proj.xcoord={};
proj.ycoord={};
proj.tcoord={};
proj.tshift=[];
proj.datasets={};
proj.xcdp={};
proj.ycdp={};
proj.dx=[];
proj.dy=[];
proj.depth=[];
proj.texthdr={};
proj.texthdrfmt={};
proj.segfmt={};
proj.byteorder={};
proj.binhdr={};
proj.exthdr={};
proj.tracehdr={};
proj.bindef={};
proj.trcdef={};
proj.segyrev={};
proj.kxline={};
proj.gui={};
proj.rspath=[];
proj.wspath=[];
proj.rmpath=[];
proj.wmpath=[];
proj.pifigures=[];
proj.isdeleted=[];
proj.deletedondisk=[];
proj.saveneeded=[];

end

function projnew=expandprojectstructure(proj,nnew)

projnew.name=proj.name;
projnew.projfilename=proj.projfilename;
projnew.filenames=[proj.filenames cell(1,nnew)];
projnew.projpath=proj.projpath;
projnew.paths=[proj.paths cell(1,nnew)];
projnew.datanames=[proj.datanames cell(1,nnew)];
projnew.isloaded=[proj.isloaded zeros(1,nnew)];
projnew.isdisplayed=[proj.isdisplayed zeros(1,nnew)];
projnew.xcoord=[proj.xcoord cell(1,nnew)];
projnew.ycoord=[proj.ycoord cell(1,nnew)];
projnew.tcoord=[proj.tcoord cell(1,nnew)];
projnew.tshift=[proj.tshift zeros(1,nnew)];
projnew.datasets=[proj.datasets cell(1,nnew)];
projnew.xcdp=[proj.xcdp cell(1,nnew)];
projnew.ycdp=[proj.ycdp cell(1,nnew)];
projnew.dx=[proj.dx ones(1,nnew)];
projnew.dy=[proj.dy ones(1,nnew)];
projnew.depth=[proj.depth zeros(1,nnew)];
projnew.texthdr=[proj.texthdr cell(1,nnew)];
projnew.texthdrfmt=[proj.texthdrfmt cell(1,nnew)];
projnew.segfmt=[proj.segfmt cell(1,nnew)];
projnew.byteorder=[proj.byteorder cell(1,nnew)];
projnew.binhdr=[proj.binhdr cell(1,nnew)];
projnew.exthdr=[proj.exthdr cell(1,nnew)];
projnew.tracehdr=[proj.tracehdr cell(1,nnew)];
projnew.bindef=[proj.bindef cell(1,nnew)];
projnew.trcdef=[proj.trcdef cell(1,nnew)];
projnew.segyrev=[proj.segyrev zeros(1,nnew)];
projnew.kxline=[proj.kxline cell(1,nnew)];
projnew.gui=[proj.gui cell(1,nnew)];
projnew.rspath=proj.rspath;
projnew.wspath=proj.wspath;
projnew.rmpath=proj.rmpath;
projnew.wmpath=proj.wmpath;
projnew.pifigures=[proj.pifigures cell(1,nnew)];
projnew.isdeleted=[proj.isdeleted zeros(1,nnew)];
projnew.deletedondisk=[proj.deletedondisk zeros(1,nnew)];
projnew.saveneeded=[proj.saveneeded zeros(1,nnew)];


end

function hpan=newdatapanel(dname,memflag,dispflag)
hsane=findsanefig;
hmpan=findobj(hsane,'tag','master_panel');
udat=get(hmpan,'userdata');
geom=udat{2};
hpanels=udat{1};
npanels=length(hpanels)+1;
panelwidth=geom(1);panelheight=geom(2);xnow=geom(3);ynow=geom(4);
wid=geom(5);ht=geom(6);xsep=geom(7);ysep=geom(8);
ynow=ynow-panelheight-ysep;
hpan=uipanel(hsane,'tag',['data_panel' int2str(npanels)],'units','normalized',...
    'position',[xnow ynow panelwidth panelheight],'userdata',npanels);
geom(4)=ynow;
hpanels{npanels}=hpan;
set(hmpan,'userdata',{hpanels geom []});
%dataset name
xn=0;yn=.1;
dg=.7*ones(1,3);
ht2=1.1;
uicontrol(hpan,'style','edit','string',dname,'tag','dataname','units','normalized',...
    'position',[xn yn wid ht],'callback','sane(''datanamechange'');','horizontalalignment','center');
uicontrol(hpan,'style','text','string','','units','normalized','position',...
    [xn+wid+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
%info button
xn=xn+wid+xsep;
wid2=(1-wid-3*xsep)/4.5;
uicontrol(hpan,'style','pushbutton','string','Information','tag','dataname','units','normalized',...
    'position',[xn yn .75*wid2 ht],'callback','sane(''datainfo'');');
uicontrol(hpan,'style','text','string','','units','normalized','position',...
    [xn+.75*wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
%memory
xn=xn+.75*wid2+xsep;
hbg1=uibuttongroup(hpan,'tag','memory','units','normalized','position',[xn yn wid2 ht]);
if(memflag==1)
    val1=1;val2=0;
else
    val1=0;val2=1;
end
uicontrol(hbg1,'style','radio','string','Y','units','normalized','position',[.1 .2 .35 .8],'value',val1,...
    'callback','sane(''datamemory'');','tag','memoryyes');
uicontrol(hbg1,'style','radio','string','N','units','normalized','position',[.6 .2 .35 .8],'value',val2,...
    'callback','sane(''datamemory'');','tag','memoryno');
uicontrol(hpan,'style','text','string','','units','normalized','position',...
    [xn+wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
%display
xn=xn+wid2+xsep;
hbg2=uibuttongroup(hpan,'tag','display','units','normalized','position',[xn yn wid2 ht]);
if(dispflag==1)
    val1=1;val2=0;
else
    val1=0;val2=1;
end
uicontrol(hbg2,'style','radio','string','Y','units','normalized','position',[.1 .2 .35 .8],'value',val1,...
    'callback','sane(''datadisplay'');','tag','displayyes');
uicontrol(hbg2,'style','radio','string','N','units','normalized','position',[.6 .2 .35 .8],'value',val2,...
    'callback','sane(''datadisplay'');','tag','displayno');

%delete button
uicontrol(hpan,'style','text','string','','units','normalized','position',...
    [xn+wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
xn=xn+wid2+xsep;
uicontrol(hpan,'style','pushbutton','string','Delete','tag','dataname','units','normalized',...
    'position',[xn yn .5*wid2 ht],'callback','sane(''datadelete'');');
uicontrol(hpan,'style','text','string','','units','normalized','position',...
    [xn+.5*wid2+.25*xsep yn .5*xsep ht2],'backgroundcolor',dg);
xn=xn+.65*wid2+xsep;
uicontrol(hpan,'style','radio','string','','units','normalized','position',[xn yn .5*wid2 ht],'value',0,...
    'callback','sane(''group'');','tag','group');

end

function waitsignalon
pos=get(gcf,'position');
spinnersize=[40 40];
waitspinner('start',[pos(3)-spinnersize(1), pos(4)-spinnersize(2), spinnersize]);
drawnow
end

function waitsignaloff
waitspinner('stop');
end

function hsane=findsanefig
hfigs=figs;
if(isempty(hfigs))
    hsane=[];
    return;
end
hfig=gcf;
if(strcmp(get(hfig,'tag'),'sane'))
    hsane=hfig;
    return;
elseif(strcmp(get(hfig,'tag'),'fromsane'))
    udat=get(hfig,'userdata');
    if(isgraphics(udat{2}))
        if(strcmp(get(udat{2},'tag'),'sane'))
            hsane=udat{2};
            return;
        end
    end
else
    
    isane=zeros(size(figs));
    for k=1:length(figs)
        if(strcmp(get(hfigs(k),'tag'),'sane'))
            isane(k)=1;
        end
    end
    if(sum(isane)==1)
        ind= isane==1;
        hsane=hfigs(ind);
    elseif(sum(isane)==0)
        hsane=[];
    else
        error('unable to resolve SANE figure')
    end
end
end

function loadprojectdialog
hsane=findsanefig;
hfile=findobj(hsane,'tag','file');
proj=get(hfile,'userdata');
dnames=proj.datanames;

isloaded=proj.isloaded;
isdisplayed=proj.isdisplayed;
islive=~proj.isdeleted;
nnames=length(dnames(islive));
pos=get(hsane,'position');
hdial=figure;
htfig=1.2*pos(4);
y0=pos(2)+pos(4)-htfig;
set(hdial,'position',[pos(1) y0 pos(3)*.6 htfig],'menubar','none','toolbar','none','numbertitle','off',...
    'name','SANE: Project load dialogue','closerequestfcn','sane(''cancelprojectload'');');
wid1=.5;wid2=.1;sep=.02;ht=.05;
xnow=.1;ynow=.95;
uicontrol(hdial,'style','text','units','normalized','position',[xnow,ynow,.8,ht],...
    'string','Choose which datasets to load and display. Initial settings are from last save.');
ynow=ynow-ht-sep;
uicontrol(hdial,'string','Dataset','units','normalized','position',[xnow,ynow,wid1,ht]);
hload=uicontrol(hdial,'string','Loaded','units','normalized','position',[xnow+wid1+sep, ynow,wid2, ht],...
    'tag','loaded');
hdisp=uicontrol(hdial,'string','Displayed','units','normalized',...
    'position',[xnow+wid1+2*sep+wid2, ynow,wid2, ht],'tag','display');
hpan=uipanel(hdial,'position',[xnow,.2,.8,.65]);
xn=.02;yn=.9;h=1.5*ht;
wd1=wid1;sp=.000;wd2=wid2;
xf=1.25;
hloaded=zeros(1,nnames);
hdisplayed=hloaded;
for k=1:length(dnames)
    if(~proj.isdeleted(k))
        mb=round(length(proj.xcoord{k})*length(proj.ycoord{k})*length(proj.tcoord{k})*4/10^6);
        uicontrol(hpan,'style','text','string',[dnames{k} ' (' int2str(mb) 'MB)'],'units','normalized','position',[xn yn wd1 h]);
        hloaded(k)=uicontrol(hpan,'style','popupmenu','string','No|Yes','units','normalized',...
            'position',[xn+wd1*xf+sp yn wd2 h],'value',isloaded(k)+1);
        hdisplayed(k)=uicontrol(hpan,'style','popupmenu','string','No|Yes','units','normalized',...
            'position',[xn+wd1*xf+wd2*xf+2*sp*xf yn wd2 h],'value',isdisplayed(k)+1);
        yn=yn-h-sp;
    end
end
ynow=.15;wid=.1;sep=.05;
uicontrol(hdial,'style','pushbutton','string','OK, continue','units','normalized',...
    'position',[xnow,ynow,2*wid,ht],'tooltipstring','Push to continue loading',...
    'callback','sane(''loadprojdial'')','tag','continue','backgroundcolor','b','foregroundcolor','w');
uicontrol(hdial,'style','pushbutton','string','All Yes','units','normalized',...
    'position',[xnow+2*wid+sep,ynow,wid,ht],'tooltipstring','set all responses to Yes',...
    'callback','sane(''loadprojdial'')','tag','allyes');
uicontrol(hdial,'style','pushbutton','string','All No','units','normalized',...
    'position',[xnow+3*wid+2*sep,ynow,wid,ht],'tooltipstring','set all responses to No',...
    'callback','sane(''loadprojdial'')','tag','allno');

set(hload,'userdata',hloaded);
set(hdisp,'userdata',hdisplayed);
set(hdial,'userdata',hsane);

end


function multiplesegyload(newproject)
hsane=findsanefig;
hfile=findobj(hsane,'tag','file');
proj=get(hfile,'userdata');
hreadmany=findobj(hsane,'tag','readmanysegy');
pos=get(hsane,'position');
hdial=figure;
htfig=1.2*pos(4);
widfig=pos(3);
y0=pos(2)+pos(4)-htfig;
set(hdial,'position',[pos(1) y0 widfig htfig],'menubar','none','toolbar','none','numbertitle','off',...
    'name','SANE: Multiple SEGY load dialogue','closerequestfcn','sane(''cancelmultipleload'');');
wid1=.8;wid2=.1;sep=.02;ht=.05;
xnow=.1;ynow=.95;
if(newproject)
    msg='Select the SEGY Datasets to be read. Then define the project save file.';
else
    msg='Select the SEGY Datasets to be read. Datasets will be included in existing project.';
end
uicontrol(hdial,'style','text','units','normalized','position',[xnow,ynow,wid1,ht],...
    'string',msg);
ynow=ynow-.5*ht-sep;
uicontrol(hdial,'style','text','units','normalized','position',[xnow,ynow,wid1,ht],...
    'string','','tag','dialmsg');
ynow=ynow-.5*ht-sep;
uicontrol(hdial,'style','pushbutton','string','New dataset','units','normalized','tag','new',...
    'position',[xnow,ynow,wid2,ht],'tooltipstring','push this to select a new dataset',...
    'callback','sane(''selectnewdataset'')');
if(newproject)
    uicontrol(hdial,'style','pushbutton','string','Project save file','units','normalized','tag','proj',...
        'position',[xnow+wid2+sep,ynow,wid2,ht],'tooltipstring','push this to define project save file',...
        'callback','sane(''defineprojectsavefile'')');
end
ynow=ynow-ht-sep;
uicontrol(hdial,'style','text','string','Project will be saved in:','units','normalized',...
    'position',[xnow,ynow,2*wid2,ht]);
if(isempty(proj.projfilename))
    projsavefile='Undefined';
else
    projsavefile=[proj.projpath proj.projfilename];
end
uicontrol(hdial,'style','text','string',projsavefile,'units','normalized','tag','projsavefile',...
    'position',[xnow+2*wid2,ynow,wid1-2*wid2,ht],'horizontalalignment','left');
ynow=ynow-.5*ht-sep;
uicontrol(hdial,'style','pushbutton','string','Done','units','normalized','tag','done',...
    'position',[xnow,ynow,wid2,ht],'tooltipstring','Push to begin reading datasets',...
    'callback','sane(''readmanysegy2'')');
uicontrol(hdial,'style','pushbutton','string','Cancel','units','normalized','tag','done',...
    'position',[xnow+wid2+sep,ynow,wid2,ht],'tooltipstring','Push to cancel reading datasets',...
    'callback','sane(''cancelmultipleload'')');

x0tab=widfig*.05;
y0tab=htfig*.05;
widtab=widfig*.9;
httab=htfig*.7;
col1=.38*widtab;
col2=.38*widtab;
col3=.1*widtab;
col4=.105*widtab;
ht=uitable(hdial,'position',[x0tab,y0tab,widtab,httab],...
    'columnname',{'Filename ','Dataset name','Display','Time shift'},...
    'columnwidth',{col1,col2,col3,col4},'columneditable',[false, true, true, true],'tag','table');
set(ht,'userdata',hreadmany);
set(hdial,'tag','fromsane','userdata',{[] , hsane});
end

function parmset=parmsetfilter(parmset)
% With no input: returns a filter parmset with default parameters
% With input parmset: checks the parmset for validity. Returns it if valid, if invalid returns error
%               message (string)
%
% The parmset is a cell array with a defined structure. It has length 3*nparms+1 where nparms is the
% number of parameters that must be defined for the task and all entries are strings. The first
% entry of the parmset is a string giving the name of the task, for example, 'spikingdecon' or
% 'domfreq' or 'filter'. Then for each of nparms parameters, there are three consequtive values: the
% name of the parameter (a string), the current parameter value (a string), and the tooltipstring.
% The current parameter value can either be a number in a string if the parameter is numeric or a
% string with choices such as 'yes|no'. The tooltipstring should be a hit for the user about the
% parameter and is displayed when the mpouse hovers over the GUI component for that parameter. See
% internal function sanetask (below) for more detail.

if(nargin<1)
    nparms=5;
    parmset=cell(1,3*nparms+1);
    parmset{1}='filter';
    parmset{2}='fmin';
    parmset{3}='10';
    parmset{4}='Define low-cut frequency in Hz';
    parmset{5}='dfmin';
    parmset{6}='5';
    parmset{7}='Define low-end rolloff in Hz (use .5*fmin if uncertain)';
    parmset{8}='fmax';
    parmset{9}='100';
    parmset{10}='Define high-cut frequency in Hz';
    parmset{11}='dfmax';
    parmset{12}='10';
    parmset{13}='Define high-end rolloff in Hz (use 10 or 20 if uncertain)';
    parmset{14}='phase';
    parmset{15}={'zero' 'minimum' 1};
    parmset{16}='Choose zero or minimum phase';
else
   fmin=str2double(parmset{3});
   dfmin=str2double(parmset{6});
   fmax=str2double(parmset{9});
   dfmax=str2double(parmset{12});
   msg=[];
   if(isnan(fmin) || isnan(fmax) || isnan(dfmin) || isnan(dfmax))
       msg='Parameters must be numbers';
   else
       if(fmin<0)
           msg=[msg '; fmin cannot be negative'];
       end
       if(dfmin<0)
           msg=[msg '; dfmin cannot be negative'];
       end
       if(dfmin>fmin)
           msg=[msg '; dfmin cannot be greater than fmin'];
       end
       if(fmax<0)
           msg=[msg '; fmax cannot be negative'];
       end
       if(dfmax>250-fmax)%bad practice but I've hardwired the Nyquist for 2 mils. I'm a baaaad boy
           msg=[msg '; dfmax is too large'];
       end
       if(fmax<fmin)
           if(fmax~=0)
               msg=[msg '; fmax must be greater than fmin'];
           end
       end
   end
   if(~isempty(msg))
        if(msg(1)==';')
           msg(1)='';
       end
       parmset=msg;
   end
end
end

function parmset=parmsetdecon(parmset,t)
% With no input: returns a filter parmset with default parameters
% With input parmset: checks the parmset for validity. Returns it if valid, if invalid returns error
%               message (string)
%
% The parmset is a cell array with a defined structure. It has length 3*nparms+1 where nparms is the
% number of parameters that must be defined for the task and all entries are strings. The first
% entry of the parmset is a string giving the name of the task, for example, 'spikingdecon' or
% 'domfreq' or 'filter'. Then for each of nparms parameters, there are three consequtive values: the
% name of the parameter (a string), the current parameter value (a string), and the tooltipstring.
% The current parameter value can either be a number in a string if the parameter is numeric or a
% string with choices such as 'yes|no'. The tooltipstring should be a hit for the user about the
% parameter and is displayed when the mpouse hovers over the GUI component for that parameter. See
% internal function sanetask (below) for more detail.

if(nargin<1)
    hsane=findsanefig;
    hfile=findobj(hsane,'tag','file');
    proj=get(hfile,'userdata');
    t=proj.tcoord{1};
    nparms=9;
    parmset=cell(1,3*nparms+1);
    parmset{1}='decon';
    parmset{2}='oplen';
    parmset{3}='0.1';
    parmset{4}='Decon operator length in seconds';
    parmset{5}='stab';
    parmset{6}='0.001';
    parmset{7}='Stability constant, between 0 and 1';
    parmset{8}='topgate';
    parmset{9}=time2str(t(1)+.25*(t(end)-t(1)));
    parmset{10}='top of design gate (seconds)';
    parmset{11}='botgate';
    parmset{12}=time2str(t(1)+.75*(t(end)-t(1)));
    parmset{13}='bottom of design gate (seconds)';
    parmset{14}='fmin';
    parmset{15}='5';
    parmset{16}='Post-decon filter low-cut frequency in Hz';
    parmset{17}='dfmin';
    parmset{18}='2.5';
    parmset{19}='Define low-end rolloff in Hz (use .5*fmin if uncertain)';
    parmset{20}='fmax';
    parmset{21}='100';
    parmset{22}='Post-decon filter high-cut frequency in Hz';
    parmset{23}='dfmax';
    parmset{24}='10';
    parmset{25}='Define high-end rolloff in Hz (use 10 or 20 if uncertain)';
    parmset{26}='phase';
    parmset{27}={'zero' 'minimum' 1};
    parmset{28}='Post-decon filter: Choose zero or minimum phase';
else
   oplen=str2double(parmset{3});
   stab=str2double(parmset{6});
   topgate=str2double(parmset{9});
   botgate=str2double(parmset{12});
   fmin=str2double(parmset{15});
   dfmin=str2double(parmset{18});
   fmax=str2double(parmset{21});
   dfmax=str2double(parmset{24});
   msg=[];
   if(isnan(fmin) || isnan(fmax) || isnan(dfmin) || isnan(dfmax) || isnan(oplen) ||isnan(stab)...
           || isnan(topgate) || isnan(botgate))
       msg='Parameters must be numbers';
   else
       if(oplen<0)
           msg=[msg '; oplen cannot be negative'];
       end
       if(oplen>1)
           msg=[msg '; oplen must be less than 1 (second)'];
       end
       if(stab<0)
           msg=[msg '; stab cannot be negative'];
       end
       if(stab>1)
           msg=[msg '; stab must be less than 1'];
       end
       if(topgate<t(1))
           msg=[msg ['; topgate must be greater than ' time2str(t(1)) 's']];
       end
       if(topgate>t(end))
           msg=[msg ['; topgate must be less than ' time2str(t(end)) 's']];
       end
       if(botgate<topgate)
           msg=[msg '; botgate must be greater than topgate' ];
       end
       if(botgate>t(end))
           msg=[msg ['; botgate must be less than ' time2str(t(end)) 's']];
       end
       if(fmin<0)
           msg=[msg '; fmin cannot be negative'];
       end
       if(dfmin<0)
           msg=[msg '; dfmin cannot be negative'];
       end
       if(dfmin>fmin)
           msg=[msg '; dfmin cannot be greater than fmin'];
       end
       if(fmax<0)
           msg=[msg '; fmax cannot be negative'];
       end
       if(dfmax>250-fmax)%bad practice but I've hardwired the Nyqiost for 2 mils
           msg=[msg '; dfmax is too large'];
       end
       if(fmax<fmin)
           if(fmax~=0)
               msg=[msg '; fmax must be greater than fmin'];
           end
       end
   end
   if(~isempty(msg))
       if(msg(1)==';')
           msg(1)='';
       end
       parmset=msg;
   end
end
end

function parmset=parmsetwavenumber(parmset)
% With no input: returns a filter parmset with default parameters
% With input parmset: checks the parmset for validity. Returns it if valid, if invalid returns error
%               message (string)

if(nargin<1)
    nparms=2;
    parmset=cell(1,3*nparms+1);
    parmset{1}='wavenumber lowpass';
    parmset{2}='sigmax';
    parmset{3}='0.125';
    parmset{4}='Define high-cut crossline wavenumber as a fraction of Nyquist';
    parmset{5}='sigmay';
    parmset{6}='0.125';
    parmset{7}='Define high-cut inline wavenumber as a fraction of Nyquist';
else
   sigmax=str2double(parmset{3});
   sigmay=str2double(parmset{6});
   msg=[];
   if(isnan(sigmax) || isnan(sigmay) )
       msg='Parameters must be numbers';
   else
       if(sigmax<0)
           msg=[msg '; sigmax cannot be negative'];
       end
       if(sigmax>1)
           msg=[msg '; sigmax cannot be greater than 1'];
       end
       if(sigmay<0)
           msg=[msg '; sigmay cannot be negative'];
       end
       if(sigmay>1)
           msg=[msg '; sigmay cannot be greater than 1'];
       end
   end
   if(~isempty(msg))
        if(msg(1)==';')
           msg(1)='';
       end
       parmset=msg;
   end
end
end

function parmset=parmsetfdom(parmset,t)
% With no input: returns a fdom parmset with default parameters
% With input parmset: checks the parmset for validity. Returns it if valid, if invalid returns error
%               message (string)
%
% The parmset is a cell array with a defined structure. It has length 3*nparms+1 where nparms is the
% number of parameters that must be defined for the task and all entries are strings. The first
% entry of the parmset is a string giving the name of the task, for example, 'spikingdecon' or
% 'domfreq' or 'filter'. Then for each of nparms parameters, there are three consequtive values: the
% name of the parameter (a string), the current parameter value (a string), and the tooltipstring.
% The current parameter value can either be a number in a string if the parameter is numeric or a
% string with choices such as 'yes|no'. The tooltipstring should be a hit for the user about the
% parameter and is displayed when the mpouse hovers over the GUI component for that parameter. See
% internal function sanetask (below) for more detail.

if(nargin<1)
    nparms=4;
    parmset=cell(1,3*nparms+1);
    parmset{1}='Dominant Frequency';
    parmset{2}='twin';
    parmset{3}='0.01';
    parmset{4}='Gaussian window half-width in seconds';
    parmset{5}='ninc';
    parmset{6}='2';
    parmset{7}='Increment between adjacent windows as an integer times the sample rate';
    parmset{8}='Fmax';
    parmset{9}='100';
    parmset{10}='Define high-cut frequency in Hz';
    parmset{11}='tfmax';
    parmset{12}='';
    parmset{13}='Time (seconds) at which Fmax occurs, leave blank if time invariant';
else
   twin=str2double(parmset{3});
   ninc=str2double(parmset{6});
   fmax=str2double(parmset{9});
   tfmax=str2double(parmset{12});
   msg=[];
   if(isnan(twin) || isnan(ninc) || isnan(fmax) )
       msg='Parameters must be numbers';
   else
       trange=t(end)-t(1);
       dt=t(2)-t(1);
       tinc=ninc*dt;
       if(twin<0)
           msg=[msg '; twin cannot be negative'];
       end
       if(twin>.1*trange)
           msg=[msg '; twin is too large'];
       end
       if(tinc<0)
           msg=[msg '; tinc=ninc*dt cannot be negative'];
       end
       if(tinc>twin)
           msg=[msg '; tinc=ninc*dt cannot be greater than twin'];
       end
       if(fmax<0)
           msg=[msg '; fmax cannot be negative'];
       end
       fnyq=.5/dt;
       if(fmax>fnyq)
           msg=[msg '; fmax is too large'];
       end
       if(isnan(tfmax))
           tfmax=[];
       end
       if(tfmax<0)
           msg=[msg '; tfmax cannot be negative'];
       end
       if(tfmax>t(end))
           msg=[msg '; tfmax too large'];
       end
   end
   if(~isempty(msg))
        if(msg(1)==';')
           msg(1)='';
       end
       parmset=msg;
   end
end
end

function val=contains(str1,str2)
ind=strfind(str1,str2);
if(isempty(ind))  %#ok<STREMP>
    val=false;
else
    val=true;
end

end

function sanetask(datasets,parmset,task)
% 
% This function is called by SANE to initiate a data-processing task on one of its datasets. It puts
% up a dialog window in which the dataset is chosen and the parameters are specified. The first two
% inputs are the list of possible datasets and the parameter set, or parmset. The list of possible
% datasets is just a cell array of strings. The parmset is also a cell array but with a defined
% structure. It has length 3*nparms+1 where nparms is the number of parameters that must be defined
% for the task and all entries are strings. The first entry of the parmset is a string giving the
% name of the task, for example, 'spikingdecon' or 'domfreq' or 'filter'. Then for each of nparms
% parameters, there are three consequtive values: the name of the parameter (a string), the current
% parameter value (a string), and the tooltipstring. The current parameter value can either be a
% number in a string if the parameter is numeric or a string with choices such as 'yes|no'.
% 
% datasets ... cell array of dataset names
% parmset ... parameter set for the task
% figsize ... length 2 vector specifying the width and height of the dialog window as a fraction of
%        the size of the SANE window.
% 

hsane=gcf;
hmsg=findobj(hsane,'tag','message');
taskname=parmset{1};
htask=figure('toolbar','none','menubar','none');
set(hmsg,'string',['Please complete parameter specifications in dialog window for task ' taskname])
% set(htask,'name',['SANE: Computation task ' task]);
pos=get(hsane,'position');
% wid=figsize(1)*pos(3);
% ht=figsize(2)*pos(4);
nparms=(length(parmset)-1)/3;
fwid=500;%width in pixels of figure
ht=30;%height of a single item in pixels
ysep=5;%y separation in pixels
xsep=5;
fht=(nparms+5)*(ht+ysep);%fig ht in pixels
%make upper left corners same
xul=pos(1);
yul=pos(2)+pos(4);
yll=yul-fht;
set(htask,'position',[xul yll fwid fht],'name',['SANE: ' taskname ' dialog'],'closerequestfcn','sane(''canceltask'')');
xnot=.1*fwid;ynot=fht-2*ht;
xnow=xnot;ynow=ynot;
wid=fwid*.8;
uicontrol(htask,'style','text','string',['Computation task: ' taskname],'units','pixels',...
    'position',[xnow ynow wid ht],'tag','task','userdata',{task parmset},'fontsize',12);
ynow=ynow-ht-ysep;
wid=fwid*.2;
uicontrol(htask,'style','text','string','Input dataset>>','units','pixels',...
    'position',[xnow ynow wid ht]);
xnow=xnow+wid+xsep;
wid=fwid*.6;
uicontrol(htask,'style','popupmenu','string',datasets,'units','pixels','tag','datasets',...
    'position',[xnow ynow wid ht],'tooltipstring','Choose the input dataset');

ynow=ynow-ht-ysep;
wid=fwid*.2;
xnow=xnot;
uicontrol(htask,'style','text','string','Output dataset>>','units','pixels',...
    'position',[xnow ynow wid ht]);
xnow=xnow+wid+xsep;
wid=fwid*.3;
outopts={'Save SEGY','Save SEGY and display','Replace input in project','Save in project as new'};
uicontrol(htask,'style','popupmenu','string',outopts,'units','pixels','tag','outputs',...
    'position',[xnow ynow wid ht],'tooltipstring','Choose the output option','value',4);

for k=1:nparms
    xnow=xnot;
    ynow=ynow-ht-ysep;
    wid=.2*fwid;
    uicontrol(htask,'style','text','string',parmset{3*(k-1)+2},'units','pixels','position',...
        [xnow,ynow,wid,ht],'fontsize',12);
    xnow=xnow+wid+xsep;
    wid=.5*fwid;
    parm=parmset{3*(k-1)+3};
    if(~iscell(parm))
        uicontrol(htask,'style','edit','string',parm,'units','pixels','position',...
            [xnow,ynow,wid,ht],'tooltipstring',parmset{3*(k-1)+4},'tag',parmset{3*(k-1)+2},'fontsize',12);
    else
        uicontrol(htask,'style','popupmenu','string',parm(1:end-1),'units','pixels','position',...
            [xnow,ynow,wid,ht],'tooltipstring',parmset{3*(k-1)+4},'tag',parmset{3*(k-1)+2},...
            'value',parm{end},'fontsize',12);
    end
end

%done and cancel buttons

xnow=.25*fwid;
wid=.3*fwid;
ynow=ynow-ht-ysep;
uicontrol(htask,'style','pushbutton','string','Done','units','pixels','tag','done',...
    'position',[xnow,ynow,wid,ht],'tooltipstring','Click to initiate the Task',...
    'callback','sane(''dotask'')');

xnow=xnow+wid+xsep;
uicontrol(htask,'style','pushbutton','string','Cancel','units','pixels','tag','cancel',...
    'position',[xnow,ynow,wid,ht],'tooltipstring','Click to cancel the Task',...
    'callback','sane(''canceltask'')');
end

function parmset=getparmset(task)
hsane=findsanefig;
hfile=findobj(hsane,'tag','file');
proj=hfile.UserData;
if(isfield(proj,'parmsets'))
    parmsets=proj.parmsets;
else
    parmsets=[];
end
%parmset=[];
for k=1:length(parmsets)
   thisparmset=parmsets{k};
   if(strcmp(thisparmset{1},task))
       parmset=thisparmset;
       return;
   end
end
%if we reach here then there is no stored parmset so we get the default one
switch task
    case 'filter'
        parmset=parmsetfilter;
    case 'spikingdecon'
        parmset=parmsetdecon;
    case 'phasemap'
        parmset=parmsetphasemap;
    case 'fdom'
        parmset=parmsetfdom;
    case 'wavenumber'
        parmset=parmsetwavenumber;
end
return

end

function setparmset(parmset)
hsane=findsanefig;
hfile=findobj(hsane,'tag','file');
proj=hfile.UserData;
if(isfield(proj,'parmsets'))
    parmsets=proj.parmsets;
    task=parmset{1};
    done=false;
    for k=1:length(parmsets)
        thisparmset=parmsets{k};
        if(strcmp(thisparmset{1},task))
            parmsets{k}=parmset;
            done=true;
            break;
        end
    end
    if(~done)
        parmsets{length(parmsets)+1}=parmset;
    end
    proj.parmsets=parmsets;
else
    proj.parmsets={parmset};
end
set(hfile,'userdata',proj);
end

function val=getparm(parmset,parm)
    nparms=(length(parmset)-1)/3;
    for k=1:nparms
        if(strcmp(parm,parmset{3*(k-1)+2}))
            parmdat=parmset{3*(k-1)+3};
            if(~iscell(parmdat))
               val=str2double(parmdat);
            else
               val=parmdat{parmdat{end}}; 
            end
        end
    end
end

function memorybuttonoff(idat)
%set the memory button off on the datapanel for dataset #idat
hsane=findsanefig;
hmpan=findobj(hsane,'tag','master_panel');
udat=get(hmpan,'userdata');
hdatapans=udat{1};
hpan=hdatapans{idat};
hno=findobj(hpan,'tag','memoryno');
hyes=findobj(hpan,'tag','memoryyes');
set(hno,'value',1)
set(hyes,'value',0)
end

function memorybuttonon(idat)
%set the memory button off on the datapanel for dataset #idat
hsane=findsanefig;
hmpan=findobj(hsane,'tag','master_panel');
udat=get(hmpan,'userdata');
hdatapans=udat{1};
hpan=hdatapans{idat};
hno=findobj(hpan,'tag','memoryno');
hyes=findobj(hpan,'tag','memoryyes');
set(hno,'value',0)
set(hyes,'value',1)
end
