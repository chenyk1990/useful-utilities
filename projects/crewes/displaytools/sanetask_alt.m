function sanetask(datasets,parmset,figsize,callback)
% 
% This function is called by SANE to initiate a data-processing task on one of its datasets The
% first two inputs are the list of possible datasets and the parameter set,or parmset. The list of
% possible datasets is just a cell array of strings. The parmset is also a cell array but with a
% defined structure. It has length 3*nparms+1 where nparms is the number of parameters that must be
% defined for the task. The first entry of the parmset is always a string giving the name of the
% task, for  example, 'gabordecon' or 'deconw' or 'filter'. Then for each of nparms parameters,
% there are three consequtive values: the name of the parameter (a string), the current parameter
% value, and the tooltipstring. The current parameter value can either be a number in a string if
% the parameter is numeric or a string with choices such as 'yes|no'.
% 
% datasets ... cell array of dataset names
% parmset ... parameter set for the task
% figsize ... length 2 vector specifying the width and height of the dialog window as a fraction of
%        the size of the SANE window.
% callback ... strong specifying the callback to execute when the user clicks the done button. For
%       example 'sane(''gabordecon'')' .
% 

if(iscell(datatsets))
    action='init';
else
    action='datasets';
end

if(strcmp(action,'init'))
hsane=gcf;
taskname=parmset{1};
htask=figure;
pos=get(hsane,'position');
wid=figsize(1)*pos(3);
ht=figsize(2)*pos(4);
%make upper left corners same
xul=pos(1);
yul=pos(2)+pos(4);
yll=yul-ht;
set(htask,'position',[xul yll wid ht],'name',['SANE: ' taskname ' dialog']);
xnot=.05;ynot=.9;
xnow=xnot;ynow=ynot;
ht=.05;wid=.1;
xsep=.02;ysep=.02;
uicontrol(htask,'style','text','string','Input dataset>>','units','normalized',...
    'position',[xnow ynow wid ht]);
xnow=xnow+wid+xsep;
uicontrol(htask,'style','popupmenu','string',datasets,'units','normalized','tag','datasets',...
    'position',[xnow ynow wid ht],'tooltipstring','Choose the input dataset');
nparms=(length(parmset)-1)/3;
for k=1:nparms
    xnow=xnot;
    ynow=ynow-ht-ysep;
    wid=.3;
    uicontrol(htask,'style','text','string',parmset{3*(k-1)+2},'units','normalized','position',...
        [xnow,ynow,wid,ht]);
    xnow=xnow+wid+xsep;
    wid=.4;
    parm=parmset{3*(k-1)+3};
    if(~contains(parm,'|'))
        uicontrol(htask,'style','edit','string',num2str(parm),'units','normalized','position',...
            [xnow,ynow,wid,ht],'tooltipstring',parmset{3*(k-1)+4},'tag',parmset{3*(k-1)+2});
    else
        uicontrol(htask,'style','popupmenu','string',parm,'units','normalized','position',...
            [xnow,ynow,wid,ht],'tooltipstring',parmset{3*(k-1)+4},'tag',parmset{3*(k-1)+2});
    end
end

%done and cancel buttons

xnow=xnot;
wid=.3;
ynow=ynow-ht-ysep;
uicontrol(htask,'style','pushbutton','string','Done','units','normalzed','tag','done',...
    'position',[xnow,ynow,wid,ht],'userdata',callback,'tooltipstring','Click to initiate the Task',...
    'callback','sanetask(''done'')');

xnow=xnow+wid+sep;
uicontrol(htask,'style','pushbutton','string','Cancel','units','normalzed','tag','cancel',...
    'position',[xnow,ynow,wid,ht],'tooltipstring','Click to cancel the Task',...
    'callback','sanetask(''cancel'')','userdata',hsane);

elseif(strcmp(action,'done'))
    callback=get(gcbo,'userdata');
    eval(callback);
elseif(strcmp(action,'cancel'))
    hsane=get(gcbo,'userdata');
    close(gcf);
    figure(hsane);
    sane('taskcancelled');
end


