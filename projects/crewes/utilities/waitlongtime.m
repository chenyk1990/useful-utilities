function waitlongtime(action,msg,pct,hwait)
global jObj
if(nargin<3)
    pct=25;
end
if(nargin<2)
    msg='';
end

if(strcmp(action,'start'))
    pos=get(gcf,'position');
    xcenter=round(pos(1)+.5*pos(3));
    ycenter=round(pos(2)+.5*pos(4));
    hwait=figure;
    wid=pos(3)*pct/100;
    ht=pos(4)*pct/100;
    set(hwait,'position',[xcenter-round(.5*wid) ycenter-round(.5*ht) wid ht],'toolbar','none',...
        'menubar','none','name','Wait awhile','numbertitle','off','tag','waitlongtime');
    
    try
        % R2010a and newer
        iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
        iconsSizeEnums = javaMethod('values',iconsClassName);
        SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
        jObj = com.mathworks.widgets.BusyAffordance(SIZE_32x32, '');  % icon, label
    catch
        % R2009b and earlier
        redColor   = java.awt.Color(1,0,0);
        blackColor = java.awt.Color(0,0,0);
        jObj = com.mathworks.widgets.BusyAffordance(redColor, blackColor);
    end
    jObj.setPaintsWhenStopped(false);  % default = false
    jObj.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
    javacomponent(jObj.getComponent, [round(.5*wid)-40,round(ht/3-40),80,80], gcf);
    jObj.start;
    uicontrol(hwait,'style','text','string',msg,'units','normalized','position',[.1 .9 .8 .1],...
        'fontsize',12)
elseif(strcmp(action,'end')||strcmp(action,'stop'))
    if(nargin<4)
        hwait=gcf;
    end
    jObj.stop;
    jObj.setBusyText('All done!');
    tag=get(hwait,'tag');
    if(strcmp(tag,'waitlongtime'))
        close(hwait);
    end
    clear jObj
end
    
    