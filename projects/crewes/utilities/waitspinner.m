function waitspinner(action,posn)
global jObj

if(strcmp(action,'start'))
    if(length(posn)<4)
        posn=[posn(1:2) 80 80];
    end
    
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
    javacomponent(jObj.getComponent, posn, gcf);
    jObj.start;
 
elseif(strcmp(action,'end')||strcmp(action,'stop'))
    jObj.stop;
    jObj.setBusyText('All done!');
    clear jObj
end
    
    