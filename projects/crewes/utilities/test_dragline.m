function test_dragline(action)
%
% test_dragline
%
% Execute this function without inputs to see an example of using dragline
% with uncontrained motion. Execute test_dragline('xonly') to see motion
% constrained in x only and similarly test_dragline('yonly'). Finally test_dragline('motioncb')
% demonstrates using the motion callback option to cause the line to "flash" while being moved.
%
% The function puts up a window with two lines plotted. One is a
% reflectivity and the other is the corresponding seismic trace. Left
% clicking on either line allows it to be dragged independently.
% Right-clicking on either line allows both lines to be dragged together.
%
global DRAGLINE_MOTION DRAGLINE_XLIMS DRAGLINE_YLIMS DRAGLINE_SHOWPOSN DRAGLINE_POSNTXT DRAGLINE_CALLBACK DRAGLINE_PAIRED DRAGLINE_MOTIONCALLBACK
if(nargin<1)
    action='init';
    motion='free';
end

if(strcmp(action,'xonly')||strcmp(action,'yonly')||strcmp(action,'free')||strcmp(action,'motioncb'))
    if(strcmp(action,'motioncb'))
        motion='free';
        mcb=1;
    else
        motion=action;
        mcb=0;
    end
    action='init';
end

if(strcmp(action,'init'))
    [r,t]=reflec(1,.001);
    [w,tw]=ricker(.001,30);
    s=convz(r,w);
    s=s*norm(r)/norm(s);
    
    figure
    hh=trplot(t,[r s],'yaxis','on');%hh contains the handles of both lines in a cell array
    title({['motion ' motion],'try left-click and drag and then right-click and drag'})
    set(hh{1},'buttondownfcn','test_dragline(''doit'')');%note the assignment to both lines
    set(hh{2},'buttondownfcn','test_dragline(''doit'')');
    DRAGLINE_PAIRED=[hh{1} hh{2}];%this sets up paired dragging with a right click
    set(gca,'userdata',{motion,mcb})%this is just putting the motion option where it can be retrieved in a callback
    ylim([-.5 .5]);
    xlim([-.5 2]);
    prepfiga
elseif(strcmp(action,'doit'))%this is called when a click occurs on one of the lines
    udat=get(gca,'userdata');
    motion=udat{1};
    mcb=udat{2};
    DRAGLINE_MOTION=motion;
    if(mcb==1)
        DRAGLINE_MOTIONCALLBACK='test_dragline(''flash'')';
    else
        DRAGLINE_MOTIONCALLBACK='';
    end
    dragline('click')
elseif(strcmp(action,'flash'))
    hline=gco;
    ls=get(hline,'linestyle');
    if(strcmp(ls,'-'))
        set(hline,'linestyle','none','marker','o');
    else
        set(hline,'linestyle','-','marker','none');
    end
        
end