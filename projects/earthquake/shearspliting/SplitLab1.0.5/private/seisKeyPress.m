function seisKeyPress(src,evnt,seis)
% handle keypress within Splitlab seismmogram plot

global thiseq
% local filter is temporary filter structure, filt structure in calling
% function will be assigned to the new value
f1 = thiseq.filter(1);
f2 = thiseq.filter(2);



if strcmp(evnt.Key,'shift')
    return

elseif strcmp(evnt.Key,'return')
    if ~isfield(thiseq, 'a')|| (isempty(thiseq.a) && isempty(thiseq.f)) 
        errordlg('no time window picked... Sorry, can''t split','Error')
        return
    else
        preSplit
        return
    end
%% STILL A bug in seismogram order 08.05.06    
% elseif strcmp(upper(evnt.Key),'L')
%     %lock aspect
%     seis  = findobj('Tag','seismo');
%     lockbut = findobj('Type','uitoggletool','Tag','LockButton');
%     switch get(lockbut,'State')
%         case 'on'%lock to common Y-limits
%             for k=1:3
%                 yyy(k,:)   =get(get(seis(k),'Parent'),'YLim');
%             end
%             yyy =[min(yyy(:,1)) max(yyy(:,2))];
%             for k=1:3
%                 set(get(seis(k),'Parent'), 'YLim',yyy,'YLimMode','manual')
%             end
%             set(lockbut ,'State','off','Cdata',getfield(get(lockbut ,'UserData'),'locked'))
%         case 'off'%unlock Y-limits
%             for k=1:3
%                 set(get(seis(k),'Parent'),'YLimMode','auto')
%             end
%             set(lockbut,'State','on','Cdata',getfield(get(lockbut,'UserData'),'unlocked'))
%     end
%     
%     
elseif strcmp(evnt.Key,'home')|strcmp(evnt.Key,'escape')
    %jump close to selected phase
    val  = get(findobj('Tag','PhaseSelector'),'Value');
    t_home = floor(thiseq.phase.ttimes(val)/10)*10 - 30; %~30 seconds before phase; at full 10 seconds
    xlim([t_home t_home+150]) % timewindow of 150 sec

elseif strcmp(evnt.Key,'pagedown')
    xx = xlim;
    xlim(xx-diff(xx)/5)

elseif strcmp(evnt.Key,'pageup')
    xx=xlim;
    xlim(xx+diff(xx)/5)
elseif strcmp(evnt.Key,'backspace')
    xlim('auto')
% elseif strcmp(evnt.Key,'m')
%     Aniso_multi;
else

    switch evnt.Character
        case 'p'
            polarisation2006
            return
        case 'f'
            [f1, f2] = filterdialog(thiseq.filter);
        case '0'
            f1 = 0;
            f2 = inf;
        case '1'
            f1 = 0.01;
            f2 = 0.1;
        case '2'
            f1 = 0.02;
            f2 = 0.2;
        case '3'
            f1 = 0.02;
            f2 = 0.3;
        case '4'
            f1 = 0.01;
            f2 = 0.3;
        case '5'
            f1 = 0.01;
            f2 = 0.4;
        case '6'
            f1 = 0.02;
            f2 = 1.0;
        case '7'
            f1 = 0.01;
            f2 = 0.15;
         case '8'
            f1 = 0.02;
            f2 = 0.25;
        case '9'
            f1 = 0.01;
            f2 = 1.0;

        case '+'
            f1  = f1 + 0.002;
            if f1  >= f2
                f1 = f2-0.002;
            elseif f2  <= f1
                f1 = f1 - 0.002;
            end
        case '-'
            f1  = f1 - 0.002;
            if f1  < 0
                f1 = 0;
            end
        case '*'
            f2  = f2 + 0.02;
        case '/'
            f2  = f2 - 0.02;
            if f2  < 0
                f2 = 0;
            elseif f2  <= f1
                f2 = f1 + 0.002;
            end

        case ' ' %space
            % rotate system
            button=findobj('Tag','SystemButton');
            if thiseq.system=='ENV';
                thiseq.system='LTQ';
                set(button, 'State','On')
            elseif thiseq.system=='LTQ';
                thiseq.system='ENV';
                set(button, 'State','Off')
            end
%             SL_updatefiltered(flipud(findobj('Tag','seismo')))
            return
    end %switch
end %else

thiseq.filter = [f1 f2]; %update variable filt in caller function
set(gcbf,'KeyPressFcn',{@seisKeyPress,seis})

%% filtering and display...
SL_updatefiltered(seis);
