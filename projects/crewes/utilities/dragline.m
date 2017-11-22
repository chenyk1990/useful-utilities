function dragline(action)
% dragline ... allows constrained dragging of a line with the mouse
%
% dragline('click')
%
% dragline is designed to allow the user to click on a line object in a
% figure window and drag it around. Behaviour is controlled by globals. The
% motion can be 'free' in both x and y or constrained to x or y. Limits on
% the motion can also be prescribed. To use this, set the 'buttondownfcn'
% of the line to call your code. Then in your callback, define the globals
% to get the desired behaviour and then your callback should call
% dragline('click'). When the drag is finished, you can have your code
% called again by providing the callback global.
%
% Globals: 
%         DRAGLINE_MOTION: either 'xonly','yonly' or 'free'
%         DRAGLINE_XLIMS: [xmin xmax] being the min and max x coordinates
%                           allowed in the drag
%         DRAGLINE_YLIMS: [ymin ymax] being the min and max y coordinates
%                           allowed in the drag
%         DRAGLINE_SHOWPOSN: set to 'on' to have a graphical display of the
%                           mouse position
%         DRAGLINE_CALLBACK: callback to execute at the conclusion of the
%                           drag. This should be a matlab command as a string as it is passed
%                           to eval.
%         DRAGLINE_MOTIONCALLBACK: callback to execute while the line is moving.  This should
%                           be a matlab command as a string as it is passed to eval.
%         DRAGLINE_PAIRED: handle(s) of related lines that are to be
%                       dragged along with the clicked line. This behaviour
%                       is invoked by a right-click only. If the actual
%                       line clicked is included in this array, then it is
%                       ignored (i.e. not shifted a second time).
%
% NOTE: This function uses the userdata of the line so you should not.
% NOTE2: Good practice is to set any globals that you are not using to null so that there is no
%       conflict with other programs that may be using dragline. For example, suppose that
%       program A uses both DRAGLINE_CALLBACK and DRAGLINE_MOTIONCALLBACK and sets them both
%       while program B uses only DRAGLINE_CALLBACK. If program B neglects to set
%       DRAGLINE_MOTIONCALLBACK to null (i.e. to '') then, if program A has been run first,
%       DRAGLINE_MOTIONCALLBACK may still be set to call program A when program B is used.
%       Program B should set this global to null to avoid this conflict.
%
% For an example of using this, run the function test_dragline with no
% inputs. Edit test_dragline to see the code.
%
% G.F. Margrave, Devon, 2016
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

%

% if(nargin<1) action='init'; end
% if(strcmp(action,'init'))
%   set(gcf,'windowbuttondownfcn','dragline(''click'')');
%   return;
% end

global DRAGLINE_MOTION DRAGLINE_XLIMS DRAGLINE_YLIMS DRAGLINE_SHOWPOSN DRAGLINE_POSNTXT DRAGLINE_CALLBACK DRAGLINE_PAIRED DRAGLINE_MOTIONCALLBACK
% DRAGLINE_POSNTXT is used internally and should not be set by user.


if(strcmp(action,'click'))
    hline=gco;
%     obj=get(gco,'type');
%     if(~strcmp(obj,'line'))
%        msgbox('Click on a line!!!')
%        return;
%     end

    pt=get(gca,'currentpoint');
    set(hline,'userdata',pt(1,1:2));
    %set(hline,'erasemode','xor','linestyle','.');

    set(gcf,'windowbuttonmotionfcn','dragline(''move'')');
    set(gcf,'windowbuttonupfcn','dragline(''fini'')');
    return;
end
if(strcmp(action,'move'))
    hline=gco;
    motion=DRAGLINE_MOTION;
    if(isempty(motion))
        motion='free';
    end
    xlims=DRAGLINE_XLIMS;
    if(isempty(xlims))
        xlims=get(gca,'xlim');
    end
    ylims=DRAGLINE_YLIMS;
    if(isempty(ylims))
        ylims=get(gca,'ylim');
    end
    showposn=DRAGLINE_SHOWPOSN;
    if(isempty(showposn))
        showposn='off';
    end
    ptxt=DRAGLINE_POSNTXT;
    mcb=DRAGLINE_MOTIONCALLBACK;
    pt1=get(hline,'userdata');
    pt2=get(gca,'currentpoint');
    pt2=pt2(1,1:2);
    %check bounds
    if(pt2(1)<xlims(1))
        pt2(1)=xlims(1);
    end
    if(pt2(1)>xlims(2))
        pt2(1)=xlims(2);
    end
    if(pt2(2)<ylims(1))
        pt2(2)=ylims(1);
    end
    if(pt2(2)>ylims(2))
        pt2(2)=ylims(2);
    end
    %compute displacements
    delx=pt2(1)-pt1(1);
    dely=pt2(2)-pt1(2);
    %contrained motion
    if(strcmp(motion,'free'))
        x=get(hline,'xdata');
        y=get(hline,'ydata');
        set(hline,'xdata',x+delx,'ydata',y+dely,'userdata',pt2);
        pstring=['(' num2str(pt2(1)) ',' num2str(pt2(2)) ')'];
        xshift=delx;yshift=dely;
    elseif(strcmp(motion,'xonly'))
        x=get(hline,'xdata');
        pt2=[pt2(1) pt1(2)];
        set(hline,'xdata',x+delx,'userdata',pt2);
        pstring=['(' num2str(pt2(1)) ')'];
        xshift=delx;yshift=0;
    elseif(strcmp(motion,'yonly'))
        y=get(hline,'ydata');
        pt2=[pt1(1) pt2(2)];
        set(hline,'ydata',y+dely,'userdata',pt2); 
        pstring=['(' num2str(pt2(2)) ')'];
        xshift=0;yshift=dely;
    end
    if(strcmp(showposn,'on'))
        if(isempty(ptxt))
            ptxt=text(pt2(1),pt2(2),pstring,...
                'backgroundcolor','w');
            DRAGLINE_POSNTXT=ptxt;
        else
            set(ptxt,'position',[pt2 0],'string',pstring);
        end
    end
    seltype=get(gcf,'selectiontype');
    if(strcmp(seltype,'alt'))
        nlines=length(DRAGLINE_PAIRED);
        for k=1:nlines
            if(isgraphics(DRAGLINE_PAIRED(k)))
                if(DRAGLINE_PAIRED(k)~=hline)
                    x=get(DRAGLINE_PAIRED(k),'xdata');
                    y=get(DRAGLINE_PAIRED(k),'ydata');
                    set(DRAGLINE_PAIRED(k),'xdata',x+xshift,'ydata',y+yshift);
                end
            end
        end
    end
    if(~isempty(mcb))
        eval(mcb);
    end
    return;
end

if(strcmp(action,'fini'))
    
    set(gcf,'windowbuttondownfcn','');
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    ptxt=DRAGLINE_POSNTXT;
    if(~isempty(ptxt))
        delete(ptxt);
        DRAGLINE_POSNTXT=[];
    end
    cb=DRAGLINE_CALLBACK;
    if(~isempty(cb))
        eval(cb);
    end

    return;
end	