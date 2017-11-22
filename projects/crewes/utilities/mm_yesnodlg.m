function answer=mm_yesnodlg(statemnt,question,title,default,parent)
%function answer=mm_yesnodlg(statemnt,question,title,default,parent)
%
% Where:
%   statemnt   = statement of problem (char)
%   question   = question to be asked (char)
%   title      = title to display at top of dialog (char)
%   default    = default answer, 'Yes' or 'No'
%   parent     = handle to parent figure. Empty implies the parent is the
%                Matlab desktop
%
% This function uses mm_adjust to cause a modified msgbox dialog to
% display on the same monitor as the parent figure
%
%
%  Authors: Kevin Hall, 2017
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

narginchk(3,5)

if nargin<4 || isempty(default)
    default='No';
end
if nargin<5
    parent = []; %Matlab desktop
end

if isempty(parent) || isa(parent,'handle')
    answer='No';
    
    h=msgbox([statemnt ' ' question],title,'warn','modal');
    mm_adjust(h,parent);
    
    set(h,'DeleteFcn',@close_request)
    set(h,'CloseRequestFcn',@close_request)
    set(h,'KeyPressFcn',@keypress_callback)
    
    c = get(h,'children');
    btn1 = c(1);
    
    btn2 = copy(btn1);
    
    btn1.String='Yes';
    btn1.Position(1) = btn1.Position(1)-btn1.Position(3)/2-2;
    set(btn1,'Callback',@button1_callback);
    btn2.String='No';
    btn2.Position(1) = btn2.Position(1)+btn2.Position(3)/2+2;
    set(btn2,'Callback',@button2_callback);
    
    btn2.Parent=h;
    
    switch lower(default(1))
        case 'y'
            h.setDefaultButton(btn1); %undocumented! may not work in future
            btn1.Value=1;
        case 'n'
            h.setDefaultButton(btn2); %undocumented! may not work in future
            btn1.Value=0;
    end
    
    uiwait(h); %wait for button or key press
    
else
    disp(statemnt)
    answer = input([question ' '],'s');    
end
    
%%begin nested functions
    function keypress_callback(hObject,eventdata)
        switch get(btn1,'Value')
            case 0
                answer='No';
            case 1
                answer='Yes';
        end
        close_request(hObject,eventdata)
    end
    function button1_callback(hObject,eventdata)
        answer = 'Yes';
        close_request(hObject,eventdata)
    end
    function button2_callback(hObject,eventdata)
        answer = 'No';
        close_request(hObject,eventdata)
    end
    function close_request(hObject,eventdata)
        uiresume; %required for output
        delete(gcbf); %close msgbox
    end
%end nested functions

end %end function mm_yesnodlg
