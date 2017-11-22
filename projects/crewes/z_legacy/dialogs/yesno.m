function yesno(action)
%
% This is part of a 3 piece tool. The other two parts are yesnoinit and yesnofini. See the hel
% for these two for a description.
%
% G.F. Margrave Feb 1994
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

if(strcmp(action,'init'))
		hax=gca;
		qs=get(hax,'userdata');
        pos=get(gcf,'position');
        w=pos(3);%width of calling figure
        h=pos(4);%height of calling figure
        xc=pos(1)+.5*w;
        yc=pos(2)+.5*h;%these are center of calling figure
        wd=400;
        hd=100;%these are width and height of the dialog
		%build the dialog box
		figure('visible','on','menubar','none','name','Please answer this question',...
            'numbertitle','off','position',[xc-.5*wd yc-.5*hd wd hd],'windowstyle','modal');
		%
		% now the question
		%
		q=qs{1};

		xnow=.1;ynow=.4;
        width=.8;
        height=.5;
		uicontrol('style','text','string',q,'units','normalized',...
					'position',[xnow ynow width height],'fontsize',12);

		% the yes/no buttons

		width=1/5;
        sep=1/20;
        xnow=width-sep;
		ynow=.1;
        height=.2;

		uicontrol('style','pushbutton','string','Yes','units','normalized','position',...
			[xnow,ynow,width,height],'callback','yesno(''button'')');
		
        xnow=xnow+sep+width;
		uicontrol('style','pushbutton','string','No','units','normalized','position',...
			[xnow,ynow,width,height],'callback','yesno(''button'')');

		
		% cancel button

		xnow=xnow+sep+width;
		uicontrol('style','pushbutton','string','Cancel','units','normalized','position',...
			[xnow,ynow,width,height],'callback','yesno(''button'')');


		transfer=qs{2};
		set(gcf,'userdata',{hax transfer});


		return;

	end

	%
	% handle a button
	%
	if(strcmp(action,'button'))
		h=get(gcf,'userdata');

		hbutton=gco;

		flag=get(hbutton,'string');

		if(strcmp(flag,'Yes'))
			reply=1;
		elseif(strcmp(flag,'No'))
			reply=0;
		else
			reply=-1;
		end
		
		set(h{1},'userdata',reply);

		close(gcf);

	% call the transfer expression
		transfer=h{2};

		eval(transfer);

		return;

	end