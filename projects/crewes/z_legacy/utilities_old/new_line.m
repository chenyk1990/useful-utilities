function new_line(action,arg2,arg3)
% new_line(action,transfer)
% 
% NEW_LINE is a general tool for creating new lines with the mouse in the current
% figure. Capabilities include: 
%		1) drawing in a line with simple mouse clicks
%		2) duplicating an existing line
%		3) breaking an existing line into two independent segments
%		4) deleting an entire line with a single mouse click
%
% NEW_LINE dose not include the ability to modify the created lines in any way. For
% this purpose use the EDITLINES facility.
%
% This is not usually called directly but rather initiated through NEWLINEINIT
%
% mouse action definitions:
% 
% button 1 click ... draw in a line. The first click defines the line and subsequent
%		clicks append new points to the previous ones.
% button 2 click ... delete an entire line
% button 3 click ... duplicate an existing line
%
% NEW_LINE is intended to be incorporated within larger tools with specific 
% tasks. Someone using it it that context need read no further.
% For further information type: help newlineinit
%  
% by G.F. Margrave December 1993
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
% mouse button assignments
% button1 up ... a point is added to the end of the current line
%	causes initiation of drawing mode
% button2 up ... if not in drawing mode, then the line clicked is deleted
%	if in drawing mode, then the current line is completed as a polygon
% button3 up ... if not in drawing mode, then the line clicked is duplicated
%	if in drawing mode, the the current line is finished
%
% The program 'terminates' everytime one of the following actions occurs:
%	1) Button 3 is clicked while in draw mode. Then the current line is finished and its 
%	handle returned
%	2) Button 3 is clicked while not in draw mode. Then the handle of the newly duplicated
%	line is returned
%	3) Button 2 is clicked while not in draw mode. Then the negative of the handle of the
%	deleted line is returned
% At termination, the handle of the new line and its x and y coordinates are stored in
% the current axes userdata. Running newlinefini is the easiset way to retrieve this
%
if( strcmp(action,'init') )
	if(nargin==1)
		transfer=[];
	else
		transfer=arg2;
	end
	set(gca,'userdata',[]);
	% make an invisible storage bucket 
	uicontrol('style','text','visible','off','tag','yrag_new',...
		'userdata',[]);
	uicontrol('style','text','visible','off','tag','yrag_trans',...
		'userdata',transfer);
	% set the windowbuttondownfcn
	set(gcf,'windowbuttondownfcn','newline(''buttondown'')');
	set(gcf,'windowbuttonupfcn','');
	set(gcf,'windowbuttonmotionfcn','');
	return;
end
if(strcmp(action,'fini') )
	% clean up graphics if needed
	% Get the storage buckets
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				if( strcmp(get(h(k),'string'),'yrag_trans') )
% 					htrans=h(k);
% 					found=found+1;
% 				end
% 				
% 				if( found== 2)
% 					break;
% 				end
% 			end
% 		end
% 	end
    hstor=findobj(gcf,'tag','yrag_new');
    htrans=findobj(gcf,'tag','yrag_trans');
	returnstuff=get(hstor,'userdata');
	set(gca,'userdata',returnstuff);
	%
	% the userdata of hstor is
	% [key hline x(1:npts) y(1:npts)]
	%
	% set the windowbuttondownfcn
	set(gcf,'windowbuttondownfcn','');
	set(gcf,'windowbuttonupfcn','');
	
	% transfer control if requested
	transfer=get(htrans,'userdata');
	if(~isempty(transfer))
		eval(transfer);
	end
	% make sure cursor is an arrow
	set(gcf,'pointer','arrow');
end
if( strcmp(action,'buttondown') )
	% Get the storage bucket
    hstor=findobj(gcf,'tag','yrag_new');
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				if( found== 1)
% 					break;
% 				end
% 			end
% 		end
% 	end
	% if there is information in the userdata of hstor, then we have a
	% line currently active
	lineinfo=get(hstor,'userdata');
	new=0;
	if(isempty(lineinfo))
		new=1;
	end
	% determine which button was pressed
	button=1;
	flag=get(gcf,'selectiontype');
	if(strcmp(flag,'extend') )
		button=2;
	elseif(strcmp(flag,'alt'))
		button=3;
	end
	% assign the windowbuttonupfcn
	if(new&&button==2)
		set(gcf,'windowbuttonupfcn','new_line(''kill'')');
	elseif(~new&&button==2)
		set(gcf,'windowbuttonupfcn','new_line(''poly'')');
	elseif(new&&button==3)
		set(gcf,'windowbuttonupfcn','new_line(''duplicate'')');
	elseif(~new&&button==3)
		set(gcf,'windowbuttonupfcn','new_line(''fini'')');
	else
		set(gcf,'windowbuttonupfcn','new_line(''draw'')');
	end
	return
end
% line drawing section
if(strcmp(action,'draw'))
	set(gcf,'windowbuttonupfcn','');
	
	% Get the storage buckets
    hstor=findobj(gcf,'tag','yrag_new');
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				
% 				if( found== 1)
% 					break;
% 				end
% 			end
% 		end
% 	end
	%get the current points
	lineinfo=get(hstor,'userdata');
	% get the current point
	pt=get(gca,'currentpoint');
	xpt=pt(1,1);
	ypt=pt(1,2);
	%determine whether to extend and existing line or start a new one
	if(isempty(lineinfo))
		hline=line(xpt,ypt,'color',[.5 .5 .5]);
		
		% set the lineinfo in userdata
		set(hstor,'userdata',{1 hline xpt ypt 0});
	else
		% save the undo information
		set(gca,'userdata',lineinfo);
		%xy=lineinfo(3:length(lineinfo));
		x=lineinfo{3};
		y=lineinfo{4};
		x=[x xpt];
		y=[y ypt];
        set(lineinfo{2},'xdata',x,'ydata',y);
        lineinfo{3}=x;
        lineinfo{4}=y;
        set(hstor,'userdata',lineinfo);
	end
	
	return;
end
% terminate as a polygon
if(strcmp(action,'poly'))
	set(gcf,'windowbuttonupfcn','');
	
	% Get the storage buckets
    hstor=findobj(gcf,'tag','yrag_new');
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				
% 				if( found== 1)
% 					break;
% 				end
% 			end
% 		end
% 	end
	%get the current points
	lineinfo=get(hstor,'userdata');
	x=lineinfo{3};
	y=lineinfo{4};
	x=[x x(1)];
	y=[y y(1)];
	set(lineinfo{2},'xdata',x,'ydata',y);
    lineinfo{3}=x;
    lineinfo{4}=y;
	set(hstor,'userdata',lineinfo);
	new_line('fini');
	return;
end
if(strcmp(action,'undo') )
% Get the storage buckets
    hstor=findobj(gcf,'tag','yrag_new');
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				
% 				if( found== 1)
% 					break;
% 				end
% 			end
% 		end
% 	end
	
	if(isempty(hstor))
		if( nargin==3 )
			transfer=arg3;
		else
			transfer=[];
		end
		new_line('init',transfer);
        hstor=findobj(gcf,'tag','yrag_new');
% 		h=get(gcf,'children');
% 		found=0;
% 		for k=1:length(h)
% 			if( strcmp(get(h(k),'type'),'uicontrol') )
% 				if( strcmp(get(h(k),'style'),'text') )
% 					if( strcmp(get(h(k),'string'),'yrag_new') )
% 						hstor=h(k);
% 						found=found+1;
% 					end
% 				
% 					if( found== 1)
% 						break;
% 					end
% 				end
% 			end
% 		end
	end
	
	if( nargin< 2 )
		% get the undo information
		undostuff=get(gca,'userdata');
	else
		undostuff = arg2;
	end
	
	if(isempty(undostuff))
		return;
	end
	
	% see if we have color information. This means we are restoring a line that
	% was deleted
	%ind=findthing(undostuff,nan);
	if(isnan(undostuff{5}))
% 		colorinfo=undostuff(ind(1)+1:length(undostuff));
		kolor=undostuff{6};
		ls=undostuff{7};
% 		if(isnan(ls(2))) ls=ls(1); end
% 		ls=setstr(ls);
		lw=undostuff{8};
		%em=setstr(colorinfo(7:length(colorinfo)));
		undostuff{5}=0;
	end
	
	key=undostuff{1};
	hline=undostuff{2};
	x=undostuff{3};
	y=undostuff{4};
	
	if( key == 4)
		hline=line(x,y,'color',kolor,'linestyle',ls,'linewidth',lw);
		set(hstor,'userdata',{3 hline x y});
	elseif(key==1)
		set(hline,'xdata',x,'ydata',y);
		set(hstor,'userdata',{1 hline x y});
	end
	
	
	
	% if we are restoring a deleted line, then we must terminate the newline
	% mode so that the line handle is returned
	if(~isempty(colorinfo))
		new_line('fini');
	end
	return;
end
% line duplication
if(strcmp(action,'duplicate') )
	set(gcf,'windowbuttonupfcn','');
	
	% Get the storage buckets
    hstor=findobj(gcf,'tag','yrag_new');
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				
% 				if( found== 1)
% 					break;
% 				end
% 			end
% 		end
% 	end
% determine if the user clicked on a line
	hline=get(gcf,'currentobject');
	flag=get(hline,'type');
	if(~strcmp(flag,'line'))
		return;
	end
	x=get(hline,'xdata');
	y=get(hline,'ydata');
	ylim=get(gca,'ylim');
	yshift=abs(ylim(1)-ylim(2))*.05;
	y=y+yshift;
	hline=line(x,y,'color',[.5 .5 .5]);
	set(hstor,'userdata',{2 hline x y 0});
	set(gca,'userdata',[]);% undo for duplicate is just delete
	new_line('fini');
return;
end
% kill an entire line
if(strcmp(action,'kill'))
	
	set(gcf,'windowbuttonupfcn','');
	
	% Get the storage buckets
    hstor=findobj(gcf,'tag','yrag_new');
% 	h=get(gcf,'children');
% 	found=0;
% 	for k=1:length(h)
% 		if( strcmp(get(h(k),'type'),'uicontrol') )
% 			if( strcmp(get(h(k),'style'),'text') )
% 				if( strcmp(get(h(k),'string'),'yrag_new') )
% 					hstor=h(k);
% 					found=found+1;
% 				end
% 				
% 				if( found== 1)
% 					break;
% 				end
% 			end
% 		end
% 	end
	
% determine if the user clicked on a line
	hline=get(gcf,'currentobject');
	flag=get(hline,'type');
	if(~strcmp(flag,'line'))
		return;
	end
	x=get(hline,'xdata');
	y=get(hline,'ydata');
	% get the lines color, style, etc
	
	kolor=get(hline,'color');
	ls=abs(get(hline,'linestyle'));
	if(length(ls)==1); ls=[ls nan]; end
	lw=get(hline,'linewidth');
	%em=get(hline,'erasemode');
	% set the undo information
	undostuff={4 hline x y NaN kolor ls lw};
	set(gca,'userdata',undostuff);
	% delete the line
	delete(hline);
	% set the information in hstor
	set(hstor,'userdata',undostuff);
	% newlinefini
	new_line('fini');
	return;
end