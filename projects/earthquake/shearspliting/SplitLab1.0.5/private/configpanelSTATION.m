%% Station data
h.panel(1) = uipanel('Units','pixel','Title','Station',...
    'FontSize',10,'Position',pos-[0 -10 0 10],'Visible','off', ...
    'BackgroundColor', [224   223   227]/255 );
set(cfig,'color',get(h.panel(1),'BackgroundColor'));

x = 60;
y = pos(4)-70;
w = 60;
v = 20;
%% field descriptions text
txt = {'Station Code:';
    'Network code';
    'Station latitude';
    'Station longitude'};

for i=0:3
    uicontrol('Parent',h.panel(1),'Units','pixel',...
        'Style','text',...
        'Position',[x y-(i*v*1.75) 250 v],...
        'String', txt{i+1},...
        'HorizontalAlignment','Left');
end

uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Pushbutton',...
    'Position',[pos(3)-160 y+3 120 25],...
    'String', 'IRIS station query',...
    'Tooltip','http://www.iris.edu/SeismiQuery/station.htm',...
    'Callback','web http://www.iris.edu/SeismiQuery/station.htm -browser');
uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Pushbutton',...
    'Position',[pos(3)-160 y-32 120 25],...
    'String', 'Orfeus station query',...
    'Tooltip','http://orfeus.knmi.nl/wg/wg1.htm ',...
    'Callback','web http://orfeus.knmi.nl/wg/wg1.htm -browser');
uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Pushbutton',...
    'Position',[pos(3)-160 y-67 120 25],...
    'String', 'NEIC station book',...
    'Tooltip','http://wwwneic.cr.usgs.gov/neis/station_book/station_book.html',...
    'Callback','web http://wwwneic.cr.usgs.gov/neis/station_book/station_book.html -browser');

google=uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Pushbutton',...
    'Position',[pos(3)-160 y-102 120 25],...
    'String', 'GoogleEarth Link',...
    'Tooltip','Open current station in GoogleEarth (Windows/MAC only)',...
    'Callback',['googleEarthlink(config.slat,config.slong,',...
    '[config.stnname, '' ('' config.netw '')''],',...
    'fullfile(config.savedir,[config.project(1:end-4), ''.kml'']),'...
    'config.comment);'] );





%% edit fields
x = 160;
y = y+5;
h.station(1) = uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Edit',...
    'BackgroundColor','w',...
    'Position',[x y w v],...
    'String', config.stnname,...
    'FontName','FixedWidth',...
    'Callback','config.stnname=get(gcbo,''String'');');
h.station(2) = uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Edit',...
    'BackgroundColor','w',...
    'Position',[x y-1.75*v w v],...
    'String', config.netw,...
    'FontName','FixedWidth',...
    'Callback','config.netw=get(gcbo,''String'');');
h.station(3) = uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Edit',...
    'BackgroundColor','w',...
    'Position',[x y-3.5*v w v],...
    'String', config.slat,...
    'FontName','FixedWidth',...
    'Callback','config.slat=str2num(get(gcbo,''String''));');
h.station(4) = uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Edit',...
    'BackgroundColor','w',...
    'Position',[x y-5.25*v w v],...
    'String', config.slong,...
    'FontName','FixedWidth',...
    'Callback','config.slong=str2num(get(gcbo,''String''));');


%% Correction
uicontrol('Parent',h.panel(1),'Units','pixel',...
        'Style','text',...
        'Position',[60 120 250 v],...
        'String', 'Misorientation                                 deg clock wise',...
        'HorizontalAlignment','Left'); 
h.station(5) = uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','Edit',...
    'BackgroundColor','w',...
    'Position',[x 125 w v],...
    'String', config.rotation,...
    'FontName','FixedWidth',...
    'Callback','config.rotation=str2num(get(gcbo,''String'')); rotate_seismographIcon') ;


uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','CheckBox',...
    'value',config.SwitchEN,...
    'Position',[60 95 200 18],...
    'String', 'Switch East and North components',...
    'Callback','config.SwitchEN=get(gcbo,''value''); rotate_seismographIcon') ;
 uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','CheckBox',...
    'value',config.signE==-1,...%make logical variable from comparison with -1
    'Position',[60 70 200 18],...
    'Tooltip','East component points West?',...
    'String', 'Change East component sign',...
    'Callback','if get(gcbo,''value''), config.signE=-1; else,  config.signE=1; end; rotate_seismographIcon') ;

uicontrol('Parent',h.panel(1),'Units','pixel',...
    'Style','CheckBox',...
    'value',config.signN==-1,...%make logical variable from comparison with -1
    'Position',[60 45 200 18],...
    'Tooltip','North component points South?',...
    'String', 'Change North component sign',...
    'Callback','if get(gcbo,''value''), config.signN=-1; else, config.signN=1; end; rotate_seismographIcon') ;

ax= axes( 'parent',h.panel(1), ...
    'Units','pixel',...
    'position', [pos(3)-125 25 80 80]);


plot([-1 1 nan 0 0],[0 0 nan 1 -1],'k:');
hold on
text([0 .7], [.7 0], {'N','E'},...
    'color','b',...
    'FontSize',get(gcf,'DefaultUicontrolFontSize')-1,...
    'BackgroundColor','w',...
    'Tag','Seimograph Orientation Labels');
plot([.65 0 0],[0 0 .65],'b', 'Tag','Seimograph Orientation Arrows');
hold off


set(ax,...
    'Ytick',[],'Xtick',[],...
    'Ylim',[-1 1],'Xlim',[-1 1])
xlabel('West - East', 'FontSize',get(gcf,'DefaultUicontrolFontSize')-1);
ylabel('South - North','FontSize',get(gcf,'DefaultUicontrolFontSize')-1);
rotate_seismographIcon

% get(findobj('Tag','Orientation Arrows'),'Xdata')


axes( 'parent',h.panel(1), 'Units','pixel','position', [20 160 60 60] )
image(icon.compas)
axis off
axes( 'parent',h.panel(1), 'Units','pixel','position', [320 160 48 48] )
image(icon.check)
axis off