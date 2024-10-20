function splitdiagnosticplot(Q, T, extime, L, E, N, inc, bazi,sampling, maxtime, pol,...
    phiRC, dtRC, Cmatrix, corFSrc, QTcorRC,...
    phiSC, dtSC, Ematrix, corFSsc, QTcorSC,...
    phiEV, dtEV, LevelSC, LevelRC, LevelEV, splitoption)

% display the results of a Rotation-Correaltion and a minimum Energy
% splitting procedure in a single plot
% Inputs are expected in the following order:
%     Q, T
%     E, N
%     inclination, backazimuth, sampling [sec]
%     phiRC
%     dtRC
%     Cmatrix
%     corFastSlowparticleRC - corrected RC particle motion [F S]
%     phiSC
%     dtSC
%     pol: initial polarisation
%     Ematrix
%     sampling
%     corFastSlowparticleSC - corrected SC particle motion [F S]
%     Phi_errorSC   - SC fast axis estimation error interval
%     dt_errorSC    - SC delay time estimation error interval
%     Level         - confidence level for Silver&Chan Energy map

% Andreas W�stefeld, 12.03.06
global thiseq config

Synfig = findobj('name', 'Diagnostic Viewer','type','figure');
if isempty(Synfig)
    S = get(0,'Screensize');
    Synfig = figure('name', 'Diagnostic Viewer',...
        'Renderer',        'painters',...
        'Color',           'w',...
        'NumberTitle',     'off',...
        'MenuBar',         'none',...
        'PaperType',       config.PaperType,...
        'PaperOrientation','landscape',...
        'PaperUnits',      'centimeter',...
        'position',        [.01*S(3) .1*S(4) .98*S(3) .75*S(4)]);
else
    figure(Synfig)
    clf
    set(Synfig,'PaperOrientation','landscape',...
        'PaperType',       config.PaperType)
end
orient landscape
colormap(gray)

fontsize = get(gcf,'DefaultAxesFontsize')-1;
titlefontsize = fontsize+2;





[axH, axRC, axSC, axSeis] = splitdiagnosticLayout(Synfig);
splitdiagnosticSetHeader(axH, phiRC, dtRC, phiSC, dtSC, phiEV, dtEV, pol, splitoption)


switch splitoption
    case 'Minimum Energy'
        Ematrix = Ematrix(:,:,1);
        optionstr ='Minimum Energy';
        phi = phiSC(2);
        dt  = dtSC(2);
        Level = LevelSC;
        Maptitle = 'Energy Map of T';
    case 'Eigenvalue: max(lambda1 / lambda2)'
        Ematrix = Ematrix(:,:,2);
        optionstr ='Maximum   \lambda_1 / \lambda_2';
        phi = phiEV(2);
        dt  = dtEV(2);
        Level =LevelEV;
        Maptitle = 'Map of Eigenvalues \lambda_1 / \lambda_2';
    case 'Eigenvalue: min(lambda2)'
        Ematrix = Ematrix(:,:,2);
        optionstr ='Minimum  \lambda_2';
        phi = phiEV(2);
        dt  = dtEV(2);
        Level =LevelEV;
        Maptitle = 'Map of Eigenvalue \lambda_2';
        
    case 'Eigenvalue: max(lambda1)'
        Ematrix = Ematrix(:,:,2);
        optionstr ='Maximum  \lambda_1';
        phi = phiEV(2);
        dt  = dtEV(2);
        Level =LevelEV;
        Maptitle = 'Map of Eigenvalue \lambda_1';

    case 'Eigenvalue: min(lambda1 * lambda2)'
        Ematrix = Ematrix(:,:,2);
        optionstr ='Minimum   \lambda_1 * \lambda_2';
        phi = phiEV(2);
        dt  = dtEV(2);
        Level =LevelEV;
        Maptitle = 'Map of Eigenvalues \lambda_1 * \lambda_2';
end


%% rotate seismograms for plots (backwards == counter-clockwise => use transposed matrix M)
M = rot3D(inc, bazi);

ZEN = M' *[L,  QTcorRC]';
Erc = ZEN(2,:); 
Nrc = ZEN(3,:);

ZEN = M' *[L,  QTcorSC]';
Esc = ZEN(2,:); 
Nsc = ZEN(3,:);

s = size(QTcorRC,1);%selection length




%% x-values for seismogram plots
t = (0:(s-1))*sampling;

%%  Rotation-Correlation Method% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast/slow seismograms
axes(axRC(1))
sumFS1 = sum(abs( corFSrc(:,1) -corFSrc(:,2)));
sumFS2 = sum(abs(-corFSrc(:,1) -corFSrc(:,2)));
if ( sumFS1 < sumFS2 )
    sig = 1;
else
    sig = -1;
end
m1 = max(abs( corFSrc(:,1)));
m2 = max(abs( corFSrc(:,2)));
plot(t, corFSrc(:,1)/m1,'b--',   t,sig*corFSrc(:,2)/m2,'r-','LineWidth',1);
xlim([t(1) t(end)])
title(['corrected Fast (' char([183 183]) ') & Slow(-)'],'FontSize',titlefontsize);
set(gca,'Ytick' , [-1 0 1])
ylabel('Rotation-Correlation','FontSize',titlefontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corrected seismograms
axes(axRC(2))
plot(t, QTcorRC(:,1),'b--',    t, QTcorRC(:,2) ,'r-','LineWidth',1);
title([' corrected Q(' char([183 183]) ') & T(-)'],'FontSize',titlefontsize);
xlim([t(1) t(end)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface Particle motion
axes(axRC(3))
plot(E, N, 'b--', Erc, Nrc,'r-','LineWidth',1);
xlabel('\leftarrowW - E\rightarrow', 'Fontsize',fontsize-1);
ylabel('\leftarrowS - N\rightarrow', 'Fontsize',fontsize-1);
title(['Particle motion before (' char([183 183]) ') & after (-)'],'FontSize',titlefontsize);
axis equal

tmp = max([abs(xlim) abs(ylim)]);%set [0 0] to centre of plot
set(gca, 'xlim', [-tmp tmp], 'ylim', [-tmp tmp], 'XtickLabel',[], 'YtickLabel',[])
set(gca, 'Ytick', get(gca,'Xtick'))
hold on
X = sin(bazi/180*pi)*tmp;
Y = cos(bazi/180*pi)*tmp;
plot( [-X X], [-Y Y], 'k:' )
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Correlation Map
axes(axRC(4))   %
hold on
f  = size(Cmatrix);
ts = linspace(0,maxtime,f(2));
ps = linspace(-90,90,f(1));

maxi = max(Cmatrix(:));% allways <=  1 since correlation coeffcient (^5)
mini = min(Cmatrix(:));% allways >= -1
maxmin = abs(mini - maxi)/2;% allways between 0 and 1

nb_contours = 12;floor((1 - maxmin)*9);
[C, h] = contourf('v6',ts,ps,-Cmatrix,-[LevelRC LevelRC]);
contour(ts, ps, Cmatrix, nb_contours);



B = mod(bazi,90);
plot([0 0]+sampling, [B B-90],'k>','markersize',5,'linewidth',1,'MarkerFaceColor','k' )
plot([maxtime maxtime]-sampling, [B B-90],'k<','markersize',5,'linewidth',1,'MarkerFaceColor','k' )
line([dtRC(2) dtRC(2)],[-90 90],'Color',[0 0 1])
line([0 maxtime], [phiRC(2) phiRC(2)],'Color',[0 0 1])
title('Map of Correlation Coefficient','FontSize',titlefontsize);
%xlabel('dt [s]', 'Fontsize',fontsize-1);
ylabel('fast axis', 'Fontsize',fontsize-1)
label = ['0' sprintf('|%u',1:maxtime) 'sec'];
set(gca, 'Xtick',[0:1:maxtime], 'XtickLabel', label ,'Ytick',[-90:30:90],'xMinorTick','on','yminorTick','on')
axis([ts(1) ts(end) -90 90])
set(h,'FaceColor',[1 1 1]*.90,'EdgeColor','k','linestyle','-','linewidth',1)


hold off



%%  Silver & Chan
% fast/slow seismograms
axes(axSC(1))
sumFS1 = sum(abs( corFSsc(:,1) -corFSsc(:,2)));
sumFS2 = sum(abs(-corFSsc(:,1) -corFSsc(:,2)));
if ( sumFS1 < sumFS2 )
    sig = 1;
else
    sig = -1;
end
m1 = max(abs( corFSsc(:,1)));
m2 = max(abs( corFSsc(:,2)));
plot(  t, corFSsc(:,1)/m1,'b--',    t, sig*corFSsc(:,2)/m2 ,'r-','LineWidth',1);
xlim([t(1) t(end)])
ylim([-1 1])
title(['corrected Fast (' char([183 183]) ') & Slow(-)'],'FontSize',titlefontsize);
set(gca,'Ytick' , [-1 0 1])
ylabel(optionstr,'FontSize',titlefontsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corrected seismograms (in Ray-system)
axes(axSC(2))
plot(t, QTcorSC(:,1),'b--',    t, QTcorSC(:,2) ,'r-','LineWidth',1);
title([' corrected Q(' char([183 183]) ') & T(-)'],'FontSize',titlefontsize);
xlim([t(1) t(end)])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface Particle motion
axes(axSC(3))
hold on
plot(E, N, 'b--', Esc, Nsc,'r-','LineWidth',1);
xlabel('\leftarrowW - E\rightarrow', 'Fontsize',fontsize-1);
ylabel('\leftarrowS - N\rightarrow', 'Fontsize',fontsize-1);
title(['Particle motion before (' char([183 183]) ') & after (-)'],'FontSize',titlefontsize);
axis equal

tmp = max([abs(xlim) abs(ylim)]);%set [0 0] to centre of plot
set(gca, 'xlim', [-tmp tmp], 'ylim', [-tmp tmp], 'XtickLabel',[], 'YtickLabel',[])
set(gca, 'Ytick', get(gca,'Xtick'))
hold on
X = sin(bazi/180*pi)*tmp;
Y = cos(bazi/180*pi)*tmp;
plot( [-X X], [-Y Y], 'k:' )
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energy Map
axes(axSC(4))
hold on
f  = size(Ematrix);
ts = linspace(0,maxtime,f(2));
ps = linspace(-90,90,f(1));


maxi = max(abs(Ematrix(:)));
mini = min(abs(Ematrix(:)));
nb_contours = floor((1 - mini/maxi)*10);
[C, h] = contourf('v6',ts,ps,-Ematrix,-[Level Level]);
contour(ts, ps, Ematrix, nb_contours);




B = mod(bazi,90);%backazimuth lines
plot([0 0]+sampling, [B B-90],'k>','markersize',5,'linewidth',1,'MarkerFaceColor','k' )
plot([maxtime maxtime]-sampling, [B B-90],'k<','markersize',5,'linewidth',1,'MarkerFaceColor','k' )
line([0 maxtime], [phi phi],'Color',[0 0 1])
line([dt dt],[-90 90],'Color',[0 0 1])



hold off
axis([0 maxtime -90 90])
set(gca, 'Xtick',[0:1:maxtime], 'XtickLabel', label ,'Ytick',[-90:30:90],'xMinorTick','on','yminorTick','on')
%xlabel('dt [s]', 'Fontsize',fontsize-1);
ylabel('fast axis', 'Fontsize',fontsize-1)
title(Maptitle,'FontSize',titlefontsize);
set(h,'FaceColor',[1 1 1]*.90,'EdgeColor','k','linestyle','-','linewidth',1)




%% plot Initial seismograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axSeis(1))
t2 = (0:length(Q)-1)*sampling - extime;
xx  = [0 0 s s]*sampling;
yy  = [0 0 1 1 ];
tmp = fill(xx, yy, [1 1 1]*.90, 'EdgeColor','None');% Selection marker

hold on
plot(t2, Q, 'b--', t2, T, 'r-','LineWidth',1)

tt = thiseq.phase.ttimes;
A  = thiseq.a-extime;
F  = thiseq.f+extime;
tt = tt(A<=tt& tt<=F); %phase arrival within selection
T  = [tt;tt];
T  = T-thiseq.a;
yy = repmat(ylim',size(tt));
plot(T,yy,'k:')

hold off

% title({['Before correction: Q(' char([183 183]) ') & T(-)']},'FontSize',fontsize);
xlim([t2(1) t2(end)])
yy = [ylim fliplr(ylim)];
set(tmp,'yData',yy)
set(axSeis,'Layer','Top')


%% plot stereoplot of current measurement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%,'EraseMode','Xor'
if license('checkout', 'MAP_Toolbox')
    axes(axSeis(2))
    b = repmat(bazi,1,3);
    I = repmat(thiseq.tmpInclination,1,3);
    [h, m] = stereoplot(b, I , [phiRC(2) phiSC(2) phiEV(2)], [dtRC(2) dtSC(2) dtEV(2)]);
    set(h(1), 'Color',[0 .6 0])
    set(h(2), 'Color',[1 0 0])
    set(h(3), 'Color',[0 0 1])
    L = axis;

    text(0, L(4),['  Inc = \bf' num2str(thiseq.tmpInclination,'%4.1f') char(186)], ...
        'Fontname','Fixedwidth', 'VerticalAlignment','top','HorizontalAlignment','center')
else
    delete(axSeis(2))
end
%% EOF %%