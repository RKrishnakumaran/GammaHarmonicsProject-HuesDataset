%% New Fig 3: PDs vs harmonic and gamma peak powers (PDstats_Power.m)
destination = './Figures/';
mkdir(destination);
%%
PDstats;
%%
data = deg2rad(alpaPDs_hueALL);
elecind = alpaelecind_elecALL;
hueind = alpahueind_elecALL;
alpaSelectedElectrodes = [];

alpa_lowpower_phases = [];
alpa_highpower_phases = [];
alpa_lowpower_mean_phases = [];
alpa_highpower_mean_phases = [];

for nh = 1:numel(alpa_ghRatio)
    prPD = [];
    % prctile05medPD = [];
    % prctile25medPD = [];
    % prctile50medPD = [];
    % prctile75medPD = [];
    % prctile95medPD = [];
    for elec = 1:numel(alpa_huemean_pds)
        pds = data((elecind == elec)&(hueind == nh));
        prPD = [prPD, circ_rtest(pds)];
        %     prctile05medPD = [prctile05medPD, prctile(medstat(:),05)];
        %     prctile25medPD = [prctile25medPD, prctile(medstat(:),25)];
        %     prctile50medPD = [prctile50medPD, median(medstat(:))];
        %     prctile75medPD = [prctile75medPD, prctile(medstat(:),75)];
        %     prctile95medPD = [prctile95medPD, prctile(medstat(:),95)];
    end
    prPD = ones(size(prPD(:)))==1;% all electrodes selected% prPD(:) < 0.01;
    alpaSelectedElectrodes = [alpaSelectedElectrodes, prPD(:)];
end
%%
data = deg2rad(tutuPDs_hueALL);
elecind = tutuelecind_elecALL;
hueind = tutuhueind_elecALL;
tutuSelectedElectrodes = [];

tutu_lowpower_phases = [];
tutu_highpower_phases = [];
tutu_lowpower_mean_phases = [];
tutu_highpower_mean_phases = [];

for nh = 1:numel(tutu_ghRatio)
    prPD = [];
    % prctile05medPD = [];
    % prctile25medPD = [];
    % prctile50medPD = [];
    % prctile75medPD = [];
    % prctile95medPD = [];
    for elec = 1:numel(tutu_huemean_pds)
        pds = data((elecind == elec)&(hueind == nh));
        prPD = [prPD, circ_rtest(pds)];
        %     prctile05medPD = [prctile05medPD, prctile(medstat(:),05)];
        %     prctile25medPD = [prctile25medPD, prctile(medstat(:),25)];
        %     prctile50medPD = [prctile50medPD, median(medstat(:))];
        %     prctile75medPD = [prctile75medPD, prctile(medstat(:),75)];
        %     prctile95medPD = [prctile95medPD, prctile(medstat(:),95)];
    end
    prPD = ones(size(prPD(:)))==1;% all electrodes selected % prPD(:) < 0.01;
    tutuSelectedElectrodes = [tutuSelectedElectrodes, prPD(:)];
end
%%
%%
% binwidth = 2;
%% load data and/or source figures
% fold = figure;
f6 = figure;
colors = [hsv(36); 0.65*[1,1,1]];
markers = 'xo';
linestyles = ':-';
% nboot = 250; bootfun = @(x)median(x);

% M1
gpowers = alpaGpower_hue;
gpowers_select = [];
data = alpaPD;
hueind = alpahueind_elecALL;
meanPD = [];
stdPD = [];
prPD = [];
CImeanPD = [];
% prctile05medPD = [];
% prctile25medPD = [];
% prctile50medPD = [];
% prctile75medPD = [];
% prctile95medPD = [];
allpds = [];
allhueind = [];
for hue = 1:numel(alpa_ghRatio)
%     pds = data(hueind == hue);
%     pds = (data{hue}); pds = squeeze(circ_mean(pds,[],1)); pds = pds(:,alpaSelectedElectrodes(:,hue)==1); pds = pds(:);
    pds = (data{hue}); pds = squeeze(circ_mean(pds(floor(end/4):ceil(3*end/4),:,:),[],1)); pds = pds(:,alpaSelectedElectrodes(:,hue)==1); pds = pds(:);
    allpds = [allpds; pds(:)];
    allhueind = [allhueind; hue*ones(size(pds(:)))];
    if any(alpaSelectedElectrodes(:,hue)==1) %~isempty(pds)
    meanPD = [meanPD, wrapTo360(rad2deg(circ_mean(pds)))];
    stdPD = [stdPD, rad2deg(circ_std(pds))];
    prPD = [prPD, circ_rtest(pds)];
    CImeanPD = [CImeanPD, rad2deg(circ_confmean(pds(:),0.05))];
        if prPD(end) & (gpowers(hue)<10 & gpowers(hue)>=5)
            alpa_lowpower_phases = [alpa_lowpower_phases ; pds(:)];
        elseif prPD(end) & (gpowers(hue)<15 & gpowers(hue)>=10)
            alpa_highpower_phases = [alpa_highpower_phases ; pds(:)];
        end
%     prctile05medPD = [prctile05medPD, prctile(medstat(:),05)];
%     prctile25medPD = [prctile25medPD, prctile(medstat(:),25)];
%     prctile50medPD = [prctile50medPD, median(medstat(:))];
%     prctile75medPD = [prctile75medPD, prctile(medstat(:),75)];
%     prctile95medPD = [prctile95medPD, prctile(medstat(:),95)];
    else
        meanPD = [meanPD, wrapTo360(rad2deg(circ_mean(pds)))];
        stdPD = [stdPD, rad2deg(circ_std(pds))];
        prPD = [prPD, nan];
		CImeanPD = [CImeanPD, rad2deg(circ_confmean(pds(:),0.05))];
    end
end
prPD = prPD(:) < 0.01;
[sum(prPD), numel(prPD)]
alpaprPD = prPD;
subplot(2,1,1);
% scatter(gpowers(:), medGH(:),[],colors,'marker','o');
% scatter(gpowers(:), meanPD(:),[],colors*1/2 + ((prPD(:)*1/2).*colors),'marker',markers(prPD(:)*1 + 1));
scatter(gpowers(prPD), meanPD(prPD),[],colors(prPD,:),'filled','marker',markers(2));
hold on;
% scatter(gpowers(~prPD), meanPD(~prPD),[],colors(~prPD,:),'marker',markers(1));
scatter(gpowers(~prPD), meanPD(~prPD),[],colors(~prPD,:),'marker',markers(2));

alpa_lowpower_mean_phases = meanPD(prPD' & (gpowers<10 & gpowers >=5));
alpa_highpower_mean_phases = meanPD(prPD' & (gpowers<15 & gpowers >=10));

% ydata = semedGH(:)*[-1,1] + medGH(:);
% ydata = [prctile25medGH(:), prctile75medGH(:)];
% ydata = [prctile05medPD(:), prctile95medPD(:)];
% ydata = stdPD(:)*[-1,1] + meanPD(:);
% for hueind = 1:numel(alpa_ghRatio)
%     x = gpowers(hueind)*[1,1]; y = ydata(hueind,:);
%     line(x,y,'color', colors(hueind,:),'linestyle',linestyles(prPD(hueind)*1+1));
% end
% title('M1')
% ylim([1, max([ylim,3])]);
axis tight;
% xlim([min(xlim)-2, max(xlim)+2]);
ylim([-15, 375]);
yticks(0:60:360); yticklabels(0:60:360);

xlim([0, 25]);
hold on;
ll = plot([0;max([35,xlim])],[180;180], 'color',[0.5,0.5,0.5],'linestyle',':');
% ll = plot([0;max([35,xlim])],[-180;-180], 'color',[0.5,0.5,0.5],'linestyle',':');

%%
% meanPDiff = [];
% stdPD = [];
% CImeanPD = [];
% % prPDiff = [];
% powcenter = [];
% for pow = 0:binwidth:35
% %     hueind = 1:37;
% %     huesel = any(hueind(:) == find((alpaGpower_hue(:)' >= pow) & (alpaGpower_hue(:)' < pow+binwidth)),2);
% % %     ghs = data(huesel);
% %     phdiff = [];
% %     hues = find(huesel);
% %     for nh = hues(:)'
% %         pds = data{nh}; pds = squeeze(circ_mean(pds,[],1)); pds = pds(:,alpaSelectedElectrodes(:,nh)==1);
% %         phdiff = [phdiff, pds(:)'];
% %     end
% %     hueindex = prPD(:) & (gpowers(:) >= pow) & (gpowers(:) < pow+binwidth);
% %     phdiff = deg2rad(meanPD(hueindex));
%     hueindex = find(prPD(:) & (gpowers(:) >= pow) & (gpowers(:) < pow+binwidth));
%     phdiff = (allpds(any(allhueind(:) == hueindex(:)',2)));
%     if ~isempty(phdiff)
%         powcenter = [powcenter, pow+binwidth/2];
%         meanPDiff = [meanPDiff, wrapTo360(rad2deg(circ_mean(phdiff(:))))];
%         stdPD = [stdPD, rad2deg(circ_std(phdiff(:)))];
%         CImeanPD = [CImeanPD, rad2deg(circ_confmean(phdiff(:),0.05))];
% %         prPDiff = [prPDiff, circ_rtest(phdiff(:))];
%     end
% end
% prPDiff = prPDiff < 0.01;
% d1M1 = scatter(powcenter, meanPDiff,[],'k','filled','marker','d');
% d1M1 = scatter(powcenter(prPDiff), meanPDiff(prPDiff),[],'k','filled','marker','d');
% hold on;
% d2M1 = scatter(powcenter(~prPDiff), meanPDiff(~prPDiff),[],'k','marker','d');
ydata = CImeanPD(:)*[-1,1] + meanPD(:);
for ind = 1:numel(gpowers)
    x = gpowers(ind)*[1,1]; y = ydata(ind,:);
%     line(x,y,'color', colors(ind,:));
    line(x,y,'color', 'k','linestyle','-');
end
%%
xlabel('Electrode-averaged Gamma power (dB)')
ylabel('Phase Difference (\circ)');
curax = subplot(2,1,1); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M1', 'fontweight','bold','fontsize',20,'edgecolor','none');
%% estimating ghratio and testing them for harmonic
alpameanPD = wrapTo360(rad2deg(circ_mean(allpds(:))));
alpaCImeanPD = wrapTo360(rad2deg(circ_confmean(allpds(:),0.05)));
alpasemeanPD = wrapTo360(rad2deg(circ_std(allpds(:))));
alpaPDrtest = circ_rtest(allpds(:));

meanPD = meanPD(:);
alpameanPDhues_untested = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD))));
alpaCImeanPDhues_untested = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD),0.05)));
[sortedgpowers,sortedgpowerIndices] = sort(gpowers,'descend');
alpameanPD5hues_untested = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD(sortedgpowerIndices(1:5))))));
gpowers2 = gpowers;
gpowers2(gpowers>10) = -inf;
[sortedgpowers2,sortedgpowerIndices2] = sort(gpowers2,'descend');
alpameanPD5hues2_untested = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD(sortedgpowerIndices2(1:5))))));

alpaCImeanPD5hues_untested = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD(sortedgpowerIndices(1:5))),0.05)));
alpaCImeanPD5hues2_untested = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD(sortedgpowerIndices2(1:5))),0.05)));

meanPD = meanPD(prPD); meanPD = meanPD(:);
alpameanPDhues = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD))));
alpaCImeanPDhues = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD),0.05)));
selgpowers = gpowers(prPD);
[sortedgpowers,sortedgpowerIndices] = sort(selgpowers,'descend');
alpameanPD5hues = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD(sortedgpowerIndices(1:5))))));
alpaCImeanPD5hues = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD(sortedgpowerIndices(1:5))),0.05)));


alpasemeanPDhues = wrapTo360(rad2deg(circ_std(deg2rad(meanPD))));
alpaPDrtesthues = wrapTo360(circ_rtest(deg2rad(meanPD)));
% meanPD = meanPDiff(:);
% alpameanPDbinnedhues = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD))));
% alpaCImeanPDbinnedhues = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD),0.05)));
% alpasemeanPDbinnedhues = wrapTo360(rad2deg(circ_std(deg2rad(meanPD))));
% alpaPDrtestbinnedhues = wrapTo360(circ_rtest(deg2rad(meanPD)));

%%
% M2
gpowers = tutuGpower_hue;
data = tutuPD;
hueind = tutuhueind_elecALL;
meanPD = [];
prPD = [];
stdPD = [];
CImeanPD = [];
% prctile05medPD = [];
% prctile25medPD = [];
% prctile50medPD = [];
% prctile75medPD = [];
% prctile95medPD = [];
allpds = [];
allhueind = [];
for hue = 1:numel(tutu_ghRatio)
%     pds = data(hueind == hue);
%     pds = (data{hue}); pds = squeeze(circ_mean(pds,[],1)); pds = pds(:,tutuSelectedElectrodes(:,hue)==1); pds = pds(:);
    pds = (data{hue}); pds = squeeze(circ_mean(pds(floor(end/4):ceil(3*end/4),:,:),[],1)); pds = pds(:,tutuSelectedElectrodes(:,hue)==1); pds = pds(:);
    allpds = [allpds; pds(:)];
    allhueind = [allhueind; hue*ones(size(pds(:)))];
    if any(tutuSelectedElectrodes(:,hue)==1)%~isempty(pds)
        meanPD = [meanPD, wrapTo360(rad2deg(circ_mean(pds)))];
        stdPD = [stdPD, rad2deg(circ_std(pds))];
        prPD = [prPD, circ_rtest(pds)];
        CImeanPD = [CImeanPD, rad2deg(circ_confmean(pds(:),0.05))];
        
        if prPD(end) & (gpowers(hue)<17.5 & gpowers(hue)>=12.5)
            tutu_lowpower_phases = [tutu_lowpower_phases ; pds(:)];
        elseif prPD(end) & (gpowers(hue)<22.5 & gpowers(hue)>=17.5)
            tutu_highpower_phases = [tutu_highpower_phases ; pds(:)];
        end
        
    else
        meanPD = [meanPD, wrapTo360(rad2deg(circ_mean(pds)))];
        stdPD = [stdPD, rad2deg(circ_std(pds))];
        prPD = [prPD, nan];
		CImeanPD = [CImeanPD, rad2deg(circ_confmean(pds(:),0.05))];
    end
end
prPD = prPD(:) < 0.01;
[sum(prPD), numel(prPD)]
subplot(2,1,2);
% scatter(gpowers(:), medGH(:),[],colors,'marker','o');
% scatter(gpowers(:), meanPD(:),[],colors*1/2 + ((prPD(:)*1/2).*colors),'marker',markers(prPD(:)*1+1));
scatter(gpowers(prPD), meanPD(prPD),[],colors(prPD,:),'filled','marker',markers(2));
hold on;
% scatter(gpowers(~prPD), meanPD(~prPD),[],colors(~prPD,:),'marker',markers(1));
scatter(gpowers(~prPD), meanPD(~prPD),[],colors(~prPD,:),'marker',markers(2));

tutu_lowpower_mean_phases = meanPD(prPD' & (gpowers<17.5 & gpowers >=12.5));
tutu_highpower_mean_phases = meanPD(prPD' & (gpowers<22.5 & gpowers >=17.5));


% ydata = semedGH(:)*[-1,1] + medGH(:);
% ydata = [prctile25medGH(:), prctile75medGH(:)];
% ydata = [prctile05medPD(:), prctile95medPD(:)];
% ydata = stdPD(:)*[-1,1] + meanPD(:);
% for hueind = 1:numel(tutu_ghRatio)
%     x = gpowers(hueind)*[1,1]; y = ydata(hueind,:);
%     line(x,y,'color', colors(hueind,:), 'linestyle',linestyles(prPD(hueind)*1+1));
% end
% title('M2')
% ylim([1, max([ylim,3])]);
axis tight;
% xlim([min(xlim)-2, max(xlim)+2]);
ylim([-15, 375]);
yticks(0:60:360); yticklabels(0:60:360);

xlim([0, 25]);
hold on;
ll = plot([0;max([35,xlim])],[180;180], 'color',[0.5,0.5,0.5],'linestyle',':');
% ll = plot([0;max([35,xlim])],[-180;-180], 'color',[0.5,0.5,0.5],'linestyle',':');
% %%
% meanPDiff = [];
% stdPD = [];
% CImeanPD = [];
% % prPDiff = [];
% powcenter = [];
% for pow = 0:binwidth:35
% %     hueind = 1:37;
% %     huesel = any(hueind(:) == find((tutuGpower_hue(:)' >= pow) & (tutuGpower_hue(:)' < pow+binwidth)),2);
% % %     ghs = data(huesel);
% %     phdiff = [];
% %     hues = find(huesel);
% %     for nh = hues(:)'
% %         pds = data{nh}; pds = squeeze(circ_mean(pds,[],1)); pds = pds(:,tutuSelectedElectrodes(:,nh)==1);
% %         phdiff = [phdiff, pds(:)'];
% %     end
% %     hueindex = prPD(:) & (gpowers(:) >= pow) & (gpowers(:) < pow+binwidth);
% %     phdiff = deg2rad(meanPD(hueindex));
%     hueindex = find(prPD(:) & (gpowers(:) >= pow) & (gpowers(:) < pow+binwidth));
%     phdiff = (allpds(any(allhueind(:) == hueindex(:)',2)));
%     if ~isempty(phdiff)
%         powcenter = [powcenter, pow+binwidth/2];
%         meanPDiff = [meanPDiff, wrapTo360(rad2deg(circ_mean(phdiff(:))))];
%         stdPD = [stdPD, rad2deg(circ_std(phdiff(:)))];
% 		CImeanPD = [CImeanPD, rad2deg(circ_confmean(phdiff(:),0.05))];
% 
% %         prPDiff = [prPDiff, circ_rtest(phdiff(:))];
%     end
% end
% prPDiff = prPDiff < 0.01;
% d1M2 = scatter(powcenter, meanPDiff,[],'k','filled','marker','d');

% d1M2 = scatter(powcenter(prPDiff), meanPDiff(prPDiff),[],'k','filled','marker','d');
% hold on;
% d2M2 = scatter(powcenter(~prPDiff), meanPDiff(~prPDiff),[],'k','marker','d');
ydata = CImeanPD(:)*[-1,1] + meanPD(:);
for ind = 1:numel(gpowers)
    x = gpowers(ind)*[1,1]; y = ydata(ind,:);
%     line(x,y,'color', colors(ind,:));
    line(x,y,'color', 'k','linestyle','-');
end
%%
xlabel('Electrode-averaged Gamma power (dB)')
ylabel('Phase Difference (\circ)');
curax = subplot(2,1,2); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M2', 'fontweight','bold','fontsize',20,'edgecolor','none');
%%
sgtitle('Phase Difference distributions across stimuli','fontsize',20,'fontweight','bold');
%% estimating ghratio and testing them for harmonic
tutumeanPD = wrapTo360(rad2deg(circ_mean(allpds(:))));
tutuCImeanPD = wrapTo360(rad2deg(circ_confmean(allpds(:),0.05)));
tutusemeanPD = wrapTo360(rad2deg(circ_std(allpds(:))));
tutuPDrtest = circ_rtest(allpds(:));

meanPD = meanPD(:);
tutumeanPDhues_untested = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD))));
tutuCImeanPDhues_untested = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD),0.05)));
[sortedgpowers,sortedgpowerIndices] = sort(gpowers,'descend');

tutumeanPD5hues_untested = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD(sortedgpowerIndices(1:5))))));
gpowers2 = gpowers;
gpowers2(gpowers>15) = -inf;
[sortedgpowers2,sortedgpowerIndices2] = sort(gpowers2,'descend');
tutumeanPD5hues2_untested = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD(sortedgpowerIndices2(1:5))))));

tutuCImeanPD5hues_untested = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD(sortedgpowerIndices(1:5))),0.05)));
tutuCImeanPD5hues2_untested = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD(sortedgpowerIndices2(1:5))),0.05)));

meanPD = meanPD(prPD); meanPD = meanPD(:);
tutumeanPDhues = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD))));
tutuCImeanPDhues = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD),0.05)));
selgpowers = gpowers(prPD);
[sortedgpowers,sortedgpowerIndices] = sort(selgpowers,'descend');
tutumeanPD5hues = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD(sortedgpowerIndices(1:5))))));
tutuCImeanPD5hues = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD(sortedgpowerIndices(1:5))),0.05)));


tutusemeanPDhues = wrapTo360(rad2deg(circ_std(deg2rad(meanPD))));
tutuPDrtesthues = wrapTo360(circ_rtest(deg2rad(meanPD)));
% meanPD = meanPDiff(:);
% tutumeanPDbinnedhues = wrapTo360(rad2deg(circ_mean(deg2rad(meanPD))));
% tutuCImeanPDbinnedhues = wrapTo360(rad2deg(circ_confmean(deg2rad(meanPD),0.05)));
% tutusemeanPDbinnedhues = wrapTo360(rad2deg(circ_std(deg2rad(meanPD))));
% tutuPDrtestbinnedhues = wrapTo360(circ_rtest(deg2rad(meanPD)));
%%
f = f6; f.WindowState='maximized';   postformatFig;
% d1M1.SizeData = round(d1M1.SizeData*3/4);
% d2M1.SizeData = round(d2M1.SizeData*3/4);

% d1M2.SizeData = round(d1M2.SizeData*3/4);
% d2M2.SizeData = round(d2M2.SizeData*3/4);
%%
savefig(gcf, [destination, 'Fig6_Draft8_blackErrorbars.fig']);
print([destination, 'Fig6_Draft8_blackErrorbars.png'],'-dpng','-r480');
mkdir([destination, '/Draft 9/']);
print([destination, '/Draft 9/Fig5.tif'],'-dtiff','-r480');

%%
disp('Rayleigh tested')
disp(['M1 mean of mean-phases/hue: ',num2str(alpameanPDhues),'+/-', num2str(alpaCImeanPDhues)])%,'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 mean of mean-phases/hue: ',num2str(tutumeanPDhues),'+/-', num2str(tutuCImeanPDhues)])%,'; p=',num2str(tutuGHratio_signRank_hues)])
% 
% disp(['M1 median of median-FreqRatio/top5hue: ',num2str(alpamedianGHratio_5Powerfulhues),'+/-', num2str(alpasemedianGHratio_5Powerfulhues)])
% disp(['M2 median of median-FreqRatio/top5hue: ',num2str(tutumedianGHratio_5Powerfulhues),'+/-', num2str(tutusemedianGHratio_5Powerfulhues)])

disp(['M1 mean of mean-phases/top10hue: ',num2str(alpameanPD5hues),'+/-', num2str(alpaCImeanPD5hues)])%,'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 mean of mean-phases/top10hue: ',num2str(tutumeanPD5hues),'+/-', num2str(tutuCImeanPD5hues)])%,'; p=',num2str(tutuGHratio_signRank_hues)])
% disp(['M1 median of median-FreqRatio/top15hue: ',num2str(alpamedianGHratio_15Powerfulhues),'+/-', num2str(alpasemedianGHratio_15Powerfulhues),'; p=',num2str(alpaGHratio_signRank_15Powerfulhues)])
% disp(['M2 median of median-FreqRatio/top15hue: ',num2str(tutumedianGHratio_15Powerfulhues),'+/-', num2str(tutusemedianGHratio_15Powerfulhues),'; p=',num2str(tutuGHratio_signRank_15Powerfulhues)])
%%
disp('Untested')
disp(['M1 mean of mean-phases/hue: ',num2str(alpameanPDhues_untested),'+/-', num2str(alpaCImeanPDhues_untested)])%,'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 mean of mean-phases/hue: ',num2str(tutumeanPDhues_untested),'+/-', num2str(tutuCImeanPDhues_untested)])%,'; p=',num2str(tutuGHratio_signRank_hues)])
% 
% disp(['M1 median of median-FreqRatio/top5hue: ',num2str(alpamedianGHratio_5Powerfulhues),'+/-', num2str(alpasemedianGHratio_5Powerfulhues)])
% disp(['M2 median of median-FreqRatio/top5hue: ',num2str(tutumedianGHratio_5Powerfulhues),'+/-', num2str(tutusemedianGHratio_5Powerfulhues)])

disp(['M1 mean of mean-phases/top5hue: ',num2str(alpameanPD5hues_untested),'+/-', num2str(alpaCImeanPD5hues_untested)])%,'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 mean of mean-phases/top5hue: ',num2str(tutumeanPD5hues_untested),'+/-', num2str(tutuCImeanPD5hues_untested)])%,'; p=',num2str(tutuGHratio_signRank_hues)])
disp(['M1 mean of mean-phases/next_5hue: ',num2str(alpameanPD5hues2_untested),'+/-', num2str(alpaCImeanPD5hues2_untested)])%,'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 mean of mean-phases/next_5hue: ',num2str(tutumeanPD5hues2_untested),'+/-', num2str(tutuCImeanPD5hues2_untested)])%,'; p=',num2str(tutuGHratio_signRank_hues)])
% disp(['M1 median of median-FreqRatio/top15hue: ',num2str(alpamedianGHratio_15Powerfulhues),'+/-', num2str(alpasemedianGHratio_15Powerfulhues),'; p=',num2str(alpaGHratio_signRank_15Powerfulhues)])
% disp(['M2 median of median-FreqRatio/top15hue: ',num2str(tutumedianGHratio_15Powerfulhues),'+/-', num2str(tutusemedianGHratio_15Powerfulhues),'; p=',num2str(tutuGHratio_signRank_15Powerfulhues)])




%%
disp('############################')

%alpa
alpa_lpp = [wrapTo360(rad2deg(circ_mean(alpa_lowpower_phases))), wrapTo360(rad2deg(circ_confmean(alpa_lowpower_phases)))]; 
alpa_hpp = [wrapTo360(rad2deg(circ_mean(alpa_highpower_phases))), wrapTo360(rad2deg(circ_confmean(alpa_highpower_phases)))];

alpa_lpmp = [wrapTo360(rad2deg(circ_mean(deg2rad(alpa_lowpower_mean_phases')))), wrapTo360(rad2deg(circ_confmean(deg2rad(alpa_lowpower_mean_phases'))))];
alpa_hpmp = [wrapTo360(rad2deg(circ_mean(deg2rad(alpa_highpower_mean_phases')))), wrapTo360(rad2deg(circ_confmean(deg2rad(alpa_highpower_mean_phases'))))];

disp(['M1 mean of all-phases/lowpower(5-10 dB): ',num2str(alpa_lpp(1)),'+/-', num2str(alpa_lpp(2))]);
disp(['M1 mean of all-phases/highpower(10-15 dB): ',num2str(alpa_hpp(1)),'+/-', num2str(alpa_hpp(2))]);

%tutu
tutu_lpp = [wrapTo360(rad2deg(circ_mean(tutu_lowpower_phases))), wrapTo360(rad2deg(circ_confmean(tutu_lowpower_phases)))]; 
tutu_hpp = [wrapTo360(rad2deg(circ_mean(tutu_highpower_phases))), wrapTo360(rad2deg(circ_confmean(tutu_highpower_phases)))];

tutu_lpmp = [wrapTo360(rad2deg(circ_mean(deg2rad(tutu_lowpower_mean_phases')))), wrapTo360(rad2deg(circ_confmean(deg2rad(tutu_lowpower_mean_phases'))))];
tutu_hpmp = [wrapTo360(rad2deg(circ_mean(deg2rad(tutu_highpower_mean_phases')))), wrapTo360(rad2deg(circ_confmean(deg2rad(tutu_highpower_mean_phases'))))];

disp(['M2 mean of all-phases/lowpower(12.5-17.5 dB): ',num2str(tutu_lpp(1)),'+/-', num2str(tutu_lpp(2))]);
disp(['M2 mean of all-phases/highpower(17.5-22.5 dB): ',num2str(tutu_hpp(1)),'+/-', num2str(tutu_hpp(2))]);

disp('############################')

disp(['M1 mean of mean-phases/lowpower(5-10 dB): ',num2str(alpa_lpmp(1)),'+/-', num2str(alpa_lpmp(2))]);
disp(['M1 mean of mean-phases/highpower(10-15 dB): ',num2str(alpa_hpmp(1)),'+/-', num2str(alpa_hpmp(2))]);

disp(['M2 mean of mean-phases/lowpower(12.5-17.5 dB): ',num2str(tutu_lpmp(1)),'+/-', num2str(tutu_lpmp(2))]);
disp(['M2 mean of mean-phases/highpower(17.5-22.5 dB): ',num2str(tutu_hpmp(1)),'+/-', num2str(tutu_hpmp(2))]);


