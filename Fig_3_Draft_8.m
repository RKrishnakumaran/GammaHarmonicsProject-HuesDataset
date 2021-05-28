%% New Fig 3: GHratio vs harmonic and gamma peak powers (PDstats_Power.m)
destination = './Figures/';
mkdir(destination);
%% load data and/or source figures
% fold = figure;
%%
if ~exist('tutup_wwtest','var')
    PDstats_Power_bootstrapWWtest;
%     PDstats_Power_bootstrapWWtest_oldData;
    close all;
end
%% Fig 3 median and stderr of GH ratios vs Gamma power - each hue
f3 = figure;
colors = [hsv(36);0.65*[1,1,1]];
nboot = 1000; bootfun = @(x)median(x);
binwidth = 2;

% M1
gpowers = alpaGpower_hue;
data = alpa_ghRatio;
% hueind = alpahueind_elecALL;
medGHhue = [];
allGHhue = [];
hueindsGHhue = [];
semedGH = [];
p_sign = [];
p_srank = [];
% prctile05medGH = [];
% prctile25medGH = [];
% prctile50medGH = [];
% prctile75medGH = [];
% prctile95medGH = [];
for hue = 1:numel(alpa_ghRatio)
    %     ghs = data(hueind == hue);
    ghs = data{hue};
    medGHhue = [medGHhue, median(ghs)];
    allGHhue = [allGHhue; ghs(:)];
    hueindsGHhue = [hueindsGHhue; hue * ones(size(ghs(:)))];
    
    if numel(ghs) > 1
        [se, ~] = getSEMedian(ghs(:));
    else
        se = 0;
    end
    semedGH = [semedGH, se];
    pval = signtest(ghs(:),2);
    p_sign = [p_sign, pval];
    pval = signrank(ghs(:),2);
    p_srank = [p_srank, pval];
    %     [ci, medstat] = bootci(nboot, {bootfun, ghs(:)}, 'type','norm','alpha', 0.01);
    %     semedGH = [semedGH, std(medstat(:))];
    %     prctile05medGH = [prctile05medGH, prctile(medstat(:),05)];
    %     prctile25medGH = [prctile25medGH, prctile(medstat(:),25)];
    %     prctile50medGH = [prctile50medGH, median(medstat(:))];
    %     prctile75medGH = [prctile75medGH, prctile(medstat(:),75)];
    %     prctile95medGH = [prctile95medGH, prctile(medstat(:),95)];
end

% subplot(2,1,1);
% scatter(gpowers(:), medGH(:),[],colors,'marker','d');
% hold on;

% medianGH = [];
% medGH = [];
% semedGH = [];
% prctile05medGHalt = [];
% prctile95medGHalt = [];
% prctile05medGH = [];
% prctile25medGH = [];
% prctile50medGH = [];
% prctile75medGH = [];
% prctile95medGH = [];
% powcenter = [];
% for pow = 0:binwidth:35
%     %     hueind = 1:36;
%     %     huesel = any(hueind(:) == find((alpaGpower_hue(:)' >= pow) & (alpaGpower_hue(:)' < pow+2.5)),2);
%     % %     ghs = data(huesel);
%     %     ghs = [];
%     %     hues = find(huesel);
%     %     for nh = hues(:)'
%     %         ghs = [ghs, data{nh}(:)'];
%     %     end
% %     hueindex = (gpowers(:) >= pow) & (gpowers(:) < pow+binwidth);
% %     ghs = medGHhue(hueindex);
%     hueindex = find((gpowers(:) >= pow) & (gpowers(:) < pow+binwidth));
%     ghs = allGHhue(any(hueindsGHhue(:) == hueindex(:)',2));
%     if ~isempty(ghs)
%         powcenter = [powcenter, pow+binwidth/2];
%         %         medGH = [medGH, median(ghs)];
%         %         [ci, medstat] = bootci(nboot, {bootfun, ghs(:)}, 'type','norm','alpha', 0.01);
%         %         semedGH = [semedGH, std(medstat(:))];
%         
%         medianGH = [medianGH, median(ghs(:))];
%         if numel(ghs) > 1
%             [se, ~] = getSEMedian(ghs(:));
%         else
%             se = 0;
%         end
%         semedGH = [semedGH, se];
%         %         prctile05medGHalt = [prctile05medGHalt, prctile(medstat(:),05)];
%         % %         prctile25medGH = [prctile25medGH, prctile(medstat(:),25)];
%         % %         prctile50medGH = [prctile50medGH, median(medstat(:))];
%         % %         prctile75medGH = [prctile75medGH, prctile(medstat(:),75)];
%         %         prctile95medGHalt = [prctile95medGHalt, prctile(medstat(:),95)];
%         %
%         %         prctile05medGH = [prctile05medGH, prctile(ghs,05)];
%         %         prctile25medGH = [prctile25medGH, prctile(ghs,25)];
%         %         prctile50medGH = [prctile50medGH, median(ghs)];
%         %         prctile75medGH = [prctile75medGH, prctile(ghs,75)];
%         %         prctile95medGH = [prctile95medGH, prctile(ghs,95)];
%     end
% end
subplot(2,1,1);
hold on;
% scatter(gpowers(:), medGHhue(:),[],colors,'filled','marker','o');
p_s = p_srank;
scatter(gpowers(p_s>0.01), medGHhue(p_s>0.01),[],colors(p_s>0.01,:),'filled','marker','o');
scatter(gpowers(p_s<=0.01), medGHhue(p_s<=0.01),[],colors(p_s<=0.01,:),'marker','o');
hold on;
% scatter(powcenter(:), prctile50medGH(:),[],'k','filled','marker','d');
% scatter(powcenter(:), medianGH(:),[],'k','filled','marker','d');

% ydata = semedGH(:)*[-1,1] + medianGH(:);
ydata = semedGH(:)*[-1,1] + medGHhue(:);
% ydata = [prctile25medGH(:), prctile75medGH(:)];
% ydata = [prctile05medGH(:), prctile95medGH(:)];
% ydata = [prctile05medGHalt(:), prctile95medGHalt(:)];
for ind = 1:numel(gpowers(:))
    x = gpowers(ind)*[1,1]; y = ydata(ind,:);
%     line(x,y,'color', colors(ind,:));
    line(x,y,'color', 'k');
end
% title('M1')
% curax = subplot(2,1,2); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M2', 'fontweight','bold','fontsize',20,'edgecolor','none');
% ylim([1, max([ylim,3])]);
xlim([0, 25]);
ylim([1,3]);
xlabel('Electrode-averaged Gamma power (dB)')
ll = plot([0;max([21,xlim])],[2;2], 'color',[0.5,0.5,0.5],'linestyle',':');
ll = plot([0;max([21,xlim])],[2+0.05;2+0.05], 'color',[0.75,0.75,0.75],'linestyle',':');
ll = plot([0;max([21,xlim])],[2-0.05;2-0.05], 'color',[0.75,0.75,0.75],'linestyle',':');
ylabel('Peak-frequency ratios');
%% estimating ghratio and testing them for harmonic
alpamedianGHratio = median(allGHhue);
[se, ~] = getSEMedian(allGHhue);
alpasemedianGHratio = se;
alpaGHratio_SignTest = signtest(allGHhue(:), 2);
alpaGHratio_SignRank = signrank(allGHhue(:), 2);
% alpaGHratio_RankSum = ranksum(allGHhue(:), 2);

alpamedianGHratio_hues = median(medGHhue(:));
[se, ~] = getSEMedian(medGHhue(:));
alpasemedianGHratio_hues = se;
[sortedgpowers,sortedgpowerIndices] = sort(gpowers,'descend');
alpamedianGHratio_5Powerfulhues = median(medGHhue(sortedgpowerIndices(1:5)));
[se, ~] = getSEMedian(medGHhue(sortedgpowerIndices(1:5)));
alpasemedianGHratio_5Powerfulhues = se;
alpamedianGHratio_10Powerfulhues = median(medGHhue(sortedgpowerIndices(1:10)));
[se, ~] = getSEMedian(medGHhue(sortedgpowerIndices(1:10)));
alpasemedianGHratio_10Powerfulhues = se;
alpamedianGHratio_5Powerfulhues = median(medGHhue(sortedgpowerIndices(1:5)));
[se, ~] = getSEMedian(medGHhue(sortedgpowerIndices(1:5)));
alpasemedianGHratio_5Powerfulhues = se;

alpaGHratio_SignTest_hues = signtest(medGHhue(:), 2);
alpaGHratio_signRank_hues = signrank(medGHhue(:), 2);
alpaGHratio_signRank_10Powerfulhues = signrank(medGHhue(sortedgpowerIndices(1:10)), 2);
alpaGHratio_signRank_5Powerfulhues = signrank(medGHhue(sortedgpowerIndices(1:5)), 2);
% alpaGHratio_RankSum_hues = ranksum(medGHhue(:), 2);
%%
% M2
gpowers = tutuGpower_hue;
data = tutu_ghRatio;
% hueind = tutuhueind_elecALL;
medGHhue = [];
allGHhue = [];
hueindsGHhue = [];
semedGH = [];
p_sign = [];
p_srank = [];
for hue = 1:numel(tutu_ghRatio)
    %     ghs = data(hueind == hue);
    ghs = data{hue};
    medGHhue = [medGHhue, median(ghs)];
    allGHhue = [allGHhue; ghs(:)];
    hueindsGHhue = [hueindsGHhue; hue * ones(size(ghs(:)))];
    if numel(ghs) > 1
        [se, ~] = getSEMedian(ghs(:));
    else
        se = 0;
    end
%     if se == 0
%         disp('WELCOME')
%     end
    %%
    semedGH = [semedGH, se];
    pval = signtest(ghs(:),2);
    p_sign = [p_sign, pval];
    pval = signrank(ghs(:),2);
    p_srank = [p_srank, pval];
    
    %% testing achromatic distribution
% %     commented by raees
%     if hue == 37
%         oldf = gcf;
%         figure;
%         histogram(ghs(:),'binwidth',0.5);
%         title({'Distribution of frequency ratios for Achromatic gratings in M2',['p-value (Signed Rank test) =',num2str(p_srank(end))]});
%         sum(ghs(:)>2)
%         figure(oldf);
%     end
    %     [ci, medstat] = bootci(nboot, {bootfun, ghs(:)}, 'type','norm','alpha', 0.01);
    %     semedGH = [semedGH, std(medstat(:))];
    %     prctile05medGH = [prctile05medGH, prctile(medstat(:),05)];
    %     prctile25medGH = [prctile25medGH, prctile(medstat(:),25)];
    %     prctile50medGH = [prctile50medGH, median(medstat(:))];
    %     prctile75medGH = [prctile75medGH, prctile(medstat(:),75)];
    %     prctile95medGH = [prctile95medGH, prctile(medstat(:),95)];
end

% subplot(2,1,2);
% scatter(gpowers(:), medGH(:),[],colors,'marker','d');
% hold on;
% 
% medGH = [];
% medianGH = [];
% semedGH = [];
% prctile05medGHalt = [];
% prctile95medGHalt = [];
% prctile05medGH = [];
% prctile25medGH = [];
% prctile50medGH = [];
% prctile75medGH = [];
% prctile95medGH = [];
% powcenter = [];
% for pow = 0:binwidth:35
%     %     hueind = 1:36;
%     %     huesel = any(hueind(:) == find((tutuGpower_hue(:)' >= pow) & (tutuGpower_hue(:)' < pow+2.5)),2);
%     % %     ghs = data(huesel);
%     %     ghs = [];
%     %     hues = find(huesel);
%     %     for nh = hues(:)'
%     %         ghs = [ghs, data{nh}(:)'];
%     %     end
% %     hueindex = (gpowers(:) >= pow) & (gpowers(:) < pow+binwidth);
% %     ghs = medGHhue(hueindex);
%     hueindex = find((gpowers(:) >= pow) & (gpowers(:) < pow+binwidth));
%     ghs = allGHhue(any(hueindsGHhue(:) == hueindex(:)',2));
%     if ~isempty(ghs)
%         powcenter = [powcenter, pow+binwidth/2];
%         %         medGH = [medGH, median(ghs)];
%         %         [ci, medstat] = bootci(nboot, {bootfun, ghs(:)}, 'type','norm','alpha', 0.01);
%         %         semedGH = [semedGH, std(medstat(:))];
%         medianGH = [medianGH, median(ghs(:))];
%         if numel(ghs) > 1
%             [se, ~] = getSEMedian(ghs(:));
%         else
%             se = 0;
%         end
%         semedGH = [semedGH, se];
%         %         prctile05medGHalt = [prctile05medGHalt, prctile(medstat(:),05)];
%         % %         prctile25medGH = [prctile25medGH, prctile(medstat(:),25)];
%         % %         prctile50medGH = [prctile50medGH, median(medstat(:))];
%         % %         prctile75medGH = [prctile75medGH, prctile(medstat(:),75)];
%         %         prctile95medGHalt = [prctile95medGHalt, prctile(medstat(:),95)];
%         %         prctile05medGH = [prctile05medGH, prctile(ghs,05)];
%         %         prctile25medGH = [prctile25medGH, prctile(ghs,25)];
%         %         prctile50medGH = [prctile50medGH, median(ghs)];
%         %         prctile75medGH = [prctile75medGH, prctile(ghs,75)];
%         %         prctile95medGH = [prctile95medGH, prctile(ghs,95)];
%     end
% end
subplot(2,1,2);
hold on;
% scatter(gpowers(:), medGHhue(:),[],colors,'filled','marker','o');
p_s = p_srank;
scatter(gpowers(p_s>0.01), medGHhue(p_s>0.01),[],colors(p_s>0.01,:),'filled','marker','o');
scatter(gpowers(p_s<=0.01), medGHhue(p_s<=0.01),[],colors(p_s<=0.01,:),'marker','o');
hold on;
% scatter(powcenter(:), prctile50medGH(:),[],'k','filled','marker','d');
% scatter(powcenter(:), medianGH(:),[],'k','filled','marker','d');

% ydata = semedGH(:)*[-1,1] + medianGH(:);
ydata = semedGH(:)*[-1,1] + medGHhue(:);
% ydata = [prctile25medGH(:), prctile75medGH(:)];
% ydata = [prctile25medGH(:), prctile75medGH(:)];
% ydata = [prctile05medGHalt(:), prctile95medGHalt(:)];
for ind = 1:numel(gpowers)
    x = gpowers(ind)*[1,1]; y = ydata(ind,:);
%     line(x,y,'color', colors(ind,:));
    line(x,y,'color', 'k');
end
% title('M2')
% curax = subplot(2,1,2); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M2', 'fontweight','bold','fontsize',20,'edgecolor','none');
% ylim([1, max([ylim,3])]);
xlim([0, 25]);
ylim([1,3]);
xlabel('Electrode-averaged Gamma power (dB)')
ll = plot([0;max([21,xlim])],[2;2], 'color',[0.5,0.5,0.5],'linestyle',':');
ll = plot([0;max([21,xlim])],[2+0.05;2+0.05], 'color',[0.75,0.75,0.75],'linestyle',':');
ll = plot([0;max([21,xlim])],[2-0.05;2-0.05], 'color',[0.75,0.75,0.75],'linestyle',':');
ylabel('Peak-frequency ratios');
%% estimating ghratio and testing them for harmonic
tutumedianGHratio = median(allGHhue);
[se, ~] = getSEMedian(allGHhue);
tutusemedianGHratio = se;
tutuGHratio_SignTest = signtest(allGHhue(:), 2);
tutuGHratio_SignRank = signrank(allGHhue(:), 2);
% tutuGHratio_RankSum = ranksum(allGHhue(:), 2);

tutumedianGHratio_hues = median(medGHhue(:));
[se, ~] = getSEMedian(medGHhue(:));
tutusemedianGHratio_hues = se;
[sortedgpowers,sortedgpowerIndices] = sort(gpowers,'descend');
tutumedianGHratio_5Powerfulhues = median(medGHhue(sortedgpowerIndices(1:5)));
[se, ~] = getSEMedian(medGHhue(sortedgpowerIndices(1:5)));
tutusemedianGHratio_5Powerfulhues = se;
tutumedianGHratio_10Powerfulhues = median(medGHhue(sortedgpowerIndices(1:10)));
[se, ~] = getSEMedian(medGHhue(sortedgpowerIndices(1:10)));
tutusemedianGHratio_10Powerfulhues = se;
tutumedianGHratio_5Powerfulhues = median(medGHhue(sortedgpowerIndices(1:5)));
[se, ~] = getSEMedian(medGHhue(sortedgpowerIndices(1:5)));
tutusemedianGHratio_5Powerfulhues = se;
tutuGHratio_SignTest_hues = signtest(medGHhue(:), 2);
tutuGHratio_signRank_hues = signrank(medGHhue(:), 2);
tutuGHratio_signRank_10Powerfulhues = signrank(medGHhue(sortedgpowerIndices(1:10)), 2);
tutuGHratio_signRank_5Powerfulhues = signrank(medGHhue(sortedgpowerIndices(1:5)), 2);
% tutuGHratio_RankSum_hues = ranksum(medGHhue(:), 2);
%%
curax = subplot(2,1,1); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M1', 'fontweight','bold','fontsize',20,'edgecolor','none');
curax = subplot(2,1,2); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M2', 'fontweight','bold','fontsize',20,'edgecolor','none');
%%
% postformat
sgtitle('Peak-frequency ratios across stimuli','fontweight','bold','fontsize',20);
f = f3; f.WindowState = 'maximized'; postformatFig;

%%
% f = f4_1;
% postformatFig;
f = gcf;
savefig(f, [destination, 'Fig3_new_2SE_Draft_8.fig']);
print([destination, 'Fig3_new_2SE_Draft_8.png'],'-dpng','-r480');
mkdir([destination, '/Draft 9'])
print([destination, '/Draft 9/Fig2.tif'],'-dtiff','-r480');
disp('printed');

%%
elecs37 = alpahueind_elecALL(:)==1:37; elecs37 = sum(elecs37,1)/64;
disp(['M1 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
elecs37 = tutuhueind_elecALL(:)==1:37; elecs37 = sum(elecs37,1)/16;
disp(['M2 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
% %%
% elecs36 = alpahueind_elecALL(:)==1:36; elecs36 = sum(elecs36,1)/64;
% disp(['M1 Electrodes per stimulus: ',num2str(mean(elecs36)),'+-',num2str(std(elecs36))]);
% elecs36 = tutuhueind_elecALL(:)==1:36; elecs36 = sum(elecs36,1)/16;
% disp(['M2 Electrodes per stimulus: ',num2str(mean(elecs36)),'+-',num2str(std(elecs36))]);
%%
elecs37 = alpahueind_elecALL(:)==37; elecs37 = sum(elecs37,1)/64;
disp(['M1 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
elecs37 = tutuhueind_elecALL(:)==37; elecs37 = sum(elecs37,1)/16;
disp(['M2 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
%%
disp(['M1 median of median-FreqRatio/hue: ',num2str(alpamedianGHratio_hues),'+/-', num2str(alpasemedianGHratio_hues),'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 median of median-FreqRatio/hue: ',num2str(tutumedianGHratio_hues),'+/-', num2str(tutusemedianGHratio_hues),'; p=',num2str(tutuGHratio_signRank_hues)])
% % 
% disp(['M1 median of median-FreqRatio/top5hue: ',num2str(alpamedianGHratio_5Powerfulhues),'+/-', num2str(alpasemedianGHratio_5Powerfulhues)])
% disp(['M2 median of median-FreqRatio/top5hue: ',num2str(tutumedianGHratio_5Powerfulhues),'+/-', num2str(tutusemedianGHratio_5Powerfulhues)])

disp(['M1 median of median-FreqRatio/top10hue: ',num2str(alpamedianGHratio_10Powerfulhues),'+/-', num2str(alpasemedianGHratio_10Powerfulhues),'; p=',num2str(alpaGHratio_signRank_10Powerfulhues)])
disp(['M2 median of median-FreqRatio/top10hue: ',num2str(tutumedianGHratio_10Powerfulhues),'+/-', num2str(tutusemedianGHratio_10Powerfulhues),'; p=',num2str(tutuGHratio_signRank_10Powerfulhues)])

disp(['M1 median of median-FreqRatio/top5hue: ',num2str(alpamedianGHratio_5Powerfulhues),'+/-', num2str(alpasemedianGHratio_5Powerfulhues),'; p=',num2str(alpaGHratio_signRank_5Powerfulhues)])
disp(['M2 median of median-FreqRatio/top5hue: ',num2str(tutumedianGHratio_5Powerfulhues),'+/-', num2str(tutusemedianGHratio_5Powerfulhues),'; p=',num2str(tutuGHratio_signRank_5Powerfulhues)])
