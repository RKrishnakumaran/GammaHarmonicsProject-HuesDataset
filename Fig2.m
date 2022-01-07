%% New Fig 3: GHratio vs harmonic and gamma peak powers (PDstats_Power.m)
destination = './Figures/';
mkdir(destination);
%% load data and/or source figures

PDstats;
%% Fig 3 median and stderr of GH ratios vs Gamma power - each hue
f3 = figure;
colors = [hsv(36);0.65*[1,1,1]];
nboot = 1000; bootfun = @(x)median(x);
binwidth = 2;

% M1
gpowers = alpaGpower_hue;
data = alpa_ghRatio;
medGHhue = [];
allGHhue = [];
hueindsGHhue = [];
semedGH = [];
p_sign = [];
p_srank = [];
for hue = 1:numel(alpa_ghRatio)
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
end

subplot(2,1,1);
hold on;
p_s = p_srank;
scatter(gpowers(p_s>0.01), medGHhue(p_s>0.01),[],colors(p_s>0.01,:),'filled','marker','o');
scatter(gpowers(p_s<=0.01), medGHhue(p_s<=0.01),[],colors(p_s<=0.01,:),'marker','o');
hold on;
ydata = semedGH(:)*[-1,1] + medGHhue(:);
for ind = 1:numel(gpowers(:))
    x = gpowers(ind)*[1,1]; y = ydata(ind,:);
    line(x,y,'color', 'k');
end
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
%%
% M2
gpowers = tutuGpower_hue;
data = tutu_ghRatio;
medGHhue = [];
allGHhue = [];
hueindsGHhue = [];
semedGH = [];
p_sign = [];
p_srank = [];
for hue = 1:numel(tutu_ghRatio)
    ghs = data{hue};
    medGHhue = [medGHhue, median(ghs)];
    allGHhue = [allGHhue; ghs(:)];
    hueindsGHhue = [hueindsGHhue; hue * ones(size(ghs(:)))];
    if numel(ghs) > 1
        [se, ~] = getSEMedian(ghs(:));
    else
        se = 0;
    end
    %%
    semedGH = [semedGH, se];
    pval = signtest(ghs(:),2);
    p_sign = [p_sign, pval];
    pval = signrank(ghs(:),2);
    p_srank = [p_srank, pval];
    
end

subplot(2,1,2);
hold on;
p_s = p_srank;
scatter(gpowers(p_s>0.01), medGHhue(p_s>0.01),[],colors(p_s>0.01,:),'filled','marker','o');
scatter(gpowers(p_s<=0.01), medGHhue(p_s<=0.01),[],colors(p_s<=0.01,:),'marker','o');
hold on;
ydata = semedGH(:)*[-1,1] + medGHhue(:);
for ind = 1:numel(gpowers)
    x = gpowers(ind)*[1,1]; y = ydata(ind,:);
    line(x,y,'color', 'k');
end
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
%%
curax = subplot(2,1,1); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M1', 'fontweight','bold','fontsize',20,'edgecolor','none');
curax = subplot(2,1,2); annotation('textbox', [curax.OuterPosition([1,2])+[0.0125, 0.42].*curax.OuterPosition([3,4]), 0.1, 0.1],'string','M2', 'fontweight','bold','fontsize',20,'edgecolor','none');
%%
sgtitle('Peak-frequency ratios across stimuli','fontweight','bold','fontsize',20);
f = f3; f.WindowState = 'maximized'; postformatFig;

%%
f = gcf;
savefig(f, [destination, 'Fig2.fig']);
print([destination, 'Fig2.png'],'-dpng','-r480');
print([destination, 'Fig2.tif'],'-dtiff','-r480');
disp('printed');

%%
elecs37 = alpahueind_elecALL(:)==1:37; elecs37 = sum(elecs37,1)/64;
disp(['M1 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
elecs37 = tutuhueind_elecALL(:)==1:37; elecs37 = sum(elecs37,1)/16;
disp(['M2 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
%%
elecs37 = alpahueind_elecALL(:)==37; elecs37 = sum(elecs37,1)/64;
disp(['M1 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
elecs37 = tutuhueind_elecALL(:)==37; elecs37 = sum(elecs37,1)/16;
disp(['M2 Electrodes per stimulus: ',num2str(mean(elecs37)),'+-',num2str(std(elecs37))]);
%%
disp(['M1 median of median-FreqRatio/hue: ',num2str(alpamedianGHratio_hues),'+/-', num2str(alpasemedianGHratio_hues),'; p=',num2str(alpaGHratio_signRank_hues)])
disp(['M2 median of median-FreqRatio/hue: ',num2str(tutumedianGHratio_hues),'+/-', num2str(tutusemedianGHratio_hues),'; p=',num2str(tutuGHratio_signRank_hues)])

disp(['M1 median of median-FreqRatio/top10hue: ',num2str(alpamedianGHratio_10Powerfulhues),'+/-', num2str(alpasemedianGHratio_10Powerfulhues),'; p=',num2str(alpaGHratio_signRank_10Powerfulhues)])
disp(['M2 median of median-FreqRatio/top10hue: ',num2str(tutumedianGHratio_10Powerfulhues),'+/-', num2str(tutusemedianGHratio_10Powerfulhues),'; p=',num2str(tutuGHratio_signRank_10Powerfulhues)])

disp(['M1 median of median-FreqRatio/top5hue: ',num2str(alpamedianGHratio_5Powerfulhues),'+/-', num2str(alpasemedianGHratio_5Powerfulhues),'; p=',num2str(alpaGHratio_signRank_5Powerfulhues)])
disp(['M2 median of median-FreqRatio/top5hue: ',num2str(tutumedianGHratio_5Powerfulhues),'+/-', num2str(tutusemedianGHratio_5Powerfulhues),'; p=',num2str(tutuGHratio_signRank_5Powerfulhues)])
