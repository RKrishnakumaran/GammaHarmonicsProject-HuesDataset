% load('E:\RK Gamma Harmonics paper\Paper draft #1\Figures - draft #1\Raees data and figure\pedestalFitGHratioAll.mat');
load('.\Data\alpaGH.mat', 'ratioGH','powerGamma');
alpa_ghRatio = ratioGH;
alpa_Gpower = powerGamma;
load('.\Data\tutuGH.mat', 'ratioGH','powerGamma');
tutu_ghRatio = ratioGH;
tutu_Gpower = powerGamma;
% load('E:\RK Gamma Harmonics paper\Paper draft #1\Figures - draft #1\Raees data and figure\forSelectedHues.mat');
%%
bootfun = @(x) median(x,'all');

GHratios = [];
alpantrials = [];
% if ~exist('alpaColorCases','var')
%    alpaColorCases = transpose(1:36); 
% end
%%
alpaColorCases = [];
alpaColorCasesGpower = [];
for h = 1:numel(alpa_Gpower)
%     if any(alpa_Gpower{h} >= 14.5,'all')
%     if mean(log(mean(exp(alpa_Gpower{h}),1)),'all') >= 14
    alpaColorCases = [alpaColorCases;h];
    alpaColorCasesGpower = [alpaColorCasesGpower;mean(log(mean(exp(alpa_Gpower{h}),1)),'all')];
%     end
%     end
end
%%
alpaGHratiomedian_h = [];
alpaGHratioiqr_h = [];
alpaGHratiomedianSE_h = [];
alpap_h = [];
for h = alpaColorCases'
    temp = alpa_ghRatio{h};
%     tempsel = alpa_Gpower{h} >= 14.5;
    temp = temp(:); %tempsel = tempsel(:); temp = temp(tempsel);
    [aci, astat] = bootci(250, {bootfun, temp(:)}, 'type','norm','alpha', 0.05);
    alpaGHratiomedian_h = [alpaGHratiomedian_h, median(temp(:))];
    alpaGHratioiqr_h = [alpaGHratioiqr_h, iqr(temp(:))];
    alpaGHratiomedianSE_h = [alpaGHratiomedianSE_h, std(astat)];
    [aselp, aselh, asselignstat] = signrank(temp(:),2); aselp, aselh
    alpap_h = [alpap_h, aselp];
    GHratios = [GHratios; temp(:)];
    alpantrials = [alpantrials, size(alpa_ghRatio{h},2)];
end
alpaselGHratios = GHratios;
GHratios = [];
tutuntrials = [];

% if ~exist('tutuColorCases','var')
%    tutuColorCases = transpose(1:36); 
% end
%%
tutuColorCases = []; 
tutuColorCases_Gpower = [];
for h = 1:numel(tutu_Gpower)
%     if any(tutu_Gpower{h} >= 18,'all')
%     if mean(log(mean(exp(tutu_Gpower{h}),1)),'all') >= 14
    tutuColorCases = [tutuColorCases;h];
    tutuColorCases_Gpower = [tutuColorCases_Gpower;mean(log(mean(exp(tutu_Gpower{h}),1)),'all')];
%     end
%     end
end
%%
tutuGHratiomedian_h = [];
tutuGHratioiqr_h = [];
tutuGHratiomedianSE_h = [];
tutup_h = [];
GHratios = [];
for h = tutuColorCases'
    temp =  tutu_ghRatio{h};
%     tempsel = tutu_Gpower{h} >= 18;
    temp = temp(:); %tempsel = tempsel(:); temp = temp(tempsel);
    [tci, tstat] = bootci(250, {bootfun, temp(:)}, 'type','norm','alpha', 0.05);
    tutuGHratiomedian_h = [tutuGHratiomedian_h, median(temp(:))];
    tutuGHratioiqr_h = [tutuGHratioiqr_h, iqr(temp(:))];
    tutuGHratiomedianSE_h = [tutuGHratiomedianSE_h, std(tstat)];
    [tselp, tselh, tsselignstat] = signrank(temp(:),2);
    tutup_h = [tutup_h, tselp];
    GHratios = [GHratios; temp(:)];
    tutuntrials = [tutuntrials, size(tutu_ghRatio{h},2)];
end
tutuselGHratios = GHratios;
%% statistics of GH frequency ratio
[aci, astat] = bootci(10000, {bootfun, alpaselGHratios(:)}, 'type','norm','alpha', 0.05);
[tci, tstat] = bootci(10000, {bootfun, tutuselGHratios(:)}, 'type','norm','alpha',0.05);
alpaMedianGHratio = median(alpaselGHratios(:))
alpaMedianGHratioSE = std(astat)
tutuMedianGHratio = median(tutuselGHratios(:))
tutuMedianGHratioSE = std(tstat)
%%
% figure; histogram(alpaselGHratios(:))
% figure; histogram(tutuselGHratios(:))
[tselp, tselh, tsselignstat] = signrank(tutuselGHratios(:),2); tselp, tselh
[aselp, aselh, asselignstat] = signrank(alpaselGHratios(:),2); aselp, aselh
[tselp2, tselh2, tsselignstat2] = signtest(tutuselGHratios(:),2); tselp2, tselh2
[aselp2, aselh2, asselignstat2] = signtest(alpaselGHratios(:),2); aselp2, aselh2
%% statistics of population (data set analysed)
alpaMeanTrials = mean(alpantrials)
alpaSDTrials = std(alpantrials)
alpaMINTrials = min(alpantrials)
tutuMeanTrials = mean(tutuntrials)
tutuSDTrials = std(tutuntrials)
tutuMINTrials = min(tutuntrials)


%% Test Results
verdicts = {'Second Peak is harmonic!','Not Harmonic! Need to check pedestal!'};
disp('### GHRatio Sign Rank Test Results ###')
disp('M1:')
disp(['GHratio = ',num2str(alpaMedianGHratio),' +/- ',num2str(alpaMedianGHratioSE)]);
disp(['verdict: ', verdicts{1*(aselp<=0.01)+1}]);
disp('M2:')
disp(['GHratio = ',num2str(tutuMedianGHratio),' +/- ',num2str(tutuMedianGHratioSE)]);
disp(['verdict: ', verdicts{1*(tselp<=0.01)+1}]);
%% Test Results
verdicts = {'Second Peak is harmonic!','Not Harmonic! Need to check pedestal!'};
disp('### GHRatio Sign Test Results ###')
disp('M1:')
disp(['GHratio = ',num2str(alpaMedianGHratio),' +/- ',num2str(alpaMedianGHratioSE)]);
disp(['verdict: ', verdicts{1*(aselp2<=0.01)+1}]);
disp('M2:')
disp(['GHratio = ',num2str(tutuMedianGHratio),' +/- ',num2str(tutuMedianGHratioSE)]);
disp(['verdict: ', verdicts{1*(tselp2<=0.01)+1}]);