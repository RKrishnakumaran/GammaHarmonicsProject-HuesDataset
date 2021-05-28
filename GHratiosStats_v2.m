% load('E:\RK Gamma Harmonics paper\Paper draft #1\Figures - draft #1\Raees data and figure\pedestalFitGHratioAll.mat');
load('.\Data\alpaGH.mat', 'ratioGH');
alpa_ghRatio = ratioGH;
load('.\Data\tutuGH.mat', 'ratioGH');
tutu_ghRatio = ratioGH;
% load('E:\RK Gamma Harmonics paper\Paper draft #1\Figures - draft #1\Raees data and figure\forSelectedHues.mat');

%%
GHratios = [];
alpantrials = [];
if ~exist('alpaColorCases')
   alpaColorCases = transpose(1:36); 
end
for h = alpaColorCases'
    temp = alpa_ghRatio{h};
    GHratios = [GHratios; temp(:)];
    alpantrials = [alpantrials, size(alpa_ghRatio{h},2)];
end
alpaselGHratios = GHratios;
GHratios = [];
tutuntrials = [];

if ~exist('tutuColorCases')
   tutuColorCases = transpose(1:36); 
end
for h = tutuColorCases'
    temp =  tutu_ghRatio{h};
    GHratios = [GHratios; temp(:)];
    tutuntrials = [tutuntrials, size(tutu_ghRatio{h},2)];
end
tutuselGHratios = GHratios;
%% statistics of GH frequency ratio
bootfun = @(x) median(x,'all');
[aci, astat] = bootci(10000, {bootfun, alpaselGHratios(:)}, 'type','norm','alpha', 0.05);
[tci, tstat] = bootci(10000, {bootfun, tutuselGHratios(:)}, 'type','norm','alpha',0.05);
alpaMedianGHratio = median(alpaselGHratios(:))
alpaMedianGHratioSE = std(astat)
tutuMedianGHratio = median(tutuselGHratios(:))
tutuMedianGHratioSE = std(tstat)
%%
figure; histogram(alpaselGHratios(:))
figure; histogram(tutuselGHratios(:))
[tselp, tselh, tsselignstat] = signrank(tutuselGHratios(:),2); tselp, tselh
[aselp, aselh, asselignstat] = signrank(alpaselGHratios(:),2); aselp, aselh
%% statistics of population (data set analysed)
alpaMeanTrials = mean(alpantrials)
alpaSDTrials = std(alpantrials)
alpaMINTrials = min(alpantrials)
tutuMeanTrials = mean(tutuntrials)
tutuSDTrials = std(tutuntrials)
tutuMINTrials = min(tutuntrials)


%% Test Results
verdicts = {'Second Peak is harmonic!','Not Harmonic! Need to check pedestal!'}
disp('### GHRatio Test Results ###')
disp('M1:')
disp(['GHratio = ',num2str(alpaMedianGHratio),' +/- ',num2str(alpaMedianGHratioSE)]);
disp(['verdict: ', verdicts{1*(aselp<=0.01)+1}]);
disp('M2:')
disp(['GHratio = ',num2str(tutuMedianGHratio),' +/- ',num2str(tutuMedianGHratioSE)]);
disp(['verdict: ', verdicts{1*(tselp<=0.01)+1}]);