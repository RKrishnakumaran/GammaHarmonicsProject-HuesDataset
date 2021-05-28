% load('.\Data\PD_all.mat');
% load('.\Data\alpaPD.mat');
load('./Data/alpaPD.mat');
alpaPD = PD;
% load('.\Data\tutuPD.mat');
load('./Data/tutuPD.mat');
tutuPD = PD;
clear PD
%%
% load('E:\RK Gamma Harmonics paper\Paper draft #1\Figures - draft #1\Raees data and figure\forSelectedHues.mat');
load('./Data/alpaGH.mat','powerGamma','ratioGH');
alpa_ghRatio = ratioGH; alpaGpower = powerGamma;
load('./Data/tutuGH.mat','powerGamma','ratioGH');
tutu_ghRatio = ratioGH; tutuGpower = powerGamma;
%load('E:\RK Gamma Harmonics paper\Paper draft #1\Figures - draft #1\Raees data and figure\pedestalFitGHpowerAll.mat');
clear ratioGH; clear powerGamma;
% %% comment out if Gpowers are proper, i.e., already multiplied by 10
% for hue = 1:numel(alpaGpower)
%     alpaGpower{hue} = 10 * alpaGpower{hue};
% end
% for hue = 1:numel(tutuGpower)
%     tutuGpower{hue} = 10 * tutuGpower{hue};
% end
%%
% destinationFolder = '.\Figures_codes\tempFigs';
% mkdir(destinationFolder);
%%
% Mean a)gamma and b)harmonic power during each hue, across all trials and electrodes
% alpaGpower = []; alpaGpower = []; 
% tutuGpower = []; tutuGpower = [];
elecDim = 2; TrialDim = 1;
%% Sort hues by a)gamma power and b)harmonic power
alpaGpower_hue = []; % alpaHpower_hue = [];
% alpaGpower_hueSE = []; % alpaHpower_hueSE = [];
% alpaGH_hueMed = [];
% alpaGH_hueIQR = [];
% alpaGH_hueALL = []; 
alpahueind_elecALL = [];
for nhue = 1:numel(alpaGpower)
    alpaGpower_selecting_hues = mean(alpaGpower{nhue});%mean(log(mean(exp(alpaGpower{nhue}),TrialDim)),elecDim);
    alpaGpower_hue = [alpaGpower_hue, alpaGpower_selecting_hues(:)];
    
    phdiff = squeeze(alpaPD{nhue}(1,:,:));
    gh = alpa_ghRatio{nhue};
%     alpaGH_hueALL = [alpaGH_hueALL; gh(:)];
    ind = nhue*ones(size(phdiff));
    alpahueind_elecALL = [alpahueind_elecALL; ind(:)];
    
%     gh = median(alpa_ghRatio{nhue},'all');
%     alpaGH_hueMed = [alpaGH_hueMed, gh(:)];
%     gh = iqr(alpa_ghRatio{nhue},'all');
%     alpaGH_hueIQR = [alpaGH_hueIQR, gh(:)];
%     alpaGpower_hue = [alpaGpower_hue, alpaGpower_selecting_hues{nhue}];
%     alpaGpower_hueSE = [alpaGpower_hueSE, std(alpaGpower{nhue}, 0, 'all')];
%     alpaHpower_hueSE = [alpaHpower_hueSE, std(alpaHpower{nhue}, 0, 'all')];
%     alpaGpower_hue = [alpaGpower_hue, mean(alpaGpower{nhue}, 'all')];
%     alpaHpower_hue = [alpaHpower_hue, mean(alpaHpower{nhue}, 'all')];
%     alpaGpower_hueSE = [alpaGpower_hueSE, std(alpaGpower{nhue}, 0, 'all')];
%     alpaHpower_hueSE = [alpaHpower_hueSE, std(alpaHpower{nhue}, 0, 'all')];
end
% [alpahueGpower_sorted, alpahueindex_Gsorted] = sort(alpaGpower_hue);
[alpahueGpower_sorted, alpahueindex_Gsorted] = sort(alpaGpower_hue,'descend');
hues = hsv(numel(alpa_ghRatio));
% alpahueGsorted = hues(alpahueindex_Gsorted, :); %(alpahueindex_Hsorted)*360/numel(alpa_ghRatio);
alpahueGsorted = hues(alpahueindex_Gsorted, :); %(alpahueindex_Gsorted)*360/numel(alpa_ghRatio);
alpaGpower_hueALL = alpaGpower_hue(alpahueind_elecALL);
[alpahueGpower_ALLsorted, alpahueindex_GALLsorted] = sort(alpaGpower_hueALL(:),'descend');

tutuGpower_hue = []; % tutuGpower_hue = [];
% tutuGpower_hueSE = []; % tutuGpower_hueSE = [];
% tutuGH_hueMed = [];
% tutuGH_hueIQR = [];
% tutuGH_hueALL = []; 
tutuhueind_elecALL = [];
for nhue = 1:numel(tutuGpower)
    tutuGpower_selecting_hues = mean(tutuGpower{nhue});% mean(log(mean(exp(tutuGpower{nhue}),TrialDim)),elecDim);
    tutuGpower_hue = [tutuGpower_hue, tutuGpower_selecting_hues];
    
    phdiff = squeeze(tutuPD{nhue}(1,:,:));
    gh = tutu_ghRatio{nhue};
%     tutuGH_hueALL = [tutuGH_hueALL; gh(:)];
    ind = nhue*ones(size(phdiff));
    tutuhueind_elecALL = [tutuhueind_elecALL; ind(:)];
%     gh = median(tutu_ghRatio{nhue},'all');
%     tutuGH_hueMed = [tutuGH_hueMed, gh(:)];
%     gh = iqr(tutu_ghRatio{nhue},'all');
%     tutuGH_hueIQR = [tutuGH_hueIQR, gh(:)];
%     tutuGpower_hue = [tutuGpower_hue, tutuGpower_selecting_hues{nhue}];
%     tutuGpower_hueSE = [tutuGpower_hueSE, std(tutuGpower, 0, 'all')];
%     tutuHpower_hueSE = [tutuHpower_hueSE, std(tutuHpower, 0, 'all')];
%     tutuGpower_hue = [tutuGpower_hue, mean(tutuGpower, 'all')];
%     tutuHpower_hue = [tutuHpower_hue, mean(tutuHpower, 'all')];
%     tutuGpower_hueSE = [tutuGpower_hueSE, std(tutuGpower, 0, 'all')];
%     tutuHpower_hueSE = [tutuHpower_hueSE, std(tutuHpower, 0, 'all')];
end
% [tutuhueGpower_sorted, tutuhueindex_Gsorted] = sort(tutuGpower_hue);
[tutuhueGpower_sorted, tutuhueindex_Gsorted] = sort(tutuGpower_hue,'descend');
hues = hsv(numel(tutu_ghRatio));
% tutuhueGsorted = hues(tutuhueindex_Gsorted, :); %(tutuhueindex_Hsorted)*360/numel(tutu_ghRatio);
tutuhueGsorted = hues(tutuhueindex_Gsorted, :); %(tutuhueindex_Gsorted)*360/numel(tutu_ghRatio);
tutuGpower_hueALL = tutuGpower_hue(tutuhueind_elecALL);
[tutuhueGpower_ALLsorted, tutuhueindex_GALLsorted] = sort(tutuGpower_hueALL(:),'descend');

%% Hues (pairs- sorted by elec-avg,trial-avg gamma power) vs indiv-trial GH ratio
% figure; subplot(2,1,1); boxplot(alpaGH_hueALL, floor(alpaGpower_hueALL(:)/0.25)*0.25); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuGH_hueALL, floor(tutuGpower_hueALL(:)/0.25)*0.25); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('GHratios vs Hue gamma power');
%%
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_hueALL(:)/0.25)*0.25,alpaGH_hueALL,[],colors(alpahueind_elecALL,:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_hueALL(:)/0.25)*0.25,tutuGH_hueALL,[],colors(tutuhueind_elecALL,:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('GHratios vs Hue gamma power');

%% Sort electrodes by a)gamma power and b)harmonic power
alpaGpower_elecMean = [];
% alpaGpower_elecSE = [];
alpaGpower_elecMed = []; %alpaHpower_elec = [];
alpaGpower_elecIQR = []; %alpaHpower_elecSE = [];
alpaGH_elecMed = [];
alpaGH_elecIQR = [];
alpaGH_elecALL = []; alpaelecind_elecALL = []; alpaGpower_elecALL = [];
for nhue = 1:numel(alpaGpower)
%     pwr = mean(alpaGpower, TrialDim); 
%     pwr = median(alpaGpower{nhue},TrialDim);
%     alpaGpower_elecMed = [alpaGpower_elecMed, pwr(:)];
%     pwr = iqr(alpaGpower{nhue},TrialDim);
%     alpaGpower_elecIQR = [alpaGpower_elecIQR, pwr(:)];
    pwr = alpaGpower{nhue}; %log(mean(exp(alpaGpower{nhue}),TrialDim));
    alpaGpower_elecMean = [alpaGpower_elecMean, pwr(:)];
    
    
    phdiff = squeeze(alpaPD{nhue}(1,:,:));
    gh = alpa_ghRatio{nhue};
    alpaGH_elecALL = [alpaGH_elecALL; ones([size(phdiff,1),1])*gh(:)'];
    ind = ones([size(phdiff,1),1])*(1:size(phdiff,2));
    alpaelecind_elecALL = [alpaelecind_elecALL; ind(:)];
%     gh = median(alpa_ghRatio{nhue},TrialDim);
%     alpaGH_elecMed = [alpaGH_elecMed, gh(:)];
%     gh = iqr(alpa_ghRatio{nhue},TrialDim);
%     alpaGH_elecIQR = [alpaGH_elecIQR, gh(:)];
%     pwr = log(std(exp(alpaGpower{nhue}),TrialDim));
%     alpaGpower_elecSE = [alpaGpower_elecSE, pwr(:)];
    
%     pwr = mean(alpaHpower, TrialDim); 
%     pwr = alpaHpower_selecting_hues{nhue};
%     alpaHpower_elec = [alpaHpower_elec, pwr(:)];
%     pwr = std(alpaGpower, 0, TrialDim); alpaGpower_elecSE = [alpaGpower_elecSE, pwr(:)];
%     pwr = std(alpaHpower, 0, TrialDim); alpaHpower_elecSE = [alpaHpower_elecSE, pwr(:)];
end
% alpaGpower_elecMed = alpaGpower_elecMed;
% alpaGpower_elecIQR = alpaGpower_elecIQR;
% alpaGpower_elecMean = alpaGpower_elecMean;
% alpaGH_elecMed = alpaGH_elecMed;
% alpaGH_elecIQR = alpaGH_elecIQR;
% tutuGpower_elecSE = [];
% alpaGpower_elec = mean(alpaGpower_elec, 2);
% alpaGpower_elec = mean(alpaGpower_elec, 2);
% alpaGpower_elecSE = std(alpaGpower_elec, 0, 2);
% alpaGpower_elecSE = std(alpaGpower_elec, 0, 2);
% [alpaelecGpower_sorted, alpaelecindex_Gsorted] = sort(alpaGpower_elec);
[alpaelecGpower_sorted, alpaelecindex_Gsorted] = sort(alpaGpower_elecMean(:),'descend');
alpaGpower_elecMean = alpaGpower_elecMean(:);
alpaGpower_elecALL = alpaGpower_elecMean(alpaelecind_elecALL(:)+(alpahueind_elecALL(:)-1)*numel(alpaGpower_elecMean)/numel(alpaGpower));
[alpaelecGpower_ALLsorted, alpaelecindex_GALLsorted] = sort(alpaGpower_elecALL(:),'descend');
% hues = hsv(numel(alpa_ghRatio));
% alpahueHsorted = hues(alpaelecindex_Hsorted, :); %(alpahueindex_Hsorted)*360/numel(alpa_ghRatio);
% alpahueGsorted = hues(alpaelecindex_Gsorted, :); %(alpahueindex_Gsorted)*360/numel(alpa_ghRatio);
tutuGpower_elecMean = [];
% tutuGpower_elecSE = [];
% tutuGpower_elecMed = []; %tutuHpower_elec = [];
% tutuGpower_elecIQR = []; %tutuHpower_elecSE = [];
% tutuGH_elecMed = [];
% tutuGH_elecIQR = [];
tutuGH_elecALL = []; tutuelecind_elecALL = []; tutuGpower_elecALL = [];
for nhue = 1:numel(tutuGpower)
%     pwr = mean(tutuGpower, TrialDim); 
%     pwr = median(tutuGpower{nhue},TrialDim);
%     tutuGpower_elecMed = [tutuGpower_elecMed, pwr(:)];
%     pwr = iqr(tutuGpower{nhue},TrialDim);
%     tutuGpower_elecIQR = [tutuGpower_elecIQR, pwr(:)];
    pwr = tutuGpower{nhue}; %log(mean(exp(tutuGpower{nhue}),TrialDim));
    tutuGpower_elecMean = [tutuGpower_elecMean, pwr(:)];
    
    phdiff = squeeze(tutuPD{nhue}(1,:,:));
    gh = tutu_ghRatio{nhue};
    tutuGH_elecALL = [tutuGH_elecALL; ones([size(phdiff,1),1])*gh(:)'];
    ind = ones([size(phdiff,1),1])*(1:size(phdiff,2));
    tutuelecind_elecALL = [tutuelecind_elecALL; ind(:)];
%     gh = median(tutu_ghRatio{nhue},TrialDim);
%     tutuGH_elecMed = [tutuGH_elecMed, gh(:)];
%     gh = iqr(tutu_ghRatio{nhue},TrialDim);
%     tutuGH_elecIQR = [tutuGH_elecIQR, gh(:)];
%     pwr = log(std(exp(tutuGpower{nhue}),TrialDim));
%     tutuGpower_elecSE = [tutuGpower_elecSE, pwr(:)];
    
%     pwr = mean(tutuHpower, TrialDim); 
%     pwr = tutuHpower_selecting_hues{nhue};
%     tutuHpower_elec = [tutuHpower_elec, pwr(:)];
%     pwr = std(tutuGpower, 0, TrialDim); tutuGpower_elecSE = [tutuGpower_elecSE, pwr(:)];
%     pwr = std(tutuHpower, 0, TrialDim); tutuHpower_elecSE = [tutuHpower_elecSE, pwr(:)];
end
% tutuGpower_elecMed = tutuGpower_elecMed;
% tutuGpower_elecIQR = tutuGpower_elecIQR;
% tutuGpower_elecMean = tutuGpower_elecMean;
% tutuGH_elecMed = tutuGH_elecMed;
% tutuGH_elecIQR = tutuGH_elecIQR;
% tutuGpower_elecSE = tutuGpower_elecSE(:);
% tutuGpower_elec = mean(tutuGpower_elec, 2);
% tutuGpower_elec = mean(tutuGpower_elec, 2);
% tutuGpower_elecSE = std(tutuGpower_elec, 0, 2);
% tutuGpower_elecSE = std(tutuGpower_elec, 0, 2);
% [tutuelecGpower_sorted, tutuelecindex_Gsorted] = sort(tutuGpower_elec);
[tutuelecGpower_sorted, tutuelecindex_Gsorted] = sort(tutuGpower_elecMean(:),'descend');
tutuGpower_elecMean = tutuGpower_elecMean(:);
tutuGpower_elecALL = tutuGpower_elecMean(tutuelecind_elecALL(:)+(tutuhueind_elecALL(:)-1)*numel(tutuGpower_elecMean)/numel(tutuGpower));
[tutuelecGpower_ALLsorted, tutuelecindex_GALLsorted] = sort(tutuGpower_elecALL(:),'descend');
% hues = hsv(numel(tutu_ghRatio));
% tutuhueHsorted = hues(tutuhueindex_Hsorted, :); %(tutuhueindex_Hsorted)*360/numel(tutu_ghRatio);
% tutuhueGsorted = hues(tutuhueindex_Gsorted, :); %(tutuhueindex_Gsorted)*360/numel(tutu_ghRatio);

%% Elec,Hues (pairs- sorted by trial-avg gamma power) vs indiv-trial GH ratio
% figure; subplot(2,1,1); boxplot(alpaGH_elecALL(alpaelecindex_GALLsorted), floor(alpaelecGpower_ALLsorted/0.5)*0.5); grid on; title('M1');
% subplot(2,1,2); boxplot(tutuGH_elecALL(tutuelecindex_GALLsorted), floor(tutuelecGpower_ALLsorted/0.5)*0.5); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('GHratios');
% figure; subplot(2,1,1); boxplot(alpaGH_elecALL, floor(alpaGpower_elecALL(:)/0.05)*0.05); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuGH_elecALL, floor(tutuGpower_elecALL(:)/0.05)*0.05); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('GHratios');
%%
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_elecALL(:)/0.5)*0.5,alpaGH_elecALL, [],colors(alpahueind_elecALL(:),:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_elecALL(:)/0.5)*0.5,tutuGH_elecALL, [],colors(tutuhueind_elecALL(:),:), '.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('GHratios');
%% PDs
alpaPDs_hueALL = [];
tutuPDs_hueALL = [];
for nhue = 1:numel(tutuGpower)
    pd = squeeze(circ_mean(alpaPD{nhue}));
    alpaPDs_hueALL = [alpaPDs_hueALL; wrapTo360(rad2deg(pd(:)))];
    pd = squeeze(circ_mean(tutuPD{nhue}));
    tutuPDs_hueALL = [tutuPDs_hueALL; wrapTo360(rad2deg(pd(:)))];
end
%% PDs (each trial) vs Gamma power (trial-averaged, hue-electrode pair)
% figure; subplot(2,1,1); boxplot(alpaPDs_hueALL, floor(alpaGpower_elecALL(:)/0.5)*0.5); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuPDs_hueALL, floor(tutuGpower_elecALL(:)/0.5)*0.5); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
%% PDs (each trial) vs Gamma power (hues)
% figure; subplot(2,1,1); boxplot(alpaPDs_hueALL, (alpaGpower_hueALL(:)/0.5)*0.5); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuPDs_hueALL, (tutuGpower_hueALL(:)/0.5)*0.5); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');

%%
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_hueALL(:)/0.5)*0.5,alpaPDs_hueALL,[],colors(alpahueind_elecALL,:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_hueALL(:)/0.5)*0.5,tutuPDs_hueALL,[],colors(tutuhueind_elecALL,:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
%%
alpaPDsALL_offset = alpaPDs_hueALL;
for elecind = union(alpaelecind_elecALL(:),[])'
    for hueind = union(alpahueind_elecALL(:),[])'
        sel = (alpaelecind_elecALL==elecind) & (alpahueind_elecALL==hueind);
        pds = deg2rad(alpaPDsALL_offset(sel)); % radians
        pdsMean = circ_mean(pds);
        pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo360(rad2deg(pdsMean));
        alpaPDsALL_offset(sel(:)) = pds(:);
    end
end
tutuPDsALL_offset = tutuPDs_hueALL;
for elecind = union(tutuelecind_elecALL(:),[])'
    for hueind = union(tutuhueind_elecALL(:),[])'
        sel = (tutuelecind_elecALL==elecind) & (tutuhueind_elecALL==hueind);
        pds = deg2rad(tutuPDsALL_offset(sel)); % radians
        pdsMean = circ_mean(pds);
        pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo360(rad2deg(pdsMean));
        tutuPDsALL_offset(sel(:)) = pds(:);
    end
end
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_elecALL(:)/0.5)*0.5,alpaPDsALL_offset, [], colors(alpahueind_elecALL,:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_elecALL(:)/0.5)*0.5,tutuPDsALL_offset, [], colors(tutuhueind_elecALL,:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
% %%
% alpaPDsALLround_offset = alpaPDs_hueALL;
% bin = 0.5;
% for ind = union(floor(alpaGpower_hueALL(:)/bin)*bin,[])'
%     sel = (floor(alpaGpower_hueALL(:)/bin)*bin==ind);
%     pds = deg2rad(alpaPDsALLround_offset(sel)); % radians
%     pdsMean = circ_median(pds);
%     pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo360(rad2deg(pdsMean));
%     alpaPDsALLround_offset(sel(:)) = pds(:);
%     
% end
% tutuPDsALLround_offset = tutuPDs_hueALL;
% for ind = union(floor(tutuGpower_hueALL(:)/bin)*bin,[])'
%     sel = (floor(tutuGpower_hueALL(:)/bin)*bin==ind);
%     pds = deg2rad(tutuPDsALLround_offset(sel)); % radians
%     pdsMean = circ_median(pds);
%     pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo360(rad2deg(pdsMean));
%     tutuPDsALLround_offset(sel(:)) = pds(:);
%     
% end
% figure; subplot(2,1,1); boxplot(alpaPDsALLround_offset, floor(alpaGpower_hueALL(:)/bin)*bin); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuPDsALLround_offset, floor(tutuGpower_hueALL(:)/bin)*bin); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
% 
% %%
% %%
% alpaPDsALLround_offset = alpaPDs_hueALL;
% alpap_rtest_Gpowerbins = [];
% bin = (max(alpaGpower_elecMean)-min(alpaGpower_elecMean))/25;
% for ind = union(floor(alpaGpower_elecALL(:)/bin)*bin,[])'
%     sel = (floor(alpaGpower_elecALL(:)/bin)*bin==ind);
%     pds = deg2rad(alpaPDsALLround_offset(sel)); % radians
%     alpap_rtest_Gpowerbins = [alpap_rtest_Gpowerbins, circ_rtest(deg2rad(alpaPDsALLround_offset(sel)))];
%     pdsMean = circ_mean(pds);
%     pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo180(rad2deg(pdsMean));
%     alpaPDsALLround_offset(sel(:)) = pds(:);
%     
% end
% 
% bin = (max(tutuGpower_elecMean)-min(tutuGpower_elecMean))/25;
% tutuPDsALLround_offset = tutuPDs_hueALL;
% tutup_rtest_Gpowerbins = [];
% for ind = union(floor(tutuGpower_elecALL(:)/bin)*bin,[])'
%     sel = (floor(tutuGpower_elecALL(:)/bin)*bin==ind);
%     pds = deg2rad(tutuPDsALLround_offset(sel)); % radians
%     tutup_rtest_Gpowerbins = [tutup_rtest_Gpowerbins, circ_rtest(deg2rad(tutuPDsALLround_offset(sel)))];
%     pdsMean = circ_mean(pds);
%     pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo180(rad2deg(pdsMean));
%     tutuPDsALLround_offset(sel(:)) = pds(:);
%     
% end
% figure; subplot(2,1,1); boxplot(alpaPDsALLround_offset, floor(alpaGpower_elecALL(:)/bin)*bin); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuPDsALLround_offset, floor(tutuGpower_elecALL(:)/bin)*bin); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
% 
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_elecALL(:)/0.5)*0.5,alpaPDsALLround_offset, [], colors(alpahueind_elecALL,:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_elecALL(:)/0.5)*0.5,tutuPDsALLround_offset, [], colors(tutuhueind_elecALL,:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
%%
% %% unbinned
% alpaPDsALLround_offset = alpaPDs_hueALL;
% bin = 0.5;
% for ind = union((alpaGpower_elecALL(:)/bin)*bin,[])'
%     sel = ((alpaGpower_elecALL(:)/bin)*bin==ind);
%     pds = deg2rad(alpaPDsALLround_offset(sel)); % radians
%     pdsMean = circ_median(pds);
%     pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo180(rad2deg(pdsMean));
%     alpaPDsALLround_offset(sel(:)) = pds(:);
%     
% end
% tutuPDsALLround_offset = tutuPDs_hueALL;
% for ind = union((tutuGpower_elecALL(:)/bin)*bin,[])'
%     sel = ((tutuGpower_elecALL(:)/bin)*bin==ind);
%     pds = deg2rad(tutuPDsALLround_offset(sel)); % radians
%     pdsMean = circ_median(pds);
%     pds = wrapTo180(rad2deg(pds-pdsMean))+wrapTo180(rad2deg(pdsMean));
%     tutuPDsALLround_offset(sel(:)) = pds(:);
%     
% end
% figure; subplot(2,1,1); boxplot(alpaPDsALLround_offset, (alpaGpower_elecALL(:)/bin)*bin); title('M1'); grid on;
% subplot(2,1,2); boxplot(tutuPDsALLround_offset, (tutuGpower_elecALL(:)/bin)*bin); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
% 
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_elecALL(:)/0.5)*0.5,alpaPDsALLround_offset, [], colors(alpahueind_elecALL,:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_elecALL(:)/0.5)*0.5,tutuPDsALLround_offset, [], colors(tutuhueind_elecALL,:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');

% %%
% figure; subplot(2,1,1); scatter(alpaGpower_elecALL(:),alpaPDs_hueALL, '.'); title('M1'); grid on;
% subplot(2,1,2); scatter(tutuGpower_elecALL(:),tutuPDs_hueALL, '.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
%%
% figure; subplot(2,1,1); scatter(floor(alpaGpower_elecALL(:)/0.5)*0.5,alpaPDs_hueALL, '.'); title('M1'); grid on;
% subplot(2,1,2); scatter(floor(tutuGpower_elecALL(:)/0.5)*0.5,tutuPDs_hueALL, '.'); title('M2'); grid on; xlabel('Gamma power (dB)');

% alpaGpowerBin = round(alpaGpower_elecALL(:)/2)*2;
% alpaGpowerBinscaled = min(alpaGpowerBin):2:max(alpaGpowerBin+2);
% for p = alpaGpowerBinscaled
%     
% end
% tutuGpowerBin = round(tutuGpower_elecALL(:)/2)*2;
% tutuGpowerBinscaled = min(tutuGpowerBin):2:max(tutuGpowerBin+2);
selalpa145 = alpaGpower_elecALL >= prctile(alpaGpower_elecALL,10); seltutu18 = tutuGpower_elecALL >= prctile(tutuGpower_elecALL,10); seltutu145 = tutuGpower_elecALL >= prctile(tutuGpower_elecALL,10);
ap = intersect(alpahueind_elecALL(selalpa145),1:numel(alpaPD)); tt = intersect(tutuhueind_elecALL(seltutu145),1:numel(tutuPD)); tt2 = intersect(tutuhueind_elecALL(seltutu18),1:numel(tutuPD));
apmathue = (alpahueind_elecALL(selalpa145)==1:numel(alpaPD)); ttmathue = (tutuhueind_elecALL(seltutu145)==(1:numel(tutuPD))); tt2mathue = (tutuhueind_elecALL(seltutu18)==(1:numel(tutuPD)));
apmatelec = selalpa145 & (alpaelecind_elecALL(:)==1:64); ttmatelec = seltutu145 & (tutuelecind_elecALL(:)==(1:16)); tt2matelec = seltutu18 & (tutuelecind_elecALL(:)==(1:16));
alpaselelectrodes145 = cell(numel(alpaGpower),1); nalpaselelectrodes145 = [];
tutuselelectrodes145 = cell(numel(alpaGpower),1); ntutuselelectrodes145 = [];
tutuselelectrodes18 = cell(numel(alpaGpower),1); ntutuselelectrodes18 = [];
for nhue = 1:numel(alpaGpower)
    alpaselelectrodes145{nhue} = union(alpaelecind_elecALL(apmathue(:,nhue)),[]);
    al = alpaselelectrodes145{nhue};
    nalpaselelectrodes145 = [nalpaselelectrodes145; al(:)];
    tutuselelectrodes145{nhue} = union(tutuelecind_elecALL(ttmathue(:,nhue)),[]);
    tu = tutuselelectrodes145{nhue};
    ntutuselelectrodes145 = [ntutuselelectrodes145; tu(:)];
    tutuselelectrodes18{nhue} = union(tutuelecind_elecALL(tt2mathue(:,nhue)),[]);
    tu = tutuselelectrodes18{nhue};
    ntutuselelectrodes18 = [ntutuselelectrodes18; tu(:)];
end

alpaselhues145 = cell(size(apmatelec,2),1); nalpaselhues145 = [];
alpaselhues145ALL = [];
for nelec = 1:size(apmatelec,2)
    alpaselhues145{nelec} = union(alpahueind_elecALL(apmatelec(:,nelec)),[]);
    al = alpaselhues145{nelec};
    alpaselhues145ALL = [alpaselhues145ALL; al(:)];
    nalpaselhues145 =[nalpaselhues145, numel(alpaselhues145{nelec})];
end
tutuselhues145 = cell(size(ttmatelec,2),1); ntutuselhues145 = [];
tutuselhues145ALL = [];
tutuselhues18 = cell(size(ttmatelec,2),1); ntutuselhues18 = [];
tutuselhues18ALL = [];
for nelec = 1:size(ttmatelec,2)
    tutuselhues145{nelec} = union(tutuhueind_elecALL(ttmatelec(:,nelec)),[]);
    tu = tutuselhues145{nelec};
    tutuselhues145ALL = [tutuselhues145ALL; tu(:)];
    ntutuselhues145 =[ntutuselhues145, numel(tutuselhues145{nelec})];
    tutuselhues18{nelec} = union(tutuhueind_elecALL(tt2matelec(:,nelec)),[]);
    tu = tutuselhues18{nelec};
    tutuselhues18ALL = [tutuselhues18ALL; tu(:)];
    ntutuselhues18 =[ntutuselhues18, numel(tutuselhues18{nelec})];
end
%%
%%
% colors = hsv(36);
% figure; subplot(2,1,1); scatter((alpaGpower_elecALL(selalpa145)/0.5)*0.5,alpaPDs_hueALL(selalpa145), [], colors(alpahueind_elecALL(selalpa145),:),'.'); title('M1'); grid on;
% subplot(2,1,2); scatter((tutuGpower_elecALL(seltutu145)/0.5)*0.5,tutuPDs_hueALL(seltutu145), [], colors(tutuhueind_elecALL(seltutu145),:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
% %%
% colors = hsv(36);
% alral = [];
% selindices = find(selalpa145);
% for ind = selindices'
%     al = circ_rtest(alpaPDs_hueALL(ind));
%     alral = [alral; al];
% end
% figure; subplot(2,1,1); scatter((alpaGpower_elecALL(selalpa145)/0.5)*0.5,alral, [], colors(alpahueind_elecALL(selalpa145),:),'.'); title('M1'); grid on;
% 
% tural = [];
% selindices = find(seltutu145);
% for ind = selindices'
%     tu = circ_rtest(tutuPDs_hueALL(ind));
%     tural = [tural; tu];
% end
% subplot(2,1,2); scatter((tutuGpower_elecALL(seltutu145)/0.5)*0.5,tural, [], colors(tutuhueind_elecALL(seltutu145),:),'.'); title('M2'); grid on; xlabel('Gamma power (dB)');
% sgtitle('PDs');
% %%
% figure; subplot(3,1,1); histogram(alpaselhues145ALL); subplot(3,1,2); histogram(tutuselhues145ALL); subplot(3,1,3); histogram(tutuselhues18ALL);
% figure; subplot(3,1,1); histogram(nalpaselelectrodes145); subplot(3,1,2); histogram(ntutuselelectrodes145); subplot(3,1,3); histogram(ntutuselelectrodes18);

%%
alpaselhues145;
tutuselhues18;
tutuselhues145;

% alpafig = figure;
% plot_1row = 12;
alpaselelecs = [];
alpaselhues = cell(size(alpaPD));
alpap_wwtest = zeros(1,numel(alpaselhues145))/0;
alpap_rtest = zeros(1,numel(alpaselhues145))/0;
alpa_mean_pds = zeros(1,numel(alpaselhues145))/0;
alpa_median_pds = zeros(1,numel(alpaselhues145))/0;
alpa_huemean_pds = cell(1,numel(alpaselhues145));
alpa_mean_pds_ALL = [];
cellcnt = 1;
for elec = 1:numel(alpaselhues145)
    huesel = alpaselhues145{elec};
    newhuesel = [];
    for hue = huesel(:)'
        p_ral = circ_rtest(deg2rad(alpaPDs_hueALL((alpahueind_elecALL==hue) & (alpaelecind_elecALL == elec) & (selalpa145))));
        if p_ral<0.01
            newhuesel = [newhuesel, hue];
        end
    end
    huesel = newhuesel;
    if ~isempty(huesel)
        alpaselelecs = [alpaselelecs, elec];
        alpaselhues{cellcnt} = huesel;
        cellcnt = cellcnt+1;
        if numel(huesel) > 2
            disp(['elec ',num2str(elec)]);
            nboot = min(ceil(numel(huesel)/2),10);
            bootfun = @(x) wwbootfun(x, alpaPDs_hueALL, apmatelec, alpahueind_elecALL, elec);
            [pci, pvals] = bootci(nboot, {bootfun, huesel}, 'type','norm','alpha', 0.05);
            alpap_wwtest(elec) = prctile(pvals,99);
            disp(['elec ',num2str(elec)]);
        elseif numel(huesel) == 2
            alpap_wwtest(elec) = circ_wwtest(deg2rad(alpaPDs_hueALL(apmatelec(huesel,elec))), alpahueind_elecALL(apmatelec(huesel,elec)));
        else
            alpap_wwtest(elec) = inf;
        end
        alpa_mean = [];
        for hue = huesel(:)'
            alpa_mean = [alpa_mean, wrapTo360(rad2deg(circ_mean(deg2rad(alpaPDs_hueALL((alpahueind_elecALL==hue) & (alpaelecind_elecALL == elec) & (selalpa145))))))];
        end
        alpa_huemean_pds{elec} = alpa_mean;
        alpa_mean_pds_ALL = [alpa_mean_pds_ALL;alpa_mean(:)];
        alpap_rtest(elec) = circ_rtest(deg2rad(alpaPDs_hueALL(apmatelec(:,elec))));
        alpa_mean_pds(elec) = wrapTo360(rad2deg(circ_mean(deg2rad(alpaPDs_hueALL(apmatelec(:,elec))))));
        alpa_median_pds(elec) = wrapTo360(rad2deg(circ_median(deg2rad(alpaPDs_hueALL(apmatelec(:,elec))))));
    end
end
% tutufig = figure;
% plot_1row = 12;
% sgtitle('PDs');
%%
tutuselelecs = [];
tutuselhues = cell(size(tutuPD));
% tutuselhues = cell(1,numel(tutuselhues145));
tutup_wwtest = zeros(1,numel(tutuselhues145))/0;
tutup_rtest = zeros(1,numel(tutuselhues145))/0;
tutu_mean_pds = zeros(1,numel(tutuselhues145))/0;
tutu_median_pds = zeros(1,numel(tutuselhues145))/0;
tutu_huemean_pds = cell(1,numel(tutuselhues145));
tutu_mean_pds_ALL = [];
cellcnt = 1;
for elec = 1:numel(tutuselhues18)
    huesel = tutuselhues18{elec};
    newhuesel = [];
    for hue = huesel(:)'
        p_ral = circ_rtest(deg2rad(tutuPDs_hueALL((tutuhueind_elecALL==hue) & (tutuelecind_elecALL == elec) & (seltutu18))));
        if p_ral<0.01
            newhuesel = [newhuesel, hue];
        end
    end
    huesel = newhuesel;
    if ~isempty(huesel)
        tutuselelecs = [tutuselelecs, elec];
        tutuselhues{cellcnt} = huesel;
        cellcnt = cellcnt + 1;
        if numel(huesel) > 2
            disp(['elec ',num2str(elec)]);
%             tutup_wwtest(elec) = circ_wwtest(deg2rad(tutuPDs_hueALL(tt2matelec(:,elec))), tutuhueind_elecALL(tt2matelec(:,elec)));
            nboot = min(ceil(numel(huesel)/2),10);
            bootfun = @(x) wwbootfun(x, tutuPDs_hueALL, tt2matelec, tutuhueind_elecALL, elec);
            [pci, pvals] = bootci(nboot, {bootfun, huesel}, 'type','norm','alpha', 0.05);
            tutup_wwtest(elec) = prctile(pvals,99);
            
            disp(['elec ',num2str(elec)]);
        elseif numel(huesel) == 2
            tutup_wwtest(elec) = circ_wwtest(deg2rad(tutuPDs_hueALL(tt2matelec(huesel,elec))), tutuhueind_elecALL(tt2matelec(huesel,elec)));
        else    
            tutup_wwtest(elec) = inf;
        end
        tutu_mean = [];
        for hue = huesel(:)'
            tutu_mean = [tutu_mean, wrapTo360(rad2deg(circ_mean(deg2rad(tutuPDs_hueALL((tutuhueind_elecALL==hue) & (tutuelecind_elecALL == elec) & (seltutu18))))))];
        end
        tutu_huemean_pds{elec} = tutu_mean;
        tutu_mean_pds_ALL = [tutu_mean_pds_ALL;tutu_mean(:)];
        tutup_rtest(elec) = circ_rtest(deg2rad(tutuPDs_hueALL(tt2matelec(:,elec))));
        tutu_mean_pds(elec) = wrapTo360(rad2deg(circ_mean(deg2rad(tutuPDs_hueALL(tt2matelec(:,elec))))));
        tutu_median_pds(elec) = wrapTo360(rad2deg(circ_median(deg2rad(tutuPDs_hueALL(tt2matelec(:,elec))))));
    end
end
%%
% figure;subplot(221);histogram(tutu_median_pds(tutuselelecs),'binwidth',10);subplot(222);histogram(tutu_mean_pds(tutuselelecs),'binwidth',10);
% subplot(2,4,[6,7]);histogram(alpa_mean_pds(alpaselelecs),'binwidth',10);




















































































































