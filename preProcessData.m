basePath = pwd; 

% For all hues (0:10:350): i = 1 => 0 i = 36 => 350
% for i = 1:36
%     hue = (i-1)*10;
%     Hue = num2str(hue);
% %     dataPreProcessing('alpa', 'Color', Hue, basePath, 0);
%     dataPreProcessing('tutu', 'Color', Hue, basePath, 0);
% %     dataPreProcessing('kesari', 'Color', Hue, basePath, 0);
%     disp(hue);
% end

% % Achromatic (SForiAchro for alpa and '37th' hue for tutu)
% disp('Processing achromatic data')
% dataPreProcessing('tutu', 'Color', '360', basePath, 1);
% dataPreProcessing('alpa', 'SForiAchro', '360', basePath, 1);

% reds = [1:3, 34:36];
% for i = 1:length(reds)
%     hue = (reds(i)-1)*10;
%     Hue = num2str(hue);
%     dataPreProcessing('M4', 'Color', Hue, basePath, 0);
% end

for i = 1:37
    hue = (i-1)*10;
    Hue = num2str(hue);
    dataPreProcessing('M4', 'Color', Hue, basePath, i == 37);
    disp(hue);
end

fig1

function dataPreProcessing(subjectName, expType, stimType, basePath, achroFlag)

    % Get Corresponding experiment/data
    subjectID = [subjectName expType];
    isnoise = 0;
    if strcmp(subjectID,'alpaColor')
        subjectName = 'alpa';expDate = '301215'; protocolName = 'GRF_001'; % 488: Hue fullscreen
    elseif strcmp(subjectID,'tutuColor')
        subjectName = 'tutu'; expDate = '191016'; protocolName = 'GRF_001'; % 111: Hue fullscreen
    elseif strcmp(subjectID,'kesariColor')
        subjectName = 'kesari'; expDate = '050516'; protocolName = 'GRF_001'; % 307: Hue fullscreen
    elseif strcmp(subjectID,'alpaSForiAchro')
        subjectName = 'alpa'; expDate = '301215'; protocolName = 'GRF_005'; % SFOri - alpa Vinay
    elseif strcmp(subjectID,'kesariSForiAchro')
        subjectName = 'kesari'; expDate = '030516'; protocolName = 'GRF_002'; % SFOri - kesari Vinay
    
    saveFolder = fullfile(basePath, 'savedData', 'processedData', subjectName);
    mkdir(saveFolder);    
    
    elseif strcmp(subjectName,'M4')
        subjectName = 'tutu'; expDate = '191016'; protocolName = 'GRF_001'; % 111: Hue fullscreen
        isnoise = 0.8; %0.9
        factor = 0.8; % factor = .8, don't chnage, offset doesn't change with isnoise
        saveFolder = fullfile(basePath, 'savedData', 'processedData', 'M4');
        mkdir(saveFolder); 
    end
    
%     saveFolder = fullfile(basePath, 'savedData', 'processedData', subjectName);
%     mkdir(saveFolder);
    gridType = 'Microelectrode';
    folderSourceString = fullfile(basePath, expType);
    folderBase = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);
    folderLFP = fullfile(folderBase,'segmentedData','LFP');
    folderSpikes = fullfile(folderBase,'segmentedData', 'Spikes');

    stimPeriod = [0.25 0.75];% 500ms
    baselinePeriod = [-0.5 0]; %500 ms

    % TimeVals, FS and StimPos
    load(fullfile(folderLFP,'lfpInfo'),'timeVals');
    Fs = 1./(timeVals(2)-timeVals(1));
    numPoints = round(diff(stimPeriod)*Fs);

    stPos = find(timeVals>=stimPeriod(1),1) + (1:numPoints);
    blPos = find(timeVals>=baselinePeriod(1),1) + (1:numPoints);

    % Multi-taper parameters
    if ~exist('TW','var'); TW = 2; end
    tw = TW; fmax = 250;
    mt.tapers = [tw (2*tw-1)];
    mt.pad = -1; mt.Fs = Fs;
    mt.fpass = [0 fmax];
    
    addpath(genpath([basePath '/Programs']))
   
%     % load badtrials
    badTrialFile = fullfile(folderBase,'segmentedData','badTrials');
    load(badTrialFile,'badTrials');
    paramCombinationsFile = fullfile(folderBase,'extractedData','parameterCombinations');
    Params = load(paramCombinationsFile);
    
    if strcmp(subjectID,'alpaSForiAchro')
        highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes', 'alpaLFPElectrodeList.mat'));
        highRMSElectrodes = highRMSElectrodesStruct.alpaLFPElectrodeList;
        highRMSElectrodes = setdiff(highRMSElectrodes, 4); % 4 is an electrode compared to color cases (65 v 64)
        
        cVals = Params.cValsUnique; oVals = Params.oValsUnique; fVals = Params.fValsUnique;
%         c = find(cVals ==100); % contrast 100
        o = find(oVals ==90); %90 degree
%         f = find(round(fVals) == 2); % SF = 2
        goodPos = Params.parameterCombinations{1,1,1,1,o};
        goodPos = setdiff(goodPos, badTrials);
        
    elseif strcmp(subjectID,'kesariSForiAchro')
        highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes', 'kesariLFPElectrodeList.mat'));
        highRMSElectrodes = highRMSElectrodesStruct.kesariLFPElectrodeList;
        highRMSElectrodes = setdiff(highRMSElectrodes, 4); % 4 is an electrode compared to color cases (65 v 64)
        
        cVals = Params.cValsUnique; oVals = Params.oValsUnique; fVals = Params.fValsUnique;
%         c = find(cVals ==100); % contrast 100
        o = find(oVals ==90); %90 degree
%         f = find(round(fVals) == 2); % SF = 2
        goodPos = Params.parameterCombinations{1,1,1,1,o};
        goodPos = setdiff(goodPos, badTrials);
        
    else
        highRMSElectrodesStruct = load(fullfile(folderSourceString,'analyzeElectrodes',subjectName,'highRMSElectrodes'));
        highRMSElectrodes = highRMSElectrodesStruct.highRMSElectrodes;
        
        cVals = Params.cValsUnique; oVals = Params.oValsUnique;
        goodPosAll = cell(1,length(oVals));
        c = find(cVals ==100); % contrast 100
        for o = 1:length(oVals)
            goodPosAll{o} = Params.parameterCombinations{1,1,1,1,o,c,1};
            goodPosAll{o} = setdiff(goodPosAll{o},badTrials);
        end
        stimIndex = str2double(stimType)/10 + 1;
        goodPos = goodPosAll{stimIndex};
        % stimIndex - represents Color number, ranging 1 to 36
        % representing hues 0 to 350 with interval size 10
    end

    % psdST = 3D Array - (:, trialIndex, electrodeIndex) $#

    n = 4; % order of butterWorth Filtering of gamma and harmonic signals
    delta = 10; % for filtfilt signal filtering
    gammaRangeHz = [30 70]; % Hz
    maxHarmFreq = 140; %Hz
    GHwidth = 12; % Hz
   
    phaseDiff = []; Gamma_phase = [];
    psdST = []; psdBL = [];
    gammaFreq = []; harmonicFreq = []; harmonicFreq_act = [];
    baseCorrLog10PSD = [];
    
    stTimevals = timeVals(stPos);
    gphaseCollect = cell(length(highRMSElectrodes), length(goodPos));
    stspikeTimeCollect = cell(length(highRMSElectrodes), length(goodPos));
    blspikeTimeCollect = cell(length(highRMSElectrodes), length(goodPos));


    for i = 1:length(highRMSElectrodes)
        load(fullfile(folderLFP,['elec' num2str(highRMSElectrodes(i))]),'analogData');
        load(fullfile(folderSpikes,['elec' num2str(highRMSElectrodes(i)) '_SID0']),'spikeData');
        % setdiff (0, 255)
        
        if isnoise
%             stLFPP = analogData(goodPos,stPos)';
            blLFPP = analogData(goodPos,blPos)';
%               stLFPP = (isnoise)*10^1.15*randn(size(analogData(goodPos,stPos)')) + (1-isnoise)*analogData(goodPos,stPos)';
%               blLFPP = (isnoise)*10^1.15*randn(size(analogData(goodPos,stPos)')) + (1-isnoise)*analogData(goodPos,blPos)';

              stLFPP = factor*(isnoise)*blLFPP + (1-isnoise)*analogData(goodPos,stPos)';
              
%             we will have a mixing coefficient (rho) 
%               we will generate white noise s.t it's psd is < 75% of
%               minimum baseline psd
%               noise is linearly mixed accordingly with rho, 
%               new signal = rho*signal + (1-rho)*noise
%               rho will range from 1 to 0.1 in log scale, [1, 1/3, 1/6, 1/12]
%             disp('M4')
        else
            stLFPP = analogData(goodPos,stPos)';
            blLFPP = analogData(goodPos,blPos)';
        end
        
        stSpikes = cellfun(@(x) x(x>=stimPeriod(1) & x<=stimPeriod(2)), spikeData, 'UniformOutput', false);
        blSpikes = cellfun(@(x) x(x>=baselinePeriod(1) & x<=baselinePeriod(2)), spikeData, 'UniformOutput', false);
        
        [psdSTelec,freqVals] = mtspectrumc(stLFPP,mt);
        psdBLelec = mtspectrumc(blLFPP,mt);
        psdST = cat(3,psdST,psdSTelec);
        psdBL = cat(3,psdBL,psdBLelec);
        baseCorrLog10PSD = cat(3,baseCorrLog10PSD,log10(psdSTelec)-log10(psdBLelec));

        % Find gamma and and following(harmonic) peak
        gammaRangePos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<gammaRangeHz(2)));
        gammaRange = baseCorrLog10PSD(gammaRangePos,:,i);
        gammaAmpLog = max(gammaRange);
        [gammaPosInd,~] = find(gammaRange==gammaAmpLog);
        gammaPosss = gammaRangePos(gammaPosInd);
        peakGammaFreq = freqVals(gammaPosss);

        
        peakHarmonicFreq = 2*peakGammaFreq; % Exact Harmonic (2*fG)

        gammaFreq = cat(2,gammaFreq,peakGammaFreq');
        harmonicFreq = cat(2,harmonicFreq,peakHarmonicFreq');

        PD = [];
        for j =1:length(goodPos)
            spike = stSpikes{goodPos(j)};
            stspikeTimeCollect{i, j} = spike;
            
            spike_bl = blSpikes{goodPos(j)};
            blspikeTimeCollect{i, j} = spike_bl;

            [B,A] = butter(n,[peakGammaFreq(j)-delta, peakGammaFreq(j)+delta]/(Fs/2));
            [D,C] = butter(n,[peakHarmonicFreq(j)-delta, peakHarmonicFreq(j)+delta]/(Fs/2));

            gammaSignal = filtfilt(B,A,stLFPP(:,j)); %#
            harmonicSignal = filtfilt(D,C,stLFPP(:,j));

            % Phasedifference between gamma and harmonic
            
            
            G = hilbert(gammaSignal);
            H = hilbert(harmonicSignal);
            PD = cat(2, PD, (angle(H)-2*angle(G)));
            
            if ~isempty(spike)
                absdiff = abs(spike - stTimevals);
                spikeIndices = floor(find(absdiff == min(absdiff, [], 2))/numel(spike)) + 1;
                G_phase = angle(G(spikeIndices));
                gphaseCollect{i, j} = G_phase;
            end
        end
        phaseDiff = cat(3, phaseDiff, PD);
    end
    
% - Uncomment this to save
    if isnoise
       subjectName = 'M4';
    end
    
    if achroFlag == 1
        save(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'phaseDiff','gammaFreq', ...
        'psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
        'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz',...
        'gphaseCollect', 'stspikeTimeCollect', 'blspikeTimeCollect')
    else
        save(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'phaseDiff','gammaFreq', ...
            'psdST', 'psdBL', 'goodPos', 'stPos', 'highRMSElectrodes', 'timeVals', 'freqVals', ...
            'subjectName', 'expType', 'stimType', 'expDate', 'protocolName', 'mt', 'Fs', 'gammaRangeHz', ...
            'gphaseCollect', 'stspikeTimeCollect', 'blspikeTimeCollect')
    end
end
 