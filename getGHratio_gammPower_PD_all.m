% Compute GHratio, gammapower (dB) and phaseDifference for both hue cases and achromatic
% takes data from savedData/processedData 
basePath = pwd;
subjectName = 'tutu'; % Change subjectname here

[ratioGH, powerGamma, freqGamma, freqHarmonic, PD] = deal({});
% Color
for i = 1:36
    hue = (i-1)*10;
    Hue = num2str(hue);
    [ratioGH{i}, powerGamma{i}, freqGamma{i}, freqHarmonic{i}, PD{i}] = extractData(subjectName, 'Color', Hue, basePath);
    disp(hue);
end
% Achromatic
i = 37;
[ratioGH{i}, powerGamma{i}, freqGamma{i}, freqHarmonic{i}, PD{i}] = extractData(subjectName, 'Achro', '360', basePath);

save(fullfile(basePath, 'Data', [subjectName 'GH' '.mat']), 'ratioGH', 'powerGamma', 'freqGamma', 'freqHarmonic');
save(fullfile(basePath, 'Data', [subjectName 'PD' '.mat']), 'PD');

function [ratioGH, powerGamma, freqGamma, freqHarmonic, PD] = extractData(subjectName, expType, stimType, basePath)
    saveFolder = fullfile(basePath, 'savedData', 'processedData', subjectName);
    gammaRangeHz = [30 70]; % Hz
    maxHarmFreq = 140; % Hz
    GHwidth = 12; % Hz
    if strcmp(expType,'Achro')
        load(fullfile(saveFolder, [subjectName 'Achro' '.mat']), 'psdST', 'psdBL', 'freqVals', 'phaseDiff');
    else
        load(fullfile(saveFolder, [subjectName expType stimType '.mat']), 'psdST', 'psdBL', 'freqVals', 'phaseDiff');
    end
    
    [~, t, e] = size(psdST);
    [ratioGH, powerGamma, freqGamma, freqHarmonic] = deal([]);
    for k = 1:e
        powerDB = 10*(log10(mean(psdST(:,:,k),2)) - log10(mean(psdBL(:,:,k),2)));
        
        [pk, loc] = findpeaks(powerDB);
        pkG = findpeaks(powerDB(gammaRangeHz(1)/2+1:gammaRangeHz(2)/2+1));

        pGG = max(pkG);
        fG = freqVals(loc(find(pk == pGG)));
        lG = loc(find(pk == pGG));
        pkH = pk(intersect(find(loc > lG+GHwidth/2), find(loc < maxHarmFreq/2+1)));
        pHH = max(pkH);
        fH = freqVals(loc(find(pk == pHH)));

        ratioGH(k) = fH/fG;
        freqGamma(k) = fG;
        freqHarmonic(k) = fH;
        powerGamma(k) = pGG;
    end
    
    PD = phaseDiff;
end
