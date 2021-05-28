load('./Data/alpaPD.mat');
alpaPD = PD;
load('./Data/tutuPD.mat');
tutuPD = PD;
clear PD
%%
load('./Data/alpaGH.mat','powerGamma','ratioGH');
alpa_ghRatio = ratioGH; alpaGpower = powerGamma;
load('./Data/tutuGH.mat','powerGamma','ratioGH');
tutu_ghRatio = ratioGH; tutuGpower = powerGamma;
clear ratioGH; clear powerGamma;
%%
% generate variables for Mean a)gamma and b)harmonic power during each hue, across all trials and electrodes

elecDim = 2; TrialDim = 1;
%% Sort hues by a)gamma power and b)harmonic power
alpaGpower_hue = [];

alpahueind_elecALL = [];
for nhue = 1:numel(alpaGpower)
    alpaGpower_selecting_hues = mean(alpaGpower{nhue});%mean(log(mean(exp(alpaGpower{nhue}),TrialDim)),elecDim);
    alpaGpower_hue = [alpaGpower_hue, alpaGpower_selecting_hues(:)];
    
    phdiff = squeeze(alpaPD{nhue}(1,:,:));
    gh = alpa_ghRatio{nhue};

    ind = nhue*ones(size(phdiff));
    alpahueind_elecALL = [alpahueind_elecALL; ind(:)];
end

[alpahueGpower_sorted, alpahueindex_Gsorted] = sort(alpaGpower_hue,'descend');
hues = hsv(numel(alpa_ghRatio));

alpahueGsorted = hues(alpahueindex_Gsorted, :); 
alpaGpower_hueALL = alpaGpower_hue(alpahueind_elecALL);
[alpahueGpower_ALLsorted, alpahueindex_GALLsorted] = sort(alpaGpower_hueALL(:),'descend');

tutuGpower_hue = []; 

tutuhueind_elecALL = [];
for nhue = 1:numel(tutuGpower)
    tutuGpower_selecting_hues = mean(tutuGpower{nhue});% mean(log(mean(exp(tutuGpower{nhue}),TrialDim)),elecDim);
    tutuGpower_hue = [tutuGpower_hue, tutuGpower_selecting_hues];
    
    phdiff = squeeze(tutuPD{nhue}(1,:,:));
    gh = tutu_ghRatio{nhue};
    ind = nhue*ones(size(phdiff));
    tutuhueind_elecALL = [tutuhueind_elecALL; ind(:)];
end
[tutuhueGpower_sorted, tutuhueindex_Gsorted] = sort(tutuGpower_hue,'descend');
hues = hsv(numel(tutu_ghRatio));

tutuhueGsorted = hues(tutuhueindex_Gsorted, :); %(tutuhueindex_Gsorted)*360/numel(tutu_ghRatio);
tutuGpower_hueALL = tutuGpower_hue(tutuhueind_elecALL);
[tutuhueGpower_ALLsorted, tutuhueindex_GALLsorted] = sort(tutuGpower_hueALL(:),'descend');


%% Sort electrodes by a)gamma power and b)harmonic power
alpaGpower_elecMean = [];
alpaGpower_elecMed = [];
alpaGpower_elecIQR = [];
alpaGH_elecMed = [];
alpaGH_elecIQR = [];
alpaGH_elecALL = []; alpaelecind_elecALL = []; alpaGpower_elecALL = [];
for nhue = 1:numel(alpaGpower)
    pwr = alpaGpower{nhue}; 
    alpaGpower_elecMean = [alpaGpower_elecMean, pwr(:)];
    
    
    phdiff = squeeze(alpaPD{nhue}(1,:,:));
    gh = alpa_ghRatio{nhue};
    alpaGH_elecALL = [alpaGH_elecALL; ones([size(phdiff,1),1])*gh(:)'];
    ind = ones([size(phdiff,1),1])*(1:size(phdiff,2));
    alpaelecind_elecALL = [alpaelecind_elecALL; ind(:)];
end
[alpaelecGpower_sorted, alpaelecindex_Gsorted] = sort(alpaGpower_elecMean(:),'descend');
alpaGpower_elecMean = alpaGpower_elecMean(:);
alpaGpower_elecALL = alpaGpower_elecMean(alpaelecind_elecALL(:)+(alpahueind_elecALL(:)-1)*numel(alpaGpower_elecMean)/numel(alpaGpower));
[alpaelecGpower_ALLsorted, alpaelecindex_GALLsorted] = sort(alpaGpower_elecALL(:),'descend');
tutuGpower_elecMean = [];

tutuGH_elecALL = []; tutuelecind_elecALL = []; tutuGpower_elecALL = [];
for nhue = 1:numel(tutuGpower)
    pwr = tutuGpower{nhue}; 
    tutuGpower_elecMean = [tutuGpower_elecMean, pwr(:)];
    
    phdiff = squeeze(tutuPD{nhue}(1,:,:));
    gh = tutu_ghRatio{nhue};
    tutuGH_elecALL = [tutuGH_elecALL; ones([size(phdiff,1),1])*gh(:)'];
    ind = ones([size(phdiff,1),1])*(1:size(phdiff,2));
    tutuelecind_elecALL = [tutuelecind_elecALL; ind(:)];
end
[tutuelecGpower_sorted, tutuelecindex_Gsorted] = sort(tutuGpower_elecMean(:),'descend');
tutuGpower_elecMean = tutuGpower_elecMean(:);
tutuGpower_elecALL = tutuGpower_elecMean(tutuelecind_elecALL(:)+(tutuhueind_elecALL(:)-1)*numel(tutuGpower_elecMean)/numel(tutuGpower));
[tutuelecGpower_ALLsorted, tutuelecindex_GALLsorted] = sort(tutuGpower_elecALL(:),'descend');

%% PDs
alpaPDs_hueALL = [];
tutuPDs_hueALL = [];
for nhue = 1:numel(tutuGpower)
    pd = squeeze(circ_mean(alpaPD{nhue}));
    alpaPDs_hueALL = [alpaPDs_hueALL; wrapTo360(rad2deg(pd(:)))];
    pd = squeeze(circ_mean(tutuPD{nhue}));
    tutuPDs_hueALL = [tutuPDs_hueALL; wrapTo360(rad2deg(pd(:)))];
end

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
%%
%%

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
alpaselhues145;
tutuselhues18;
tutuselhues145;

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


















































































































