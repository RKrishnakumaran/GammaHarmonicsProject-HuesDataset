%% bootfun function for wwtest bootstrap
function pval = wwbootfun(x, PDs_hueALL, matelec, hueind_elecALL, elec)
    sel = any(hueind_elecALL(:) == x(:)',2);
    pval = circ_wwtest(deg2rad(PDs_hueALL(matelec(sel,elec))), hueind_elecALL(matelec(sel,elec)));
    if isnan(pval)
       disp('Zoinks!'); 
    end
end