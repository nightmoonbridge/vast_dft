function pfm_mtx = getPerformanceMetrics_v2...
    (general, array, zone, varout, ctrfilt, indata, irArray)

nmethods = varout.nmethods;

varout.ismonitor = false;
ctr = initPerformanceMetrics(general, varout);

for nidx = 1:nmethods
    fxInputSignals = getFxInputSignals(general, array, zone, ...
        indata, ctrfilt{nidx}(1,1));
    
    ctrReprodSoundField = getCReprodSoundField(general, array, ...
        zone, fxInputSignals, irArray, false);
    
    varout.nidx = nidx;
    varout.vidx = 1;
    varout.midx = 1;
    
    ctr = calculatePerformanceMetrics_v2(varout, ctrReprodSoundField, ctr);
end

pfm_mtx.ctr = ctr;

end