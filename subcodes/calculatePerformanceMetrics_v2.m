function pfm_mtx = calculatePerformanceMetrics_v2...
    (varout, ctrReprodSoundField, pfm_mtx)

nzones = varout.nzones;

nidx = varout.nidx; % nmethods
vidx = varout.vidx; % numV
midx = varout.midx; % nummu

for zidx = 1:nzones                     % the number of zones
    ctrReprodBr = ctrReprodSoundField.Br{zidx};
    ctrReprodDk = ctrReprodSoundField.Dk{zidx};
    
    % Acoustic Contrast
    pfm_mtx.ac_mtx{nidx,zidx}(vidx,midx).scores = ...
        10*log10(norm(ctrReprodBr,'fro')^2/norm(ctrReprodDk,'fro')^2);
end

end