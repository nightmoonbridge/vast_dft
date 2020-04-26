function [irVirsrc_full, irVirsrc] = getVirsrcrir(general, array, zone, irMeasured_full)
vpos = array.virsrcpos;
nrir = array.nsamples_rirs;

irVirsrc_reference = cell(zone.number,1);
for sridx = 1:zone.number
    irVirsrc_reference{sridx} = generateImpulseResponseVirtualSource(general, array, zone, sridx);
end

vir_peak_pool = zeros(zone.numCtrPts,zone.number);
for sridx = 1:zone.number
    for moidx = 1:zone.numCtrPts
        vir_peak_pool(moidx,sridx) = find(abs(irVirsrc_reference{sridx}(:,moidx)) == max(abs(irVirsrc_reference{sridx}(:,moidx))));
    end
end

lenir = general.lenRir*general.sfactor;
lidx = (1:lenir) + (vpos-1)*lenir;

lenir_rsmpl = general.lenRir;

irVirsrc_full = cellfun(@(x) zeros(lenir, zone.numCtrPts), ...
    cell(zone.number,1), 'UniformOutput', false);
irVirsrc = cellfun(@(x) zeros(lenir_rsmpl, zone.numCtrPts), ...
    cell(zone.number,1), 'UniformOutput', false);

for zidx = 1:zone.number
    for midx = 1:zone.numCtrPts
        virorg = irMeasured_full{zidx}(lidx,midx);
        virpeakidx = find(abs(virorg) == max(abs(virorg)));
        
        addzerosidx = vir_peak_pool(midx,zidx) - virpeakidx;
        
        comp_rir = [zeros(addzerosidx,1); virorg(1:end-addzerosidx)];
        
        partidx = 1:vir_peak_pool(midx,zidx)+nrir;
        irVirsrc_full{zidx}(partidx,midx) = comp_rir(partidx);
        if nargout > 1
            irVirsrc{zidx}(:,midx) = general.sfactor*resample(irVirsrc_full{zidx}(:,midx),general.fs,general.tfs);
        end
    end
end
end

