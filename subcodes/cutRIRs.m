function [irMeasured_cut, irVirsrc_cut] = cutRIRs(general, varout, cutlength, irMeasured, irVirsrc)
irMeasured_cut = cellfun(@(x) zeros(cutlength*varout.nloudspks, varout.nctrpts), ...
    cell(varout.nzones,1), 'UniformOutput', false);
irVirsrc_cut = cellfun(@(x) zeros(cutlength, varout.nctrpts), ...
    cell(varout.nzones,1), 'UniformOutput', false);

for ii = 1:varout.nctrpts
    for jj = 1:varout.nzones
        for kk = 1:varout.nloudspks
            rlidx = (kk-1)*general.lenIRmeasured + (1:general.lenIRmeasured);
            llidx = (kk-1)*cutlength + (1:cutlength);
            temprir = irMeasured{jj}(rlidx,ii);
            irMeasured_cut{jj}(llidx,ii) = temprir(1:cutlength);
        end
        irVirsrc_cut{jj}(:,ii) = irVirsrc{jj}(1:cutlength,ii);
    end
end

end