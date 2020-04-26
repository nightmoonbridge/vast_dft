function fxInputSignals = getFxInputSignals(general, array, zone, indata, ctrfilt)

fprintf('------- Linear convolution reconstruction is being used. -------\n')
% Filtered input signals and the filter will be the control filter
fxInputSignals = cellfun(@(x) zeros(general.lenInput, ...
    array.numLoudspk), cell(zone.number,1), 'UniformOutput', false);

% Time invariant system
for sridx = 1:zone.number
    linidx = 1:general.lenConFilter;
    xin = indata{sridx}.xin;
    for lidx = 1:array.numLoudspk
        fxInputSignals{sridx}(:,lidx) = filter(ctrfilt.conFilter{sridx}(linidx),1,xin);
        linidx = linidx + general.lenConFilter;
    end
end

end