function desiredSoundField = getDesiredSoundField(general, zone, indata, irVirsrc)

numPts = size(irVirsrc{1},2);

desiredSoundField = cellfun(@(x) zeros(general.lenInput,numPts), ...
    cell(zone.number,1), 'UniformOutput', false);

for sridx = 1:zone.number
    zonedata = indata{sridx}.xin;
    for midx = 1:numPts
        desiredSoundField{sridx}(:,midx) = filter(irVirsrc{sridx}(:,midx),1,zonedata);
    end
end
end