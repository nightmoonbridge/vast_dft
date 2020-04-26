function plpl
figs = findobj('Type','figure');
if isempty(figs)
    warning('No figure is found.')
end
for ii = 1:length(figs)
    figure(figs(ii).Number)
end

end