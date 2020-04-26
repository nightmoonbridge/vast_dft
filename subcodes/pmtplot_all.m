function pmtplot_all(general, metric_all, dummy, labels)
if nargin < 4
    labels.titletext = '';
    labels.figname1 = 'fVASTnarrow';
    labels.figname2 = 'fVASTbroad';
    labels.tightfig = false;
end

if ~isfield(labels,'tightfig')
    tightfig = false;
else
    tightfig = labels.tightfig;
end

activeidx = find(~cellfun('isempty',metric_all));

legendnames = general.legendnames(activeidx);

activemetrics = metric_all(activeidx);

ncandidates = length(activeidx);

titletext = labels.titletext;

colorpool = lines(ncandidates);

dF = general.fs/general.lenConFilter;

if general.incl_dcnyq
    fxtk = (0:dF:general.fs/2)';
else
    fxtk = (dF:dF:general.fs/2-dF)';
end

zonelabels = {'Zone alpha', 'Zone beta'}';

fontweight = 'Bold';
fontsize = 13;
if tightfig
    figsize = [10 10 2*[8.89 5]];
    figure('Name', titletext, ...
        'Units', 'Centimeters', 'Position', figsize)
else
    figure('Position',[450 170 900 650], 'Name', titletext)
end

for zidx = 1:2
    subplot(2,1,zidx)
    hold on
    for spidx = 1:ncandidates
        plot(fxtk, 10*log10(activemetrics{spidx}{zidx}), ...
            'DisplayName', legendnames{spidx}, 'LineStyle', '-', ...
            'Color', colorpool(spidx,:), 'LineWidth', 2)
    end
    hold off
    grid on
    ll = legend;
    xlabel('Frequency (Hz)')
    set(ll,'location','best')
    legendhitcallback
    set(gca, 'XScale', 'log', 'Box','on', 'XLim', [100 general.fs/2], ...
        'LooseInset', max(get(gca,'TightInset'), 0.02), ...
        'XTick', [100 1000 2000 4000 8000], 'XTickLabel', ...
        {'100', '1k', '2k', '4k', '8k'},...
        'FontWeight', fontweight, 'FontSize', fontsize)
    title(zonelabels{zidx})
end
end