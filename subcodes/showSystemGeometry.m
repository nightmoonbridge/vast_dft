function showSystemGeometry(array,zone,room)
h = zeros(zone.number+3, 1);

colororder = prism(zone.number);

figure('Name','System geometry', 'Position', [600 100 800 800])
h(1) = scatter(array.loudspkPositions(:,1),array.loudspkPositions(:,2),50,...
    'MarkerFaceColor','k','MarkerEdgeColor','none',...
    'DisplayName','loudspks', 'Marker', 'v');
hold on
axis equal
zonelabel = {'Zone A','Zone B'};

cursize = 5;
for zidx = 1:zone.number
    h(zidx+1) = scatter(zone.geometry(zidx).ctrPtsPositions(:,1),zone.geometry(zidx).ctrPtsPositions(:,2),cursize,...
        'DisplayName',zonelabel{zidx},'MarkerEdgeColor','none', ...
        'MarkerFaceColor',colororder(zidx,:));
end

% Plotting the locations of the virtual sources
h(end-1) = scatter(array.virsrcPosition(1,1),array.virsrcPosition(1,2), 70,...
    'markerfacecolor','k','markeredgecolor','none',...
    'displayname','virtSrc A', 'marker', 'pentagram');
h(end) = scatter(array.virsrcPosition(2,1),array.virsrcPosition(2,2), 70,...
    'markerfacecolor','b','markeredgecolor','none',...
    'displayname','virtSrc B', 'marker', 'pentagram');

if strcmpi(zone.ctrPtsType, 'import')
    roomcorners = [0 -room.xlength -room.xlength 0 0;
        0 0 room.ywidth room.ywidth 0]';
    xlimhalf = [-5.5 1];
    ylimhalf = [-0.5 5];
else
    xhalf = room.xlength*0.5;
    yhalf = room.ywidth*0.5;

    roomcorners = [xhalf -xhalf -xhalf xhalf xhalf;
        -yhalf -yhalf yhalf yhalf -yhalf]';

    xlimhalf = ceil(room.xlength)*0.5*[-1 1];
    ylimhalf = ceil(room.ywidth)*0.5*[-1 1];
end
plot(roomcorners(:,1),roomcorners(:,2),'k-', 'HandleVisibility', 'off', 'LineWidth', 2)
hold off;box on;grid on
xlim(xlimhalf),ylim(ylimhalf)
hh = legend(h, 'Location', 'southeast');
set(hh, 'Position', [0.2692 0.6767 0.1462 0.0950])
xlabel('x [m]'),ylabel('y [m]')
title('System geometry', 'FontSize', 15)
legendhitcallback
end