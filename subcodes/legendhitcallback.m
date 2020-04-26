function legendhitcallback
calllegend = legend;
calllegend.ItemHitFcn = @hitcallback;
end

function hitcallback(src,evnt)
% function hitcallback(src,evnt)
% Check MatWorks website
% https://se.mathworks.com/help/releases/R2017a/matlab/creating_plots/create-callbacks-for-interacting-with-legend-items.html

switch evnt.Region
    case 'icon'
        lwidth = 2;
        if evnt.Peer.LineWidth == lwidth
            evnt.Peer.LineWidth = 1;
        else 
            evnt.Peer.LineWidth = lwidth;
        end
    case 'label'
        evnt.Peer.LineWidth = 1;
        if strcmp(evnt.Peer.Visible,'on')
            evnt.Peer.Visible = 'off';
        else
            evnt.Peer.Visible = 'on';
        end
end

end