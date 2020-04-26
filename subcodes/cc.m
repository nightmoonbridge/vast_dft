function pts = cc(radius, numctrpts, initphase, rotdirection)
% pts = cc(1.5,20,0,'ccw');
% 
% figure
% scatter(pts(:,1),pts(:,2))
% axis equal
% box on
% xlim([-2.5 2.5]);ylim([-2.5 2.5])
% grid minor
% 

narginchk(2,4)
if nargin == 2
    initphase = 0;
    rotdirection = 'ccw';
elseif nargin == 3
    rotdirection = 'ccw';
end

if strcmpi(rotdirection, 'ccw')
    ridx = 1;
else
    ridx = -1;
end

pts = zeros(numctrpts,2);

dphi = (2*pi)/numctrpts;

midx = (0:numctrpts-1)';

exppos = radius*exp(1j*(initphase + ridx*midx*dphi));
pts(:,1) = real(exppos);
pts(:,2) = imag(exppos);

end