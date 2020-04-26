function ccpos = drawCircle(radius, center)
% This returns a circle based on 1) radius and 2) center information.
narginchk(1,2)
if nargin == 1
    center = zeros(1,2);
end
if numel(center) ~= 2
    error('The center should be a 1x2 row vector.')
elseif iscolumn(center)
    center = center';
end

circtheta = linspace(0,2*pi)';
circpos = radius*exp(1j*circtheta);
circx = real(circpos);
circy = imag(circpos);
ccpos = center + [circx,circy];
end