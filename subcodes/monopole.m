% function [m,s] = stat(x)
% n = length(x);
% m = sum(x)/n;
% s = sqrt(sum((x-m).^2/n));
% end

function srcPolar = monopole(f,rx,ry,XX,YY)
c = 343;
k = (2*pi*f)/c;

Rx = XX-rx;Ry = YY-ry;
R = sqrt(Rx.^2+Ry.^2)+eps;
% monopole: (1/r)*exp(jkr)

denomi = exp(1j*k*R);
srcPolar = (1./R).*denomi;
% srcPolar = (1./(R+0.05)).*denomi;
end