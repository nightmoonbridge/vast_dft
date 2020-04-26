function mpos = dcirc(M,radi,rado,cent,rot)
% function mpos = dcirc(M,radi,rado,cent,rot)
% Double circle microphone positions
% M        : the number of microphone positions, default is 6
% radi     : the radius of the inner circle, default is 0.2
% rado     : the radius of the outer circle, default is 0.5
% cent     : the center of two circles, concentric, default is [0 0]
% rot      : rotation angle in radian, default is 0
if nargin == 0
    M = 6;
    radi = 0.2;
    rado = 0.5;
    cent = [0 0];
    rot = 0;
elseif nargin == 1
    radi = 0.2;
    rado = 0.5;
    cent = [0 0];
    rot = 0;
elseif nargin == 2
    rado = radi+0.2;
    cent = [0 0];
    rot = 0;
elseif nargin == 3
    cent = [0 0];
    rot = 0;
elseif nargin == 4
    rot = 0;
elseif nargin == 5
else
    error('See help')
end

% radi = 0.2;
mindx = 0:M/2-1;
rtheta = radi.*exp(1j*mindx'*2*pi/(M/2)+1j*rot);

% rado = 0.5;
rtheta2 = rado.*exp(1j*(mindx*2+1)'*pi/(M/2)+1j*rot);

mpos = zeros(M,2);
mpos(:,1) = cent(1)+[real(rtheta);real(rtheta2)];
mpos(:,2) = cent(2)+[imag(rtheta);imag(rtheta2)];
end

%{
circ1 = radius.*exp(1j*(0:5:360)'*2*pi/360);
circ2 = radius2.*exp(1j*(0:5:360)'*2*pi/360);


figure
plot(real(circ1),imag(circ1))
hold on
plot(real(circ2),imag(circ2))
scatter(mpos1x,mpos1y)
scatter(mpos2x,mpos2y);axis equal
xlim([-1 1]);ylim([-1 1])
grid minor
%}
