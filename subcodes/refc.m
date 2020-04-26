function [betaavg, alphaavg, minRT60] = refc(c,roomdim,rt60)
% function [betaavg, alphaavg, minRT60] = REFC(c,roomdim,rt60)
%
% This function returns an averaged reflection coefficient from the given
% reverberation time (rt60)
% Input variables
% - c         : the speed of sound            [m/s]
% - roomdim   : the room dimension (WxLxH) in [meters]
% - rt60      : the reverberation time in     [seconds]
%
% Output variables
% - betaavg   : an averaged reflection coefficient
% - alphaavg  : an averaged absoprtion coefficient
% - minRT60   : the minimum RT60 based on the given information
%
% The details can be found from [+]
%
% [+] A. Pierce, Acoustics, Acoustical Society of America, 2019, p. 296

V = prod(roomdim);
S = 2*(prod(roomdim(1:2)) + prod(roomdim(2:3)) + prod(roomdim(1:2:3)));
minRT60 = round((24*log(10)*V)/(c*S),2);
alphaavg = (24*log(10)*V)/(c*S*rt60);
% alphaavg = minRT60/rt60;
betaavg = sqrt(1 - alphaavg);
end