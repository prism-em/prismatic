function [pot] = projPot(atomID,xr,yr)
% Projected potential function - potentials from Kirkland
% Version 02 - separate function
% potMin is now gone, instead set from outside edge

ss = 8;  % Super sampling for potential integration (should be even!!)

% Get fparams file
load('fparams.mat');

% Constants
a0 = 0.5292;
e = 14.4;
% term1 = 2*pi^2*a0*e;
% term2 = 2*pi^(5/2)*a0*e;
term1 = 4*pi^2*a0*e;
term2 = 2*pi^2*a0*e;

% Make supersampled 2D grid for integration
dx = (xr(2) - xr(1));
dy = (yr(2) - yr(1));
sub = (-(ss-1)/ss/2):(1/ss):((ss-1)/ss/2);
[x1,x2] = meshgrid(xr,sub*dx);
xv = x1(:) + x2(:);
[y1,y2] = meshgrid(yr,sub*dy);
yv = y1(:) + y2(:);
[ya,xa] = meshgrid(yv,xv);
r2 = xa.^2 + ya.^2;
r = sqrt(r2);

% Compute potential
ap = fparams(atomID,:);
potSS = term1*( ...
    ap(1)*besselk(0,2*pi*sqrt(ap(2))*r) ...
    + ap(3)*besselk(0,2*pi*sqrt(ap(4))*r) ...
    + ap(5)*besselk(0,2*pi*sqrt(ap(6))*r)) ...
    + term2*( ...
    ap(7)/ap(8)*exp(-pi^2/ap(8)*r2) ...
    + ap(9)/ap(10)*exp(-pi^2/ap(10)*r2) ...
    + ap(11)/ap(12)*exp(-pi^2/ap(12)*r2));
% Integrate!
potMid = zeros(length(xr),length(yr)*ss);
for a0 = 1:ss
    potMid = potMid + potSS(a0:ss:(end+a0-ss),:);
end
pot = zeros(length(xr),length(yr));
for a0 = 1:ss
    pot = pot + potMid(:,a0:ss:(end+a0-ss));
end
pot = pot / ss^2;
% Cut off radius for potential (to remove edge discontinuities)
% pot = pot - potMin;
% pot(pot<0) = 0;
% Compute potMin from edge
[~,xInd] = min(abs(xr));
[~,yInd] = min(abs(yr));
dx = round(sqrt(2*xInd-1));
dy = round(sqrt(2*yInd-1));
xv = [xInd-dx xInd+dx xInd-dx xInd+dx 1 1 length(xr) length(xr)];
yv = [1 1 length(yr) length(yr) yInd-dy yInd+dy yInd-dy yInd+dy];
potMin = max(pot(sub2ind(size(pot),xv,yv)));
pot = pot - potMin;
pot(pot<0) = 0;

end