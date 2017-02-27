function [emdSTEM] = PRISM02(emdSTEM)
% tic
% 02 - calculate the S matrix from the projected potentials, out to some
% maximum scattering angle for the input tensor.

% Inputs
emdSTEM.E0 = 80e3;  % Microscope voltage in volts
emdSTEM.alphaBeamMax = 24/1000; % in rads, for specifying maximum angle

% Calculate wavelength and electron interaction parameter
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
emdSTEM.lambda = h/sqrt(2*m*e*emdSTEM.E0) ...
    /sqrt(1 + e*emdSTEM.E0/2/m/c^2) * 10^10; % wavelength in A
emdSTEM.sigma = (2*pi/emdSTEM.lambda/emdSTEM.E0) ...
    *(m*c^2+e*emdSTEM.E0)/(2*m*c^2+e*emdSTEM.E0);

% Generate Fourier coordinates
emdSTEM.imageSize = [size(emdSTEM.pot,1) size(emdSTEM.pot,2)];  % reset image size
qx = makeFourierCoords(emdSTEM.imageSize(1),emdSTEM.pixelSize(1));
qy = makeFourierCoords(emdSTEM.imageSize(2),emdSTEM.pixelSize(2));
[emdSTEM.qya,emdSTEM.qxa] = meshgrid(qy,qx);
q2 = emdSTEM.qxa.^2 + emdSTEM.qya.^2;

% propagators and mask
emdSTEM.qMax = min(max(abs(qx)),max(abs(qy)))/2;
qMask = false(emdSTEM.imageSize);
qMask([(1:(emdSTEM.imageSize(1)/4)) ...
    ((1-emdSTEM.imageSize(1)/4):0)+emdSTEM.imageSize(1)],...
    [(1:(emdSTEM.imageSize(2)/4)) ...
    ((1-emdSTEM.imageSize(2)/4):0)+emdSTEM.imageSize(2)]) = true;
emdSTEM.qMask = qMask;
emdSTEM.prop = qMask ...
    .* exp((-1i*pi*emdSTEM.lambda*emdSTEM.sliceThickness)*q2);
emdSTEM.propBack = qMask ...
    .* exp((1i*pi*emdSTEM.lambda*emdSTEM.cellDim(3))*q2);


% Generate S-matrix map 
xv = makeFourierCoords(emdSTEM.imageSize(1),1/emdSTEM.imageSize(1));
yv = makeFourierCoords(emdSTEM.imageSize(2),1/emdSTEM.imageSize(2));
[ya,xa] = meshgrid(round(yv),round(xv));
mask = qMask ...
    & (q2 < (emdSTEM.alphaBeamMax/emdSTEM.lambda)^2) ...
    & (mod(xa,emdSTEM.interpolationFactor) == 0) ...
    & (mod(ya,emdSTEM.interpolationFactor) == 0);
emdSTEM.numberBeams = sum(mask(:));
emdSTEM.beams = zeros(emdSTEM.imageSize);
emdSTEM.beams(mask==true) = 1:emdSTEM.numberBeams;
emdSTEM.beamsIndex = find(emdSTEM.beams);

% Subset of non-zero pixels to export in compact S-matrix
emdSTEM.qxInd = [(1:(emdSTEM.imageSize(1)/4)) ...
    ((1-emdSTEM.imageSize(1)/4):0)+emdSTEM.imageSize(1)];
emdSTEM.qyInd = [(1:(emdSTEM.imageSize(2)/4)) ...
    ((1-emdSTEM.imageSize(2)/4):0)+emdSTEM.imageSize(2)];

% Generate compact S matrix
emdSTEM.Scompact = zeros(...
    emdSTEM.imageSize(1)/2,...
    emdSTEM.imageSize(2)/2,...
    emdSTEM.numberBeams,...
    emdSTEM.numFP,'single');
psi = zeros(emdSTEM.imageSize);
trans = exp(1i*emdSTEM.sigma*emdSTEM.pot);
comp = 0;
progressbar(comp,2);
for a0 = 1:emdSTEM.numberBeams
    for a1 = 1:emdSTEM.numFP
        % Initialize plane wave
        psi(:) = 0;
        psi(emdSTEM.beamsIndex(a0)) = 1;
        psi = ifft2(psi);
        
        % Propgate through all potential planes
        for a2 = 1:emdSTEM.numPlanes
            psi = ifft2(fft2(psi.*trans(:,:,a2,a1)).*emdSTEM.prop);
        end
%         psi(:) = fft2(psi).*emdSTEM.propBack;  % Stay in Fourier space
        %         psi(:) = ifft2(fft2(psi).*emdSTEM.propBack);
        
        % Output subset of S-matrix
        psi(:) = fft2(psi);
        emdSTEM.Scompact(:,:,a0,a1) = ...
            ifft2(psi(emdSTEM.qxInd,emdSTEM.qyInd));
        
        comp = (a1 / emdSTEM.numFP ...
            + a0 - 1) / emdSTEM.numberBeams;
        progressbar(comp,2);
    end
end
if comp < 1
    progressbar(1,2);
end

% downsample all Fourier components by x2 to match output
emdSTEM.imageSizeOutput = emdSTEM.imageSize / 2;
emdSTEM.pixelSizeOutput = emdSTEM.pixelSize * 2;
emdSTEM.qxaOutput = emdSTEM.qxa(emdSTEM.qxInd,emdSTEM.qyInd);
emdSTEM.qyaOutput = emdSTEM.qya(emdSTEM.qxInd,emdSTEM.qyInd);
% emdSTEM.qMax = emdSTEM.qMax / 2;
% emdSTEM.qMask = qMask(emdSTEM.qxInd,emdSTEM.qyInd);
% emdSTEM.prop = emdSTEM.prop(emdSTEM.qxInd,emdSTEM.qyInd);
emdSTEM.beamsOutput = emdSTEM.beams(emdSTEM.qxInd,emdSTEM.qyInd);



% 
% toc
end