function [emdSTEM] = PRISM02(emdSTEM)
tic
% 02 - calculate the S matrix from the projected potentials, out to some
% maximum scattering angle for the input tensor.

% Inputs
emdSTEM.E0 = 80e3;  % Microscope voltage in volts
emdSTEM.alphaBeamMax = 32/1000;%80/1000;  % in rads, for specifying maximum angle

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
qx = makeFourierCoords(emdSTEM.imageSize(1),emdSTEM.pixelSize(1));
qy = makeFourierCoords(emdSTEM.imageSize(2),emdSTEM.pixelSize(2));
[emdSTEM.qya,emdSTEM.qxa] = meshgrid(qy,qx);
q2 = emdSTEM.qxa.^2 + emdSTEM.qya.^2;

% propagators and mask
emdSTEM.prop = exp(-1i*pi*emdSTEM.lambda*q2*emdSTEM.sliceThickness);
% emdSTEM.propBack = exp(-1i*pi*emdSTEM.lambda*q2 ...
%     *(emdSTEM.sliceThickness*-size(emdSTEM.pot,3)));
emdSTEM.qMax = min(max(abs(qx)),max(abs(qy)))/2;
% qMask = (abs(emdSTEM.qxa) <= emdSTEM.qMax) ...
%     & (abs(emdSTEM.qya) <= emdSTEM.qMax);
qMask = false(emdSTEM.imageSize);
qMask([(1:(emdSTEM.imageSize(1)/4)) ...
    ((1-emdSTEM.imageSize(1)/4):0)+emdSTEM.imageSize(1)],...
    [(1:(emdSTEM.imageSize(2)/4)) ...
    ((1-emdSTEM.imageSize(2)/4):0)+emdSTEM.imageSize(2)]) = true;
emdSTEM.qMask = qMask;

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
[xBeams,yBeams] = ind2sub(emdSTEM.imageSize,emdSTEM.beamsIndex);

% Create bilinear interpolation kernel
v = 1 - abs((1-emdSTEM.interpolationFactor) ...
    :(emdSTEM.interpolationFactor-1)) ...
    / emdSTEM.interpolationFactor;
emdSTEM.kernel = v'*v;

% Generate compact S matrix
emdSTEM.Scompact = zeros(...
    emdSTEM.imageSize(1),...
    emdSTEM.imageSize(2),...
    emdSTEM.numberBeams,...
    emdSTEM.numFP,'single');
psi = zeros(emdSTEM.imageSize);
trans = exp(1i*emdSTEM.sigma*emdSTEM.pot);
progressbar(0,2);
for a0 = 1:emdSTEM.numberBeams
    % Make shifted antialising apertures
    shift = [xBeams(a0) yBeams(a0)]-1;
    qMaskShift = circshift(qMask,shift);
    prop = emdSTEM.prop .* qMaskShift;
    %     propBack = emdSTEM.propBack .* qMaskShift;
    
    for a1 = 1:emdSTEM.numFP
        % Initialize plane wave
        psi(:) = 0;
        psi(emdSTEM.beamsIndex(a0)) = 1;
        psi = ifft2(psi);
        
        % Propgate through all potential planes
        for a2 = 1:emdSTEM.numPlanes
            psi = ifft2(fft2(psi.*trans(:,:,a2,a1)).*prop);
        end
        
        % Shift S matrix entry to excitation origin and output result
        emdSTEM.Scompact(:,:,a0,a1) = psi;
        %         psi = circshift(psi,-shift);
        %         emdSTEM.Scompact(:,:,a0,a1) = ...
        %             reshape(psi(qMask),emdSTEM.imageSize/2);
        
        comp = (a1 / emdSTEM.numFP ...
            + a0 - 1) / emdSTEM.numberBeams;
        progressbar(comp,2);
    end
end
if comp < 1
    progressbar(1,2);
end


% % % % log(prod(emdSTEM.imageSize)*sum(mask(:))*4*2)/log(10)
% % prod(emdSTEM.imageSize)*sum(mask(:))*4*2 /1024^3
% figure(1)
% clf
% imagesc(abs(psi))
% % imagesc(fftshift(qMaskShift+psi))
% % imagesc(angle(psi))
% % imagesc(conv2(fftshift(emdSTEM.beams),ones(5),'same'))
% % Ip = emdSTEM.kernel ...
% %     + circshift(emdSTEM.kernel,[emdSTEM.interpolationFactor 0]) ...
% %     + circshift(emdSTEM.kernel,[0 emdSTEM.interpolationFactor]) ...
% %     + circshift(emdSTEM.kernel,emdSTEM.interpolationFactor+[0 0]) ;
% % imagesc(fftshift(mask)
% axis equal off 
% colormap(gray(256))
% colorbar
% set(gca,'position',[0 0 1 1])
% 
% % emdSTEM.alphaMax = emdSTEM.qMax * emdSTEM.lambda;
% % emdSTEM.alphaMax * 1000
% % qMax = emdSTEM.alphaMax / emdSTEM.lambda;
% 
% 
% 
% % sum(qMask(:))
% % sum(qMask(:))^2*4/1024^3
% % sum(qMask(:))*prod(emdSTEM.imageSize)/4/1024^3*8/64
% 
% 
% 
% % Smap = z



toc
end