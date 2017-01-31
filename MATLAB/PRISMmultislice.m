function [emdSTEM] = PRISMmultislice(emdSTEM)
tic
% Mutlislice STEM simulation to compare with PRISM.

indFP = 1;  % Do a single frozen phonon each time
emdSTEM.MULTIprobeDefocusArray = 0;
emdSTEM.MULTIprobeSemiangleArray = 20/1000;  % Rads

% Probe positions
dxy = 0.25;
% xR = [0.6 0.64]*emdSTEM.cellDim(1);
% yR = [0.6 0.64]*emdSTEM.cellDim(2);
xR = [0.1 0.9]*emdSTEM.cellDim(1);
yR = [0.1 0.9]*emdSTEM.cellDim(2);
emdSTEM.MULTIxp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
emdSTEM.MULTIyp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
% emdSTEM.MULTIxp = 50;
% emdSTEM.MULTIyp = 50;

% Detector coords
dr = 2.5 / 1000;
alphaMax = emdSTEM.qMax * emdSTEM.lambda;
emdSTEM.MULTIdetectorAngles = (dr/2):dr:(alphaMax-dr/2);
q2 = emdSTEM.qxa.^2 + emdSTEM.qya.^2;
q1 = sqrt(q2);
alpha = q1 * emdSTEM.lambda;
alphaInds = round((alpha + dr/2) / dr);
% alphaInds(alphaInds<1) = 1;
alphaInds(alphaInds>length(emdSTEM.MULTIdetectorAngles)) = 0;
alphaMask = alphaInds > 0;
alphaIndsSub = alphaInds(alphaMask);
Ndet = length(emdSTEM.MULTIdetectorAngles);

% Initial probe
qProbeMax = emdSTEM.MULTIprobeSemiangleArray / emdSTEM.lambda;
dq = mean([emdSTEM.qxa(2,1) emdSTEM.qya(1,2)]);
PsiProbeInit = (erf((qProbeMax - q1) ...
    /(0.5*dq))*0.5 + 0.5);
PsiProbeInit(:) = PsiProbeInit ...
    .* exp((-1i*pi*emdSTEM.lambda ...
    * emdSTEM.MULTIprobeDefocusArray)*q2);
PsiProbeInit(:) = PsiProbeInit(:) ...
    / sqrt(sum(abs(PsiProbeInit(:)).^2));

% Transmission grating
trans = exp(1i*emdSTEM.sigma*emdSTEM.pot(:,:,:,indFP));

% Mask @ 1/4?
% xM = [(1:(emdSTEM.imageSize(1)/4)) ...
%     ((1-emdSTEM.imageSize(1)/4):0) + emdSTEM.imageSize(1)];
% yM = [(1:(emdSTEM.imageSize(2)/4)) ...
%     ((1-emdSTEM.imageSize(2)/4):0) + emdSTEM.imageSize(2)];
% xMask = false(emdSTEM.imageSize(1),1);
% yMask = false(1,emdSTEM.imageSize(2));
% xMask(xM) = true;
% yMask(yM) = true;

% Main simulation loop
emdSTEM.MULTIstack = zeros( ...
    length(emdSTEM.MULTIxp),...
    length(emdSTEM.MULTIyp),...
    length(emdSTEM.MULTIdetectorAngles));
psi = zeros(emdSTEM.imageSize);
% propMask = emdSTEM.prop(emdSTEM.qMask);
progressbar(0,2);
for a0 = 1:length(emdSTEM.MULTIxp)
    for a1 = 1:length(emdSTEM.MULTIyp)
%         % Shifted mask
%         dxy = round([emdSTEM.MULTIxp(a0) ...
%             emdSTEM.MULTIyp(a1)] ./ emdSTEM.pixelSize);
%         %         mask = circshift(maskReal,dxy);
%         %         inds0 = indsMask
%         xMaskShift = circshift(xMask,[dxy(1) 0]);
%         yMaskShift = circshift(yMask,[0 dxy(2)]);
        
        % Make probe
        psi(:) = PsiProbeInit ...
            .* exp(-2i*pi ...
            *(emdSTEM.qxa*emdSTEM.MULTIxp(a0) ...
            + emdSTEM.qya*emdSTEM.MULTIyp(a1)));
        
        % Propgate through all potential planes
        for a2 = 1:emdSTEM.numPlanes
%             psi = ifft2(psi);
%             psi(xMaskShift,yMaskShift) = ....
%                 psi(xMaskShift,yMaskShift) ...
%                 .* trans(xMaskShift,yMaskShift,a2);
%             psi = fft2(psi);
%             psi(emdSTEM.qMask) = psi(emdSTEM.qMask) .* propMask;
            % psi = ifft2(fft2(psi.*trans(:,:,a2)).*emdSTEM.prop);
            psi = fft2(ifft2(psi).*trans(:,:,a2)).*emdSTEM.prop;
        end
        
        % Record output
        emdSTEM.MULTIstack(a0,a1,:) = ...
            accumarray(alphaIndsSub,abs(psi(alphaMask)).^2,[Ndet 1]);
        
        % Completion
        comp = (a1 / length(emdSTEM.MULTIyp) ...
            + a0 - 1) / length(emdSTEM.MULTIxp);
        progressbar(comp,2);
    end
end



% figure(1)
% clf
% imagesc(abs(ifft2(psi)).^0.25)
% axis equal off
% colormap(gray(256))


toc
end