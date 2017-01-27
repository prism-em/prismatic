function [emdSTEM] = PRISM03(emdSTEM)
tic

% 03 - compute outputs from compact S matrix into 3D STEM outputs

% Inputs
emdSTEM.probeDefocusArray = 0;  % in Angstroms      % dim 4
emdSTEM.probeSemiangleArray = 25/1000;  % rads      % dim 5
emdSTEM.probeXtiltArray = 0/1000;  % rads           % dim 6
emdSTEM.probeYtiltArray = 0/1000;  % rads           % dim 7
% Probe positions
dxy = 0.25;
xR = [0.1 0.9]*emdSTEM.cellDim(1);
yR = [0.1 0.9]*emdSTEM.cellDim(2);
emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
% emdSTEM.xp = 50;
% emdSTEM.yp = 50;
% Detector output angles
dr = 5 / 1000;
alphaMax = emdSTEM.qMax * emdSTEM.lambda;
emdSTEM.detectorAngles = (dr/2):dr:(alphaMax-dr/2);

% Realspace coordinate system
r = emdSTEM.imageSize / emdSTEM.interpolationFactor / 2;
xVec = (-r(1)):(r(1)-1);
yVec = (-r(2)):(r(2)-1);
% Probe shift to center of cut out region
xCent = emdSTEM.cellDim(1) /emdSTEM.interpolationFactor / 2;
yCent = emdSTEM.cellDim(2) /emdSTEM.interpolationFactor / 2;

% Downsampled beams
emdSTEM.beamsReduce = emdSTEM.beams( ...
    1:emdSTEM.interpolationFactor:end,...
    1:emdSTEM.interpolationFactor:end);
imageSizeReduce = size(emdSTEM.beamsReduce);
xyBeams = zeros(length(emdSTEM.beamsIndex),2);
for a0 = 1:emdSTEM.numberBeams;
    [~,ind] = min(abs(emdSTEM.beamsReduce(:) - a0));
    [xx,yy] = ind2sub(imageSizeReduce,ind);
    xyBeams(a0,:) = [xx yy];
end

% Generate detector downsampled coordinate system
qxaReduce = emdSTEM.qxa( ...
    1:emdSTEM.interpolationFactor:end,...
    1:emdSTEM.interpolationFactor:end);
qyaReduce = emdSTEM.qya( ...
    1:emdSTEM.interpolationFactor:end,...
    1:emdSTEM.interpolationFactor:end);
Ndet = length(emdSTEM.detectorAngles);
emdSTEM.alphaInd = round((sqrt(qxaReduce.^2 + qyaReduce.^2) ...
    *emdSTEM.lambda - emdSTEM.detectorAngles(1)) / dr) + 1;
emdSTEM.alphaInd(emdSTEM.alphaInd<1) = 1;
alphaMask = emdSTEM.alphaInd <= Ndet;
alphaInds = emdSTEM.alphaInd(alphaMask);
emdSTEM.alphaInd(~alphaMask) = 0;

% Initialize pieces
emdSTEM.stackSize = [ ...
    length(emdSTEM.xp) ...
    length(emdSTEM.yp) ...
    length(emdSTEM.detectorAngles) ...
    length(emdSTEM.probeDefocusArray) ...
    length(emdSTEM.probeSemiangleArray) ...
    length(emdSTEM.probeXtiltArray) ...
    length(emdSTEM.probeYtiltArray)];
emdSTEM.stack = zeros(emdSTEM.stackSize,'single');
q1 = zeros(imageSizeReduce);
q2 = zeros(imageSizeReduce);
dq = mean([qxaReduce(2,1) qyaReduce(1,2)]);
PsiProbeInit = zeros(imageSizeReduce);
psi = zeros(imageSizeReduce);
intOutput = zeros(imageSizeReduce);

% Main loops


emdSTEM.stack = PRISM03_mexwrap(emdSTEM);


% scale = prod(emdSTEM.imageSize) / emdSTEM.interpolationFactor^2 / 4;
% for a0 = 1:length(emdSTEM.probeDefocusArray)
%     for a1 = 1:length(emdSTEM.probeSemiangleArray)
%         qProbeMax = emdSTEM.probeSemiangleArray(a1) / emdSTEM.lambda;
%         
%         for a2 = 1:length(emdSTEM.probeXtiltArray)
%             for a3 = 1:length(emdSTEM.probeYtiltArray)
%                 % Build unshifted, tilted probed probe
%                 q2(:) = ((qxaReduce  ...
%                     -  emdSTEM.probeXtiltArray(a2) / emdSTEM.lambda).^2 ...
%                     + (qyaReduce  ...
%                     -  emdSTEM.probeYtiltArray(a3) / emdSTEM.lambda).^2);
%                 q1(:) = sqrt(q2);
%                 PsiProbeInit(:) = min(max((qProbeMax - q1)/dq,0),1) ...
%                     .* exp((-1i * pi * emdSTEM.lambda ...
%                     * emdSTEM.probeDefocusArray) * q2);
%                 PsiProbeInit(:) = PsiProbeInit(:) ...
%                     / sqrt(sum(abs(PsiProbeInit(:)).^2));
%                 
%                 for ax = 1:length(emdSTEM.xp)
%                     for ay = 1:length(emdSTEM.yp)
%                         % Cut out the potential segment
%                         x0 = emdSTEM.xp(ax)/emdSTEM.pixelSize(1);
%                         y0 = emdSTEM.yp(ay)/emdSTEM.pixelSize(2);
%                         x = mod(xVec + round(x0) ,...
%                             emdSTEM.imageSize(1)) + 1;
%                         y = mod(yVec + round(y0),...
%                             emdSTEM.imageSize(2)) + 1;
%                         
%                         intOutput(:) = 0;
%                         for a5 = 1:emdSTEM.numFP
%                             % Build signal from beams
%                             psi(:) = 0;
%                             for a4 = 1:length(emdSTEM.beamsIndex)
%                                 xB = xyBeams(a4,1);
%                                 yB = xyBeams(a4,2);
%                                 
%                                 if abs(PsiProbeInit(xB,yB)) > 0
%                                     q0 = [qxaReduce(xB,yB) ...
%                                         qyaReduce(xB,yB)]; % ...
%                                     phaseShift = exp(-2i*pi ...
%                                         *(q0(1)*((x(1)-1)*emdSTEM.pixelSize(1)+xCent) ...
%                                         + q0(2)*((y(1)-1)*emdSTEM.pixelSize(1)+yCent)));
%                                     psi = psi ...
%                                         + (PsiProbeInit(xB,yB) * phaseShift) ...
%                                         * emdSTEM.Scompact(x,y,a4);
%                                 end
%                             end
%                             
%                             intOutput = intOutput ...
%                                 + abs(fft2(psi)).^2;
%                         end
%                                                 
%                         % Output signal
%                         emdSTEM.stack(ax,ay,:,a0,a1,a2,a3) = ...
%                             emdSTEM.stack(ax,ay,:,a0,a1,a2,a3) ...
%                             + reshape(accumarray(alphaInds,...
%                             intOutput(alphaMask),[Ndet 1]),...
%                             [1 1 Ndet]) * scale;
% %                         
% %                         comp = (((((ay / length(emdSTEM.yp) ...
% %                             + ax - 1) / length(emdSTEM.xp) ...
% %                             + a3 - 1) / length(emdSTEM.probeYtiltArray) ...
% %                             + a2 - 1) / length(emdSTEM.probeXtiltArray) ...
% %                             + a1 - 1) / length(emdSTEM.probeSemiangleArray) ...
% %                             + a0 - 1) / length(emdSTEM.probeDefocusArray);
% %                         progressbar(comp,2);
%                     end
%                     
%                     
%                         comp = ((((ax / length(emdSTEM.xp) ...
%                             + a3 - 1) / length(emdSTEM.probeYtiltArray) ...
%                             + a2 - 1) / length(emdSTEM.probeXtiltArray) ...
%                             + a1 - 1) / length(emdSTEM.probeSemiangleArray) ...
%                             + a0 - 1) / length(emdSTEM.probeDefocusArray);
%                         progressbar(comp,2);
%                 end
%             end
%         end
%     end
% end


% figure(1)
% clf
% % imagesc(sqrt([abs(ifft2(PsiProbe)) ...
% %     abs(psi)*emdSTEM.interpolationFactor^2;
% %     abs(psi)*emdSTEM.interpolationFactor^2 ...
% %     abs(ifft2(PsiProbe))]))
% % imagesc(abs(ifft2(PsiProbe)))
% imagesc(fftshift(abs(fft2(psi))))
% % imagesc(sqrt(sum(emdSTEM.pot(x,y,:,1),3)));
% % imagesc(abs(PsiProbe))
% axis equal off
% colormap(gray(256))
% set(gca,'position',[0 0 1 1])


toc
end