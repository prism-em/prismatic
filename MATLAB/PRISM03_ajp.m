function [emdSTEM,psi] = PRISM03_ajp(emdSTEM)
tic

% 03 - compute outputs from compact S matrix into 3D STEM outputs

% Inputs
% emdSTEM.probeDefocusArray = [-80 -40 0];%[-40 -30 -20 -10 0];%[-40 -20 0];%[-100 -50 0];  % in Angstroms      % dim 4
% emdSTEM.probeSemiangleArray = [70 50 30 10]/1000;%fliplr([25 50 75])/1000;%25/1000;  % rads      % dim 5
emdSTEM.probeDefocusArray = 0;%[-40 0 40];%[-80 -60 -40 -20 0 20];%[-40 -30 -20 -10 0];%[-40 -20 0];%[-100 -50 0];  % in Angstroms      % dim 4
emdSTEM.probeSemiangleArray = 20/1000;%[30 20 10]/1000;%fliplr([25 50 75])/1000;%25/1000;  % rads      % dim 5
emdSTEM.probeXtiltArray = 0/1000;  % rads           % dim 6
emdSTEM.probeYtiltArray = 0/1000;  % rads           % dim 7
% Probe positions
dxy = 0.25 * 2;
% xR = [0.3 0.5]*emdSTEM.cellDim(1);
% yR = [0.65 0.85]*emdSTEM.cellDim(2);
% xR = [0.40 0.45]*emdSTEM.cellDim(1);
% yR = [0.70 0.75]*emdSTEM.cellDim(2);
xR = [0.1 0.9]*emdSTEM.cellDim(1);
yR = [0.1 0.9]*emdSTEM.cellDim(2);
emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
% emdSTEM.xp = 50;
% emdSTEM.yp = 50;
% Detector output angles
dr = 2.5 / 1000;
alphaMax = emdSTEM.qMax * emdSTEM.lambda;
emdSTEM.detectorAngles = (dr/2):dr:(alphaMax-dr/2);
flag_plot = 0;
flag_keep_beams = 0;

% Realspace coordinate system
r = emdSTEM.imageSizeOutput / emdSTEM.interpolationFactor / 2;
xVec = ((-r(1)):(r(1)-1));
yVec = ((-r(2)):(r(2)-1));

% Downsampled beams
emdSTEM.beamsReduce = emdSTEM.beamsOutput( ...
    1:(emdSTEM.interpolationFactor):end,...
    1:(emdSTEM.interpolationFactor):end);
imageSizeReduce = size(emdSTEM.beamsReduce);
xyBeams = zeros(length(emdSTEM.beamsIndex),2);
for a0 = 1:emdSTEM.numberBeams;
    [~,ind] = min(abs(emdSTEM.beamsReduce(:) - a0));
    [xx,yy] = ind2sub(imageSizeReduce,ind);
    xyBeams(a0,:) = [xx yy];
end

% Generate detector downsampled coordinate system
qxaReduce = emdSTEM.qxaOutput( ...
    1:emdSTEM.interpolationFactor:end,...
    1:emdSTEM.interpolationFactor:end);
qyaReduce = emdSTEM.qyaOutput( ...
    1:emdSTEM.interpolationFactor:end,...
    1:emdSTEM.interpolationFactor:end);
Ndet = length(emdSTEM.detectorAngles);

% figure(56)
% clf
% imagesc(qxaReduce)
% axis equal off

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
scale = emdSTEM.interpolationFactor^4;
if flag_keep_beams == 1
    emdSTEM.beamsOutput = zeros( ...
        imageSizeReduce(1),...
        imageSizeReduce(2),...
        length(emdSTEM.beamsIndex));
end
% 
% emdSTEM.stack = loop_wrapper(emdSTEM.Scompact, ...
%     emdSTEM.stack, emdSTEM.probeDefocusArray, emdSTEM.probeSemiangleArray, ...
%     emdSTEM.probeXtiltArray, emdSTEM.probeYtiltArray,qxaReduce, qyaReduce, ...
%     emdSTEM.xp, emdSTEM.yp, emdSTEM.beamsIndex,xyBeams,xVec, yVec, imageSizeReduce,emdSTEM.imageSizeOutput, ...
%     emdSTEM.detectorAngles,emdSTEM.cellDim,emdSTEM.pixelSizeOutput, scale, ...
%     emdSTEM.lambda, dr, dq, Ndet, ...
%      emdSTEM.numFP);
emdSTEM.Scompact = double(emdSTEM.Scompact);
emdSTEM.stack = double(emdSTEM.stack);
tic
b = PRISM03_mex(emdSTEM.Scompact, ...
    emdSTEM.stack, emdSTEM.probeDefocusArray, emdSTEM.probeSemiangleArray, ...
    emdSTEM.probeXtiltArray, emdSTEM.probeYtiltArray,qxaReduce, qyaReduce, ...
    emdSTEM.xp, emdSTEM.yp, emdSTEM.beamsIndex,xyBeams,xVec, yVec, imageSizeReduce,emdSTEM.imageSizeOutput, ...
    emdSTEM.detectorAngles,emdSTEM.cellDim,emdSTEM.pixelSizeOutput, scale, ...
    emdSTEM.lambda, dr, dq, Ndet, ...
     emdSTEM.numFP);
 toc
 emdSTEM.stack = loop_wrapper(emdSTEM.Scompact, ...
    emdSTEM.stack, emdSTEM.probeDefocusArray, emdSTEM.probeSemiangleArray, ...
    emdSTEM.probeXtiltArray, emdSTEM.probeYtiltArray,qxaReduce, qyaReduce, ...
    emdSTEM.xp, emdSTEM.yp, emdSTEM.beamsIndex,xyBeams,xVec, yVec, imageSizeReduce,emdSTEM.imageSizeOutput, ...
    emdSTEM.detectorAngles,emdSTEM.cellDim,emdSTEM.pixelSizeOutput, scale, ...
    emdSTEM.lambda, dr, dq, Ndet, ...
     emdSTEM.numFP);

 
 
% for a0 = 1:length(emdSTEM.probeDefocusArray)
%     for a1 = 1:length(emdSTEM.probeSemiangleArray)
%         qProbeMax = emdSTEM.probeSemiangleArray(a1) / emdSTEM.lambda;
%         
%         for a2 = 1:length(emdSTEM.probeXtiltArray)
%             for a3 = 1:length(emdSTEM.probeYtiltArray)
%                 qxaShift = qxaReduce  ...
%                     -  emdSTEM.probeXtiltArray(a2) / emdSTEM.lambda;
%                 qyaShift = qyaReduce  ...
%                     -  emdSTEM.probeYtiltArray(a3) / emdSTEM.lambda;
%                 q2(:) = (qxaShift.^2 + qyaShift.^2);
%                 q1(:) = sqrt(q2);
%                 
%                 % Build shifted detector coordinate array
%                 emdSTEM.alphaInd = round((q1*emdSTEM.lambda ...
%                     - emdSTEM.detectorAngles(1)) / dr) + 1;
%                 emdSTEM.alphaInd(emdSTEM.alphaInd<1) = 1;
%                 alphaMask = emdSTEM.alphaInd <= Ndet;
%                 alphaInds = emdSTEM.alphaInd(alphaMask);
%                 emdSTEM.alphaInd(~alphaMask) = 0;
%                 
%                 % Build unshifted, tilted probed probe with defocus
%                 PsiProbeInit(:) = (erf((qProbeMax - q1) ...
%                     /(0.5*dq))*0.5 + 0.5);
%                 PsiProbeInit(:) = PsiProbeInit ...
%                     .* exp((-1i*pi*emdSTEM.lambda ...
%                     * emdSTEM.probeDefocusArray(a0))*q2);
%                 PsiProbeInit(:) = PsiProbeInit(:) ...
%                     / sqrt(sum(abs(PsiProbeInit(:)).^2));
%              
%                 % Calculate additional probe shift in cut out region to
%                 % beam tilt, apply in final probe calculation.
%                 zTotal = emdSTEM.cellDim(3);  % plus defocus shift?
%                 xTiltShift = -zTotal ...
%                     * tan(emdSTEM.probeXtiltArray(a2));
%                 yTiltShift = -zTotal ...
%                     * tan(emdSTEM.probeYtiltArray(a3));
%                 
%                 for ax = 1:length(emdSTEM.xp)
%                     for ay = 1:length(emdSTEM.yp)
%                         % Cut out the potential segment
%                         x0 = emdSTEM.xp(ax)/emdSTEM.pixelSizeOutput(1);
%                         y0 = emdSTEM.yp(ay)/emdSTEM.pixelSizeOutput(2);
%                         x = mod(xVec + round(x0),...
%                             emdSTEM.imageSizeOutput(1)) + 1;
%                         y = mod(yVec + round(y0),...
%                             emdSTEM.imageSizeOutput(2)) + 1;
%                         
%                         intOutput(:) = 0;
%                         for a5 = 1:emdSTEM.numFP
%                             % Build signal from beams
%                             psi(:) = 0;
%                             for a4 = 1:length(emdSTEM.beamsIndex)
%                                 xB = xyBeams(a4,1);
%                                 yB = xyBeams(a4,2);
%                                 if abs(PsiProbeInit(xB,yB)) > 0
%                                     q0 = [qxaReduce(xB,yB) ...
%                                         qyaReduce(xB,yB)];
%                                     phaseShift = exp(-2i*pi ...
%                                         *(q0(1)*(emdSTEM.xp(ax) + xTiltShift) ...
%                                         + q0(2)*(emdSTEM.yp(ay) + yTiltShift)));
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
%                     end
%                     
%                     comp = ((((ax / length(emdSTEM.xp) ...
%                         + a3 - 1) / length(emdSTEM.probeYtiltArray) ...
%                         + a2 - 1) / length(emdSTEM.probeXtiltArray) ...
%                         + a1 - 1) / length(emdSTEM.probeSemiangleArray) ...
%                         + a0 - 1) / length(emdSTEM.probeDefocusArray);
%                     progressbar(comp,2);
%                 end
%             end
%         end
%     end
% end


if flag_plot == true
    figure(1)
    clf
    imagesc(fftshift(abs(fft2(psi)).^0.5))
    axis equal off
    colormap(jet)
    set(gca,'position',[0 0 1 1])
    
    figure(2)
    clf
    amp = abs(psi).^0.5;
    ampRange = [0 1];    
    phase = angle(psi)+pi*0.15;
    [Irgb] = colorComplex(amp/max(amp(:)),phase,ampRange);
    imagesc(Irgb)
    axis equal off
    colormap(jet)
    set(gca,'position',[0 0 1 1])
end

toc
end

function stack = loop_wrapper(Scompact,stack,probeDefocusArray, probeSemiangleArray, ...
    probeXtiltArray, probeYtiltArray,qxaReduce, qyaReduce, ...
    xp, yp, beamsIndex, xyBeams,xVec, yVec, imageSizeReduce, imageSizeOutput, ...
    detectorAngles,cellDim, pixelSizeOutput, scale, ...
    lambda, dr, dq, Ndet, numFP)

% the point of this function is just to transcribe the above loop into a
% format that is closer to what will actually be written in C++ for easier
% transcription. Mostly the parsing of the emdSTEM object is the difference

PsiProbeInit = zeros(imageSizeReduce);
q1 = zeros(imageSizeReduce);
q2 = zeros(imageSizeReduce);

for a0 = 1:length(probeDefocusArray)
    for a1 = 1:length(probeSemiangleArray)
        qProbeMax = probeSemiangleArray(a1) / lambda
        
        for a2 = 1:length(probeXtiltArray)
            for a3 = 1:length(probeYtiltArray)
                qxaShift = qxaReduce  ...
                    -  probeXtiltArray(a2) / lambda;
                qyaShift = qyaReduce  ...
                    -  probeYtiltArray(a3) / lambda;
                q2(:) = (qxaShift.^2 + qyaShift.^2);
                q1(:) = sqrt(q2);
                
                % Build shifted detector coordinate array
                alphaInd = round((q1*lambda ...
                    - detectorAngles(1)) / dr) + 1;
                alphaInd(alphaInd<1) = 1;
                alphaMask = alphaInd <= Ndet;
                alphaInds = alphaInd(alphaMask);
%                 alphaInd(~alphaMask) = 0;
                
                % Build unshifted, tilted probed probe with defocus
                PsiProbeInit(:) = (erf((qProbeMax - q1) ...
                    /(0.5*dq))*0.5 + 0.5);
                PsiProbeInit(:) = PsiProbeInit ...
                    .* exp((-1i*pi*lambda ...
                    * probeDefocusArray(a0))*q2);
                PsiProbeInit(:) = PsiProbeInit(:) ...
                    / sqrt(sum(abs(PsiProbeInit(:)).^2));
             
                % Calculate additional probe shift in cut out region to
                % beam tilt, apply in final probe calculation.
                
                %%% left off here
                zTotal = cellDim(3);  % plus defocus shift?
                xTiltShift = -zTotal ...
                    * tan(probeXtiltArray(a2));
                yTiltShift = -zTotal ...
                    * tan(probeYtiltArray(a3));
                
                for ax = 1:length(xp)
                    for ay = 1:length(yp)
                        % Cut out the potential segment
                        x0 = xp(ax)/pixelSizeOutput(1);
                        y0 = yp(ay)/pixelSizeOutput(2);
                        x = mod(xVec + round(x0),...
                            imageSizeOutput(1)) + 1;
                        y = mod(yVec + round(y0),...
                            imageSizeOutput(2)) + 1;
                        
                        intOutput(:) = 0;
                        for a5 = 1:numFP
                            % Build signal from beams
                            psi(:) = 0;
                            for a4 = 1:length(beamsIndex)
                                xB = xyBeams(a4,1);
                                yB = xyBeams(a4,2);
                                if abs(PsiProbeInit(xB,yB)) > 0
                                    q0 = [qxaReduce(xB,yB) ...
                                        qyaReduce(xB,yB)];
                                    phaseShift = exp(-2i*pi ...
                                        *(q0(1)*(xp(ax) + xTiltShift) ...
                                        + q0(2)*(yp(ay) + yTiltShift)));
                                    psi = psi ...
                                        + (PsiProbeInit(xB,yB) * phaseShift) ...
                                        * Scompact(x,y,a4);
%                                             psi = psi ...
%                                         + (PsiProbeInit(xB,yB) * phaseShift) ...
%                                         * T;
                                end
                            end
                            
                            intOutput = intOutput ...
                                + abs(fft2(psi)).^2;
                        end
                        
                        % Output signal
                        stack(ax,ay,:,a0,a1,a2,a3) = ...
                            stack(ax,ay,:,a0,a1,a2,a3) ...
                            + reshape(accumarray(alphaInds,...
                            intOutput(alphaMask),[Ndet 1]),...
                            [1 1 Ndet]) * scale;
                    end
                    
                    comp = ((((ax / length(xp) ...
                        + a3 - 1) / length(probeYtiltArray) ...
                        + a2 - 1) / length(probeXtiltArray) ...
                        + a1 - 1) / length(probeSemiangleArray) ...
                        + a0 - 1) / length(probeDefocusArray);
                    progressbar(comp,2);
                end
            end
        end
    end
end

end