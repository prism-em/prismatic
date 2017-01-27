function stack = PRISM03_mexwrap(emdSTEM)
% This wrapper takes in an emdSTEM object and separates out the relevant
% data to be passed to the mex function, which is itself a wrapper around
% the PRISM C++/CUDA backend.

% MATLAB will only shallow copy these so these next steps are inexpensive


% scalars
imageSizeReduce     = size(emdSTEM.beamsReduce);
imageSize           = emdSTEM.imageSize;
scale               = prod(emdSTEM.imageSize) / emdSTEM.interpolationFactor^2 / 4;
interpolationFactor = emdSTEM.interpolationFactor;
pixelSize           = emdSTEM.pixelSize;
lambda              = emdSTEM.lambda;

% arrays
probeDefocusArray   = emdSTEM.probeDefocusArray;
probeSemiangleArray = emdSTEM.probeSemiangleArray;
probeXtiltArray     = emdSTEM.probeXtiltArray;
probeYtiltArray     = emdSTEM.probeYtiltArray;

qxaReduce           = emdSTEM.qxa(1:interpolationFactor:end, ...
                                  1:interpolationFactor:end);
                              
qyaReduce           = emdSTEM.qya(1:interpolationFactor:end, ...
                                  1:interpolationFactor:end);
xp = emdSTEM.xp;
yp = emdSTEM.yp;



stack = PRISM03_mex(probeDefocusArray, probeSemiangleArray,probeXtiltArray, ...
                    probeYtiltArray, qxaReduce, qyaReduce, xp, yp, ...
                    imageSize, imageSizeReduce, scale, interpolationFactor, ...
                    pixelSize, lambda );
end