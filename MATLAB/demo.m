
load('step2')
emdSTEM = PRISM03_ajp(emdSTEM);
% emdSTEM = PRISM03(emdSTEM);
image = sum(emdSTEM.stack(:,:,14:18),3);
figure, imagesc(image),axis image

