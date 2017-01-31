% load inputData01
% emdSTEM = PRISM01(atoms,cellDim);
% emdSTEM = PRISM02(emdSTEM);
load step2.mat
tic
emdSTEM_col = PRISM03_gpu(emdSTEM);
toc
image = sum(emdSTEM_col.stack(:,:,14:18),3);
figure, imagesc(image),axis image

