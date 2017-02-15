load inputData01
emdSTEM = PRISM01(atoms,cellDim);
emdSTEM = PRISM02_ajp(emdSTEM);
tic
emdSTEM_tmp = PRISM03_ajp(emdSTEM);
toc

tic
emdSTEM = PRISM03(emdSTEM);
toc

% image = sum(emdSTEM_col.stack(:,:,14:18),3);
% figure, imagesc(image),axis image

