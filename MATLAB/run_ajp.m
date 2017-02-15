load inputData01
emdSTEM = PRISM01(atoms,cellDim);
emdSTEM = PRISM02_ajp(emdSTEM);
emdSTEM = PRISM03_ajp(emdSTEM);
image = sum(emdSTEM.stack(:,:,14:18),3);
figure, imagesc(image),axis image

