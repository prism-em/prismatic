load inputData01
emdSTEM = PRISM01(atoms,cellDim);
emdSTEM = PRISM02(emdSTEM);
emdSTEM_col = PRISM03(emdSTEM);
image = sum(emdSTEM_col.stack(:,:,14:18),3);
figure, imagesc(image),axis image