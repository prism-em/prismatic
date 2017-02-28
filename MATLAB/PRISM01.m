function [emdSTEM] = PRISM01(atoms,cellDim)

% 01 - calculate the projected potential from a collection of an orthogonal
%      unit cell, sliced up for multislice calculations.

% Inputs
% atoms -   N x 4 array containing fractional atomic coordinates, with rows
%           given by [x y z atomic_number], for N total atoms.
            
% cellDim - Dimensions of the unit cell [x_cell y_cell z_cell]

pixelSize = 100/1000;   % Realspace pixel size.
emdSTEM.potBound = 1;       % Radial distance to integrate atomic potentials.
emdSTEM.numFP = 8/8;          % Number of frozen phonon configurations.
emdSTEM.sliceThickness = 2; % Thickness of each potential slice.
emdSTEM.interpolationFactor = 100;

% u = ones(118,1) * 0.08;      % Debye waller coefficients.
u = ones(118,1) * 0.0;      % Debye waller coefficients.

% Keep atomic positions in struct
emdSTEM.atoms = atoms;
emdSTEM.cellDim = cellDim;

% Simulation size
f = 4*emdSTEM.interpolationFactor;
emdSTEM.imageSize = round(cellDim(1:2)/pixelSize/f)*f;
emdSTEM.pixelSize = cellDim(1:2) ./ emdSTEM.imageSize;


% Construct projected potentials
xyLeng = ceil(emdSTEM.potBound./emdSTEM.pixelSize);
xvec = -xyLeng(1):xyLeng(1);
yvec = -xyLeng(2):xyLeng(2);
xr = xvec*emdSTEM.pixelSize(1);
yr = yvec*emdSTEM.pixelSize(2);
% Lookup table for atom types
atomTypes = unique(atoms(:,4));
potLookup = zeros(length(xvec),length(yvec),length(atomTypes));
uLookup = zeros(length(atomTypes),1);
for a0 = 1:length(atomTypes)
    potLookup(:,:,a0) = projPot(atomTypes(a0),xr,yr);
    if u(atomTypes(a0)) > 0
        uLookup(a0) = u(atomTypes(a0));
    else
        disp(['Warning, no RMS displacement given for atomic number ' ...
            num2str(atomTypes(a0)) '!  Setting u = 0.05.'])
        uLookup(a0) = 0.05;
    end
end

% Convert atomic coordinates into Cartesian
x = atoms(:,1)*cellDim(1);
y = atoms(:,2)*cellDim(2);
z = atoms(:,3)*cellDim(3);
ID = atoms(:,4);

% Determine z plane slice index for all atoms
% zPlane = round((z - min(z))/emdSTEM.sliceThickness + 0.5);
zPlane = round((-z + max(z))/emdSTEM.sliceThickness + 0.5);
emdSTEM.numPlanes = max(zPlane);

% Generate projected potentials for all atoms /  frozen phonon configs
emdSTEM.pot = zeros(emdSTEM.imageSize(1),...
    emdSTEM.imageSize(2),emdSTEM.numPlanes,emdSTEM.numFP);
        potProj = zeros(emdSTEM.imageSize);
for a0 = 1:emdSTEM.numPlanes
    inds = find(zPlane == a0);
    
    for a1 = 1:emdSTEM.numFP
        potProj(:) = 0;
        
         for a2 = 1:length(inds)
            [~,indType] = min(abs(atomTypes-ID(inds(a2))));
            xp = mod(xvec+round((x(inds(a2))...
                + 0*randn*uLookup(indType)) ...
                /emdSTEM.pixelSize(1)),emdSTEM.imageSize(1))+1;
            yp = mod(yvec+round((y(inds(a2)) ...
                + 0*randn*uLookup(indType)) ...
                /emdSTEM.pixelSize(2)),emdSTEM.imageSize(2))+1;

%             xp = mod(xvec+round((x(inds(a2))...
%                 + randn*uLookup(indType)) ...
%                 /emdSTEM.pixelSize(1)),emdSTEM.imageSize(1))+1;
%             yp = mod(yvec+round((y(inds(a2)) ...
%                 + randn*uLookup(indType)) ...
%                 /emdSTEM.pixelSize(2)),emdSTEM.imageSize(2))+1;
            potProj(xp,yp) = potProj(xp,yp) + potLookup(:,:,indType);
        end
        
        emdSTEM.pot(:,:,a0,a1) = potProj;
    end
    
    comp = a0 / emdSTEM.numPlanes;
    progressbar(comp,2);
end




end