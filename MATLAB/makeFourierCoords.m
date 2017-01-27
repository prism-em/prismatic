function [q] = makeFourierCoords(N,pSize)
% This function generates Fourier coordinates 
if mod(N,2) == 0
    q = circshift(((-N/2):(N/2-1))/(N*pSize),[0 -N/2]);
else
    q = circshift(((-N/2+.5):(N/2-.5))/((N-1)*pSize),[0 -N/2+.5]);
end
end