function [newvec] = LoPaFi(oldvec,convolution)
% low-pass filter with Butterworth filtering (conv=0)
% dim=4,cutoff=1 (Reimer et al., 2014, Neuron)
% dim=3,cutoff=4 (de Gee et al., 2014, PNAS)

switch convolution
    case 0
%       dim=3;cutoff=4;
      dim=4;cutoff=1;
      [b,a] = butter(dim,2*cutoff/500);

      newvec = filter(b,a,oldvec);
    case 1
        k = 5;
        kar = ones(1,k)./k;
        newvec = conv(oldvec,kar,'same');
end