function out = hilbertwavelettransform(data,K,L,scale)
%
% hilbertwavelettransform
%
% Uses, hilbert wavelet basis set estimations by Ivan Selesnick (see hwlet.m)
%
% out= modwt transform using HWP at scale.
%
% data = 1*t vector of unfiltered data.
% K =  - number of zeros at z=-1 (for HWP estimation)
% L = L - degree of fractional delay (for HWP estimation)
% scale = HWP scale for modwt transformation.
%
% PJH 2017
%
data = data - mean(data);
[h0,g0] = hwlet(K,L);

h1=qmf(h0);
g1=qmf(g0);

h = modwt(data,h0,h1, scale);
g = modwt(data,g0,g1, scale);
out = complex(h,g)';
end