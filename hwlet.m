function [h,g] = hwlet(K,L)
% Hilbert transform pair of orthogonal wavelet bases
% h, g - scaling filters of length 2*(K+L)
% K - number of zeros at z=-1
% L - degree of fractional delay
%
%
%

n = 0:L-1;
t = 1/2;
d = cumprod([1, (L-n).*(L-n-t)./(n+1)./(n+1+t)]);
s1 = binom(2*K,0:2*K);
s2 = conv(d,d(end:-1:1));
s = conv(s1,s2);
M = K+L;
C = convmtx(s',2*M-1);
C = C(2:2:end,:);
b = zeros(2*M-1,1);
b(M) = 1;
r = (C\b)';
q = sfact(r);
f = conv(q,binom(K,0:K));
h = conv(f,d);
g = conv(f,d(end:-1:1));
end
function a = binom(n,k)
%
% a = binom(n,k)
% BINOMIAL COEFFICIENTS
%
% allowable inputs:
%       n : integer, k : integer
%       n : integer vector, k : integer
%       n : integer, k : integer vector
%       n : integer vector, k : integer vector (of equal dimension)
%
%
% Code Reproduced from: 
%
% Ivan Selesnick
% selesi@taco.poly.edu
% Polytechnic University

nv = n;
kv = k;
if (length(nv) == 1) & (length(kv) > 1)
        nv = nv * ones(size(kv));
elseif (length(nv) > 1) & (length(kv) == 1)
        kv = kv * ones(size(nv));
end
a = nv;
for i = 1:length(nv)
   n = nv(i);
   k = kv(i);
   if n >= 0
      if k >= 0
         if n >= k
            c = prod(1:n)/(prod(1:k)*prod(1:n-k));
         else
            c = 0;
        end
     else
        c = 0;
     end
   else
      if k >= 0
         c = (-1)^k * prod(1:k-n-1)/(prod(1:k)*prod(1:-n-1));
      else
         if n >= k
            c = (-1)^(n-k)*prod(1:-k-1)/(prod(1:n-k)*prod(1:-n-1));
         else
            c = 0;
         end
      end
   end
   a(i) = c;
end
end
function [b,r] = sfact(h)
% [b,r] = sfact(h)
% spectral factorization of a polynomial h.
% b: new polynomial
% r: roots of new polynomial
%
% % example:
%    g = rand(1,10);
%    h = conv(g,g(10:-1:1));
%    b = sfact(h);
%    h - conv(b,b(10:-1:1)) % should be 0

% required subprograms: seprts.m, leja.m
%
% Ivan Selesnick, Polytechnic University, Brooklyn, NY
% selesi@taco.poly.edu
%
% leja.m is by Markus Lang, and is available from the 
% Rice University DSP webpage: www.dsp.rice.edu

if length(h) == 1
	b = sqrt(h);
	r = [];
	return
end

% Get the appropriate roots.
r = seprts(h);	

% Form the polynomial from the roots
r = leja(r);
b = poly(r);	

if isreal(h)
	b = real(b);
end

% normalize
b = b*sqrt(max(h)/sum(abs(b).^2));

if max(b)+min(b) < 0
   b = -b;
end

end
% ------------------------------------------------------------


function r = seprts(p)
% r = seprts(p)
% This program is for spectral factorization.
% The roots on the unit circle must have even degree.
% Roots with high multiplicity will cause problems,
% they should be handled by extracting them prior to
% using this program.

% Ivan Selesnick, Polytechnic University, Brooklyn, NY
% selesi@taco.poly.edu

SN = 0.0001;    % Small Number (criterion for deciding if a 
                % root is on the unit circle).

rts = roots(p);  

% The roots INSIDE the unit circle
irts = rts(abs(rts)<(1-SN));    

% The roots ON the unit circle            
orts = rts((abs(rts)>=(1-SN)) & (abs(rts)<=(1+SN)));
N = length(orts);
if rem(N,2) == 1
        disp('Sorry, but there is a problem (1) in seprts.m')
        r = [];
	return
end

% Sort roots on the unit circle by angle
[a,k] = sort(angle(orts));
orts = orts(k(1:2:end));

% Make final list of roots
r = [irts; orts];
end

% ------------------------------------------------------------


function [x_out] = leja(x_in)
% function [x_out] = leja(x_in)
%
%    Input:   x_in
%
%    Output:  x_out
%
%    Program orders the values x_in (supposed to be the roots of a 
%    polynomial) in this way that computing the polynomial coefficients
%    by using the m-file poly yields numerically accurate results.
%    Try, e.g., 
%               z=exp(j*(1:100)*2*pi/100);
%    and compute  
%               p1 = poly(z);
%               p2 = poly(leja(z));
%    which both should lead to the polynomial x^100-1. You will be
%    surprised!
%
%

%File Name: leja.m
%Last Modification Date: %G%	%U%
%Current Version: %M%	%I%
%File Creation Date: Mon Nov  8 09:53:56 1993
%Author: Markus Lang  <lang@dsp.rice.edu>
%
%Copyright: All software, documentation, and related files in this distribution
%           are Copyright (c) 1993  Rice University
%
%Permission is granted for use and non-profit distribution providing that this
%notice be clearly maintained. The right to distribute any portion for profit
%or as part of any commercial product is specifically reserved for the author.
%
%Change History:
%

x = x_in(:).'; n = length(x);

a = x(ones(1,n+1),:);
a(1,:) = abs(a(1,:));

[dum1,ind] = max(a(1,1:n));  
if ind~=1
  dum2 = a(:,1);  a(:,1) = a(:,ind);  a(:,ind) = dum2;
end
x_out(1) = a(n,1);
a(2,2:n) = abs(a(2,2:n)-x_out(1));

for l=2:n-1
  [dum1,ind] = max(prod(a(1:l,l:n)));  ind = ind+l-1;
  if l~=ind
    dum2 = a(:,l);  a(:,l) = a(:,ind);  a(:,ind) = dum2;
  end
  x_out(l) = a(n,l);
  a(l+1,(l+1):n) = abs(a(l+1,(l+1):n)-x_out(l));
end
x_out = a(n+1,:);
end