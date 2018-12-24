function cp = rocN(x,y)
N = 100;
[m n] = size(x);
x = reshape(x,1,m*n);
[m n] = size(y);
y = reshape(y,1,m*n);

zlo = min([min(x(:)) min(y(:))]);
zhi = max([max(x(:)) max(y(:))]);
z = linspace(zlo,zhi,N);
fa = zeros(1,N);	% allocate the vector
hit = zeros(1,N);
for i = 1:N
  fa(N-i+1) = sum(y > z(i));
  hit(N-i+1) = sum(x > z(i));
end
[m,ny] = size(y);
fa = fa/ny;
[m,nx] = size(x);
hit = hit/nx;
fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;
cp = trapz(fa,hit);