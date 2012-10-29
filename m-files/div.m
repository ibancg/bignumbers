
an = 583849;
bn = 374817;
cn = an*bn;

a0 = [ 5 8 3 8 4 9];
b0 = [ 3 7 4 8 1 7];
n = length(a0);
a = [ zeros(1,n) a0 ]; 
b = [ zeros(1,n) b0 ];
c = [ 2 1 8 8 3 6 5 3 0 6 3 3 ] 
F = fft(eye(2*n));

Fa = F*a.';
Fb = F*b.';
Fc = Fa.*Fb
%Fc = F*c.'

c = real(F \ Fc).';

Far = [ real(Fa(1:n+1)); imag(Fa(2:n)) ];
Fbr = [ real(Fb(1:n+1)); imag(Fb(2:n)) ];
Fcr = [ real(Fc(1:n+1)); imag(Fc(2:n)) ];

R = [eye(n+1) zeros(n+1, n-1); zeros(n-1,1) rot90(eye(n-1)) zeros(n-1, n)];
I = full([zeros(1, 2*n); zeros(n-1, n+1) eye(n-1); zeros(1,2*n); zeros(n-1, n+1) -rot90(eye(n-1))]);
C = R + 1i*I;

norm(C*Far - Fa)
norm(C*Fbr - Fb)
norm(C*Fcr - Fc)

Fabr_r = (R*Far).*(R*Fbr) - (I*Far).*(I*Fbr);
Fabr_i = (R*Far).*(I*Fbr) + (I*Far).*(R*Fbr);

R2 = [eye(n+1) zeros(n+1,n-1); zeros(n-1,2*n)];
I2 = [zeros(n+1, 2*n); zeros(n-1, 1) eye(n-1) zeros(n-1, n)];

Fcr2 = R2*Fabr_r + I2*Fabr_i;

norm(Fcr2 - Fcr)

A = diag(R*Far);
B = diag(I*Far);

K = R2*(A*R - B*I) + I2*(A*I + B*R);
F\(C*(K\Fcr))




