clear;
N = 10000000;

% Normalized to make sure a valid comparison is possible
var = 1/32;

mu1 = sqrt(15/32);
mu2 = sqrt(15/32);

x1 = sqrt(var)*randn([1 N]) + mu1;
x2 = sqrt(var)*randn([1 N]) + mu2;

y = x1 + 1i*(x2);
r = abs(y);

x_step = 1/5000;
x = 0:x_step:5;
n = hist(r,x);
r1 = (n./(N*x_step));

% Theoretical rician fading
s = sqrt(mu1^2+mu2^2);
r2 = (x/(var)).*exp(-(x.^2+s^2)/(2*var)).*besseli(0,x*s/var);

% Simulation plot of rayleigh
plot(x,r1)
grid
sum(r1*x_step) % The area inside the graph adds to 1
hold
% Theoretical plot of rayleigh
plot(x,r2,'r')

