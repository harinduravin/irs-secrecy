clear;
N = 1000000;
p = 1/sqrt(2);
x1 = p*randn([1 N]);
x2 = p*randn([1 N]);

% A vector of Rayleigh distributed random variables (r)
y = x1 + 1i*(x2);
r = abs(y);

x_step = 1/500;
x = 0:x_step:5;
n = hist(r,x);

r1 = (n./(N*x_step));

% Theoretical rayleigh fading

r2 = (x/(p^2)).*exp(-(x.^2)/(2*p^2));

% Simulation plot of rayleigh
plot(x,r1)
grid
sum(r1*x_step) % The area inside the graph adds to 1
hold
% Theoretical plot of rayleigh
plot(x,r2,'r')

