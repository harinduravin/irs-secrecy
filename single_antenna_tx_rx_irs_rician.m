% Monte Carlo simulation for N scenarios
N = 10000;

% Number of elements of the IRS
L = 50;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
dti = 10;
dir = 20;
drt = 25;
pl = 2;

% Rician factor and pathloss between users
beta = 8;
plusers = 3;

s = sqrt(beta)/sqrt(1+beta);

% Rayleigh channels initialized
gir = sqrt(1/2)*(randn(L,N)+1i*randn(L,N));
grt = sqrt(1/2)*(randn(1,N)+1i*randn(1,N));

% Rician channel initialized
gti = sqrt(1/2)*(randn(L,N)+1i*randn(L,N));

hti = sqrt(Lo*dti^(-pl))*gti;
hir = sqrt(Lo*dir^(-pl))*gir;
hrt = sqrt(Lo*drt^(-pl))*grt;

% IRS reflection matrix (Initialized randomly)
wi = (2*rand(L,1)-1)*pi;

% Signal (For Shannon's capacity limit function to work)
s = sqrt(1/2)*(randn(1,N)+1i*randn(1,N));

% Transmit power
P = 1;

% Gaussian noise at the receiver
n = sqrt(10^(-10.5)*randn(1,N));
w = [wi ; 1].';
H = cat(1,hti.*hir, hrt);

sig = sqrt(P)*w*H;
sig_p = abs(sig).^2;
snr = sig_p/10^(-10.5);
% rate = log2(1 + snr);

% A histogram of the rate
histogram(snr,100)
std(snr)

% Mean SNR value for Random reflections ->> 5.1507e+04
% SD SNR value for Random reflections ->> 5.2245e+04


