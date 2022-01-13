clear;
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

% Rayleigh channels initialized
gti = sqrt(1/2)*(randn(L,1)+1i*randn(L,1));
gir = sqrt(1/2)*(randn(L,1)+1i*randn(L,1));
grt = sqrt(1/2)*(randn(1)+1i*randn(1));

hti = sqrt(Lo*dti^(-pl))*gti;
hir = sqrt(Lo*dir^(-pl))*gir;
hrt = sqrt(Lo*drt^(-pl))*grt;

% IRS reflection matrix (Initialized randomly)
wi = (2*rand(L,1)-1)*pi;

% Signal (For Shannon's capacity limit function to work)
s = sqrt(1/2)*(randn(1)+1i*randn(1));

% Transmit power
P = 1;

% Gaussian noise at the receiver
n_ = sqrt(10^(-10.5)*randn(1));
w = [wi ; 1].';
H = cat(1,hti.*hir, hrt);

C = H*H';

cvx_begin sdp
  variable X(L+1,L+1) symmetric
  minimize(trace(C*X));
  subject to
    diag(X) == 1;
  X == semidefinite(L+1);
cvx_end

[U,D,U_] = svd(X);

% To initiate Gaussian Random vector for Gaussian Randomization method

r = sqrt(1/2)*(randn(L+1,N)+1i*randn(L+1,N));
w_bar = U*sqrt(D)*r;


obj = w_bar'*C*w_bar;

obj_diag = diag(obj);


[value,argu] = max(real(obj_diag));

opt_w_bar = w_bar(:,argu);
opt_w_unnormal = opt_w_bar(1:end-1)/opt_w_bar(end);

% The optimal reflection matrix
opt_w = exp(1i*angle(opt_w_unnormal));

