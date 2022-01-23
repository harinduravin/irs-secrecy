clear;
% % Monte Carlo simulation for N scenarios
N = 10000;

% Number of elements of the IRS
L = 50;

% Number of transmitters of the IRS
M = 10;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
dti = 10;
dir = 20;
drt = 25;
pl = 2;

% Rayleigh channels initialized
gti = sqrt(1/2)*(randn(L,M)+1i*randn(L,M));
gir = sqrt(1/2)*(randn(L,1)+1i*randn(L,1));
grt = sqrt(1/2)*(randn(M,1)+1i*randn(M,1));

hti = sqrt(Lo*dti^(-pl))*gti;
hir = sqrt(Lo*dir^(-pl))*gir;
hrt = sqrt(Lo*drt^(-pl))*grt;

% Preparing variables for SDR optimization
phi = diag(hir)*hti;

TL = phi*phi';
TR = phi*hrt;
BL = hrt'*phi';

RL =  cat(1,TL,BL);
RR = cat(1,TR,0);
R = cat(2,RL,RR);

cvx_begin sdp
  variable X(L+1,L+1) symmetric
  maximize(real(trace(R*X)));
  subject to
    diag(X) == 1;
  X == semidefinite(L+1);
cvx_end

[U,D,U_] = svd(X);

% To initiate Gaussian Random vector for Gaussian Randomization method

r = sqrt(1/2)*(randn(L+1,N)+1i*randn(L+1,N));
v_bar = U*sqrt(D)*r;
obj = v_bar'*R*v_bar;

obj_diag = diag(obj);
[value,argu] = max(real(obj_diag));

opt_v_bar = v_bar(:,argu);
opt_v_unnormal = opt_v_bar(1:end-1)/opt_v_bar(end);

% The optimal reflection matrix
opt_v = exp(1i*angle(opt_v_unnormal));

