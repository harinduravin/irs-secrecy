
% Number of IRS elements
L = 40;

% Number of channel realizations
N = 200;

% Residual interference and noise values
sigma_c = dbm2watt(-105);
sigma_ab = dbm2watt(-105);
sigma_loop = dbm2watt(-100);

% Config 5
% P_a = [-30 0 0];
% P_b = [40 0 0];
% P_a2 = [-30 -15 0];
% P_b2 = [40 -15 0];
P_a = [-20 -10 0];
P_a2 = [20 -10 0];
P_b = [-20 -40 0];
P_b2 = [20 -40 0];
P_r = [0 10 0];
P_c1 = [0 5 0];
% P_c2 = [32 2 0];

P_c = P_c1;


d_cr = norm(P_r - P_c);
d_ab = norm(P_a - P_b);
d_ar = norm(P_a - P_r);
d_br = norm(P_b - P_r);
d_ac = norm(P_a - P_c);
d_bc = norm(P_b - P_c);

d_a2r = norm(P_a2 - P_r);
d_b2r = norm(P_b2 - P_r);
d_a2b2 = norm(P_a2 - P_b2);
d_a2c = norm(P_a2 - P_c);
d_b2c = norm(P_b2 - P_c);

d_a2a = norm(P_a2 - P_a);
d_a2b = norm(P_a2 - P_b);
d_b2a = norm(P_b2 - P_a);
d_b2b = norm(P_b2 - P_b);

Hab = zeros(L+1,N);
Hac = zeros(L+1,N);
Hbc = zeros(L+1,N);
Ha2b2 = zeros(L+1,N);
Ha2c = zeros(L+1,N);
Hb2c = zeros(L+1,N);
Ha2a = zeros(L+1,N);
Ha2b = zeros(L+1,N); 
Hb2a = zeros(L+1,N);
Hb2b = zeros(L+1,N); 

for i = 1:N
    Hab(:,i) = get_H(d_ar,d_br,d_ab,L,10*i);
    Hac(:,i) = get_H(d_ar,d_cr,d_ac,L,10*i-1);
    Hbc(:,i) = get_H(d_br,d_cr,d_bc,L,10*i-2);
    Ha2b2(:,i) = get_H(d_a2r,d_b2r,d_a2b2,L,10*i-3);
    Ha2c(:,i) = get_H(d_a2r,d_cr,d_a2c,L,10*i-4);
    Hb2c(:,i) = get_H(d_b2r,d_cr,d_b2c,L,10*i-5);
    Ha2a(:,i) = get_H(d_a2r,d_ar,d_a2a,L,10*i-6);
    Ha2b(:,i) = get_H(d_a2r,d_br,d_a2b,L,10*i-7); 
    Hb2a(:,i) = get_H(d_b2r,d_ar,d_b2a,L,10*i-8);
    Hb2b(:,i) = get_H(d_b2r,d_br,d_b2b,L,10*i-9); 
end

save(strcat('L',num2str(L),'channels_config2'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = get_H(dti,dir,dtr,L,rng_val)
    % Generates H. This matrix stores both direct channel and the channel through
    % the IRS. Two channels are concatenated to form a single column matrix.
    %
    % Using 3GPP Urban Micro pathloss models; A frequency of 750 MHz is assumed.


    % LOS channel : Used for Users to IRS channels
    % Lo2 = 10^(2*-2.55012);
    Lo2 = 0.001;
    % NLOS channel : Used for direct path
    % Lo3 = 10^(2*-1.94515);
    Lo3 = 0.001;

    gti = get_ricean_channel(dti, 2, L, 5*rng_val,0);
    gir = get_ricean_channel(dir, 2, L, 5*rng_val+2,0); %2.2
    gtr = get_ricean_direct_channel(dtr, 3, 5*rng_val+4,8); %3.67
    htr = sqrt(Lo3)*gtr;
    H = cat(1,Lo2*gti.*gir, htr);
end

function ricean_channel = get_ricean_channel(dist, pl, L,rng_val_,kappa)
    % Generates rayleigh fading channel through the IRS
    %

    g1 = sqrt(kappa/(1+kappa));
    g2 = sqrt(1/(2*(1+kappa)));

    rng(rng_val_);
    real_part = randn(L,1);
    rng(rng_val_+1);
    img_part = randn(L,1);
    cg = real_part+1i*img_part;
    ricean_channel = sqrt(1/(dist^pl))*(g1+g2*cg);

end

function ricean_direct_channel = get_ricean_direct_channel(dist, pl, rng_val_,kappa)
    % Generates ricean fading channel for the direct path
    %

    g1 = sqrt(kappa/(1+kappa));
    g2 = sqrt(1/(2*(1+kappa)));

    rng(rng_val_)
    cg = (randn(1)+1i*randn(1));
    ricean_direct_channel = sqrt(1/(dist^pl))*(g1+g2*cg);
end

function  removed = remove_small(X)
    indices = abs(X)<=(10^(-7));
    X(indices) = 0 + 0*1i;
    removed = X;
end

function power = dbm2watt(dbm_value)
    % Converting power values from dBm to Watt for calculations

    power = (10^(dbm_value/10))*(10^(-3));
end