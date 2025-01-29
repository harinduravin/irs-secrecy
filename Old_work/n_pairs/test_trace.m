
lambda = 0.06;
Lo = 10^(-3);
pl = 2;
L = 20;
seed = 3;

% IRS reflection matrix (Initialized randomly)
wi = exp(1i*(2*rand(L,1)-1)*pi);
w = [wi ; 1].';
W = w'*w;

HAA = get_H(34.33,8.43,33.43,lambda,L,Lo,pl,10*seed);
HBB = get_H(35.66,10.44,43.46,lambda,L,Lo,pl,10*seed+1);
H11 = get_H(34.33,35.66,9.13,lambda,L,Lo,pl,10*seed+2);
H21 = get_H(8.43,35.66,36.88,lambda,L,Lo,pl,10*seed+3);

value1 =  get_channel_term(0.04, W, HAA)+ get_channel_term(0.001, W, HBB)+ get_channel_term(0.04, W, H11)+ get_channel_term(0.001, W, H21)
value2 =  get_channel_terms_four(0.04, 0.001, 0.04, 0.001, W, HAA, HBB, H11, H21)

function channel_term = get_channel_term(power, W, H)
    H_tilda = H*H';

    channel_term = power*(real(trace(H_tilda*W)));
end

function channel_term = get_channel_term_double(power1, power2, W, H1, H2)
    H_tilda1 = H1*H1';
    H_tilda2 = H2*H2';

    channel_term = trace((power1*H_tilda1+power2*H_tilda2)*W);
end

function channel_term = get_channel_terms_four(power1, power2, power3, power4, W, H1, H2, H3, H4)

    H_tilda1 = H1*H1';
    H_tilda2 = H2*H2';
    H_tilda3 = H3*H3';
    H_tilda4 = H4*H4';

    channel_term = trace((power1*H_tilda1+power2*H_tilda2+power3*H_tilda3+power4*H_tilda4)*W);
end

function H = get_H(dti,dir,dtr,wave_l,L,Lo,pl,rng_val)

    gti = get_rayleigh_channel(dti, pl, L, 5*rng_val);
    hti = sqrt(Lo)*gti;
    gir = get_rayleigh_channel(dir, pl, L, 5*rng_val+2);
    hir = sqrt(Lo)*gir;
    gtr = get_rayleigh_direct_channel(dtr, pl, 5*rng_val+4);
    htr = sqrt(Lo)*gtr;
    H = cat(1,hti.*hir, htr);
end

function rayleigh_channel = get_rayleigh_channel(dist, pl, L,rng_val_)

    rng(rng_val_);
    real_part = randn(L,1);
    rng(rng_val_+1);
    img_part = randn(L,1);
    rayleigh_channel = sqrt(1/(2*dist^pl))*(real_part+1i*img_part);
%     rayleigh_channel = sqrt(1/(2*dist^pl))*(randn(L,1)+1i*randn(L,1));
end

function rayleigh_direct_channel = get_rayleigh_direct_channel(dist, pl, rng_val_)
    rng(rng_val_)
    rayleigh_direct_channel = sqrt(1/(2*dist^pl))*(randn(1)+1i*randn(1));
end