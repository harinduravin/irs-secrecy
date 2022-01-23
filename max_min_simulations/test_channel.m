
lambda = 0.06;
pl = 2;
Lo = 10^(-3);
L = 10;
seed = 0;

HAA = get_H(18.46,15.37,31.81,lambda,L,Lo,pl,10*seed);
HBB = get_H(17.95,15.23,25.98,lambda,L,Lo,pl,10*seed+1);
H11 = get_H(18.46,17.95,5.57,lambda,L,Lo,pl,10*seed+2);
H21 = get_H(15.37,17.95,29.24,lambda,L,Lo,pl,10*seed+3);
H12 = get_H(18.46,15.23,29.26,lambda,L,Lo,pl,10*seed+4);
H22 = get_H(15.37,15.23,5.25,lambda,L,Lo,pl,10*seed+5);
HA1C = get_H(18.46,20.87,20.87,lambda,L,Lo,pl,10*seed+6);
HA2C = get_H(15.37,20.87,22.80,lambda,L,Lo,pl,10*seed+7);
HB1C = get_H(17.95,20.87,15.47,lambda,L,Lo,pl,10*seed+8);
HB2C = get_H(15.23,20.87,17.81,lambda,L,Lo,pl,10*seed+9);

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

    fprintf(string(rng_val_))
    fprintf('\n')
    fprintf(string(rng_val_+1))
    fprintf('\n')
    rng(rng_val_);
    real_part = randn(L,1);
    rng(rng_val_+1);
    img_part = randn(L,1);
    rayleigh_channel = sqrt(1/(2*dist^pl))*(real_part+1i*img_part);
%     rayleigh_channel = sqrt(1/(2*dist^pl))*(randn(L,1)+1i*randn(L,1));
end

function rayleigh_direct_channel = get_rayleigh_direct_channel(dist, pl, rng_val_)

    fprintf(string(rng_val_))
    fprintf('\n')
    rng(rng_val_)
    rayleigh_direct_channel = sqrt(1/(2*dist^pl))*(randn(1)+1i*randn(1));
end
