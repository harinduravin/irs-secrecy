clear;
% Wavelength

lambda = 0.06;

% Number of elements of the IRS
L = 10;

% Random arrays for Gaussian Randomization
N = 10000;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
dti = 5.17;
dir = 20;
drt = 25;
pl = 2.5;

% Rayleigh channels initialized

HAA = get_H(9.94,19.63,11.66,lambda,L,Lo,pl);
HBB = get_H(8.39,18.06,12.20,lambda,L,Lo,pl);
H11 = get_H(9.94,8.39,14.89,lambda,L,Lo,pl);
 H21 = get_H(19.63,8.39,20.17,lambda,L,Lo,pl);
H12 = get_H(9.94,18.06,17.91,lambda,L,Lo,pl);
H22 = get_H(19.63,18.06,14.89,lambda,L,Lo,pl);
HA1C = get_H(9.94,55.30,47.60,lambda,L,Lo,pl);
HA2C = get_H(19.63,55.30,36.04,lambda,L,Lo,pl);
HB1C = get_H(8.39,55.30,53.37,lambda,L,Lo,pl);
HB2C = get_H(18.06,55.30,42.36,lambda,L,Lo,pl);

% IRS reflection matrix (Initialized randomly)
wi = exp(1i*(2*rand(L,1)-1)*pi);
w = [wi ; 1].';
W = w'*w;


% Noise and residual loop interference
sigma_ab = 10^(-14.5);
sigma_loop = 10^(-14);
sigma_c = 10^(-14.5);

for j = 0:10

    max_min = 0;
    max_min_ind = 0;
    for i = 0:15
        power_string = num2str(dec2bin(i,4));
        PA1 = get_power(power_string(1));
        PA2 = get_power(power_string(2));
        PB1 = get_power(power_string(3));
        PB2 = get_power(power_string(4));

        RA1 = get_P(W, H11, HAA, H12, PB1, PA2, PB2) - get_Q(W, HAA, H12,HB1C, PA2, PB2 ,PB1);
        RA2 = get_P(W, H22, HAA, H21, PB2, PA1, PB1) - get_Q(W, HAA, H21,HB2C, PA1, PB1,PB2);
        RB1 = get_P(W, H11, HBB, H21, PA1, PB2, PA2) - get_Q(W, HBB, H21,HA1C, PB2, PA2,PA1);
        RB2 = get_P(W, H22, HBB, H12, PA2, PB1, PA1) - get_Q(W, HBB, H12,HA2C, PB1, PA1,PA2);

        min_value = min([RA1,RA2,RB1,RB2]);
%         fprintf(string(min_value))
%         fprintf(" ")
        if (min_value > max_min)
            
            max_min = min_value;
            max_min_ind = i;
        end

    end
%     fprintf("done")

    power_string = num2str(dec2bin(max_min_ind,4));
    PA1 = get_power(power_string(1));
    PA2 = get_power(power_string(2));
    PB1 = get_power(power_string(3));
    PB2 = get_power(power_string(4));

    cvx_begin quiet
    cvx_solver mosek
      variable X(L+1,L+1) complex semidefinite %symmetric
      variable t

      minimize t
    
      subject to
        -log(PB1*real(trace((H11*H11')*X))+PA2*real(trace((HAA*HAA')*X))+PB2*real(trace((H12*H12')*X))+sigma_ab+sigma_loop) - real(trace(get_grad_S(W,  PA2, HAA, PB2, H12, PB1, HB1C)*(X-W))) <= t;
        -log(PB2*real(trace((H22*H22')*X))+PA1*real(trace((HAA*HAA')*X))+PB1*real(trace((H21*H21')*X))+sigma_ab+sigma_loop) - real(trace(get_grad_S(W,  PA1, HAA, PB1, H21, PB2, HB2C)*(X-W))) <= t;
        -log(PA1*real(trace((H11*H11')*X))+PB2*real(trace((HBB*HBB')*X))+PA2*real(trace((H21*H21')*X))+sigma_ab+sigma_loop) - real(trace(get_grad_S(W,  PB2, HBB, PA2, H21, PA1, HA1C)*(X-W))) <= t;
        -log(PA2*real(trace((H22*H22')*X))+PB1*real(trace((HBB*HBB')*X))+PA1*real(trace((H12*H12')*X))+sigma_ab+sigma_loop) - real(trace(get_grad_S(W,  PB1, HBB, PA1, H12, PA2, HA2C)*(X-W))) <= t;
        diag(X) == 1;
        norm((X-W),1)<=2;
    cvx_end
    fprintf(string(t))

    W = X;


end

[U,D,U_] = svd(X);

% % To initiate Gaussian Random vector for Gaussian Randomization method
% 
% r = sqrt(1/2)*(randn(L+1,N)+1i*randn(L+1,N));
% w_bar = U*sqrt(D)*r;
% obj = w_bar'*(H11*H11')*w_bar;
% 
% obj_diag = diag(obj);
% [value,argu] = max(real(obj_diag));
% 
% opt_w_bar = w_bar(:,argu);
% opt_w_unnormal = opt_w_bar(1:end-1)/opt_w_bar(end);
% 
% % The optimal reflection matrix
% opt_w = exp(1i*angle(opt_w_unnormal));

% Taking the 1st column instead of using Gaussian Randomization

w0 = U(:,1)/U(L+1,1);
w0 = w0./abs(w0);
Wx = w0*w0';

% s = ((log(1 + (opt_P2/(sig1^2+sigl1^2))*real(trace(H*Wx)))+log(1 + (opt_P1/(sig2^2+sigl2^2))*real(trace(G*Wx)))))-log(1+(opt_P1/sige^2)*real(trace(E1*Wx))+(opt_P2/sige^2)*real(trace(E2*Wx)));
% s = s/(log(2));
RA1 = get_secrecy_rate(Wx, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1)
RA2 = get_secrecy_rate(Wx, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2)
RB1 = get_secrecy_rate(Wx, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1)
RB2 = get_secrecy_rate(Wx, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2)

function H = get_H(dti,dir,dtr,wave_l,L,Lo,pl)

    gti = get_rayleigh_channel(dti, pl,L);
    hti = sqrt(Lo)*gti;
    gir = get_rayleigh_channel(dir, pl,L);
    hir = sqrt(Lo)*gir;
    gtr = get_rayleigh_direct_channel(dtr, pl);
    htr = sqrt(Lo)*gtr;
    H = cat(1,hti.*hir, htr);
end

function grad_S = get_grad_S(W, pinf1, inf1, pinf2, inf2, psec, sec)
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    sigma_c = 10^(-14.5);

    term1 = (pinf1*(inf1*inf1') + pinf2*(inf2*inf2'))/(get_A(W, pinf1, inf1, pinf2, inf2)+sigma_ab+sigma_loop);
    term2 = (psec*(sec*sec'))/(get_B(W, psec, sec)+sigma_c);

    grad_S = (-1/log(2))*(term1 + term2)';
end

function channel_term = get_channel_term(power, W, H)
    H_tilda = H*H';

    channel_term = power*(real(trace(H_tilda*W)));
end

function A = get_A(W, pinf1, inf1, pinf2, inf2)

    A = get_channel_term(pinf1, W, inf1)+get_channel_term(pinf2, W, inf2);
end

function B = get_B(W, psec, sec)

    B = get_channel_term(psec, W, sec)+get_channel_term(psec, W, sec);
end

% function rice_channel = get_rice_channel(dist, wave_l, L)
%     rice_channel = sqrt(1/32)*(randn(L,1)+1i*randn(1)) + sqrt(30/32)*exp(1i*2*pi*(dist/wave_l-floor(dist/wave_l)));
% end
% 
% function rice_direct_channel = get_rice_direct_channel(dist, wave_l)
%     rice_direct_channel = sqrt(1/32)*(randn(1)+1i*randn(1)) + sqrt(30/32)*exp(1i*2*pi*(dist/wave_l-floor(dist/wave_l)));
% end

% H =sqrt(1/2*d^loss_exponent)*(randn(1,1)+1i*randn(1,1)) SISO

function rayleigh_channel = get_rayleigh_channel(dist, pl, L)
    rayleigh_channel = sqrt(1/(2*dist^pl))*(randn(L,1)+1i*randn(1));
end

function rayleigh_direct_channel = get_rayleigh_direct_channel(dist, pl)
    rayleigh_direct_channel = sqrt(1/(2*dist^pl))*(randn(1)+1i*randn(1));
end

function c_rate = get_secrecy_rate(W, leg, inf1, inf2, eves, pleg, pinf1, pinf2, peves)

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    eves_inf = log2(1 + (get_channel_term(peves, W, eves))/(sigma_c));
    leg_inf = log2(1 + (get_channel_term(pleg, W, leg))/(get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop));
    c_rate = leg_inf - eves_inf;
%     c_rate = max(0,leg_inf - eves_inf);

end

function P = get_P(W, leg, inf1, inf2, pleg, pinf1, pinf2)

    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    P = log(get_channel_term(pleg, W, leg)+get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop);

end

function Q = get_Q(W,inf1, inf2, eves,pinf1, pinf2,peves)
    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    Q = log(get_channel_term(pinf1, W, inf1)+get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop)+log(get_channel_term(peves, W, eves) + sigma_c);

end

function power = get_power(string_value)
    pmin = 0.001;
    pmax = 0.04;
    if string_value == '0' 
        power = pmin;
    else 
        power = pmax;
    end
end