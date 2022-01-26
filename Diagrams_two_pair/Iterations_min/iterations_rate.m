clear;
% Wavelength

lambda = 0.06;

% Random arrays for Gaussian Randomization
N = 10000;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
pl = 2;

L = 20;
inner_iter = 8;
outer_iter = 3;

iter_list = 1:inner_iter*outer_iter;

A1_master_list = [];
B1_master_list = [];
A2_master_list = [];
B2_master_list = [];
min_master_list = [];

num_seeds = 100;
f = waitbar(0,'Please wait...');

for seed = 1:num_seeds
    waitbar(seed/num_seeds,f,'Simulating...');

    A1_list = [];
    B1_list = [];
    A2_list = [];
    B2_list = [];
    min_list = [];


    % Rayleigh channels initialized

    HAA = get_H(34.33,8.43,33.43,lambda,L,Lo,pl,10*seed);
    HBB = get_H(35.66,10.44,43.46,lambda,L,Lo,pl,10*seed+1);
    H11 = get_H(34.33,35.66,9.13,lambda,L,Lo,pl,10*seed+2);
    H21 = get_H(8.43,35.66,36.88,lambda,L,Lo,pl,10*seed+3);
    H12 = get_H(34.33,10.44,40.32,lambda,L,Lo,pl,10*seed+4);
    H22 = get_H(8.43,10.44,6.92,lambda,L,Lo,pl,10*seed+5);
    HA1C = get_H(34.33,26.47,9.36,lambda,L,Lo,pl,10*seed+6);
    HA2C = get_H(8.43,26.47,27.29,lambda,L,Lo,pl,10*seed+7);
    HB1C = get_H(35.66,26.47,9.60,lambda,L,Lo,pl,10*seed+8);
    HB2C = get_H(10.44,26.47,33.90,lambda,L,Lo,pl,10*seed+9);

    % IRS reflection matrix (Initialized randomly)
    wi = exp(1i*(2*rand(L,1)-1)*pi);
    w = [wi ; 1].';
    W = w'*w;

    PA1 = 0.04;
    PB1 = 0.04;
    PA2 = 0.04;
    PB2 = 0.04;

    RA1 = get_secrecy_rate(W, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C);
    RA2 = get_secrecy_rate(W, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C);
    RB1 = get_secrecy_rate(W, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C);
    RB2 = get_secrecy_rate(W, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C);

    % Noise and residual loop interference
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    sigma_c = 10^(-14.5);

    min_rate = 0;

    for j = 1:outer_iter

        max_min = 0;
        max_min_ind = 0;
        for i = 0:15
            power_string = num2str(dec2bin(i,4));
            PA1 = get_power(power_string(1));
            PA2 = get_power(power_string(2));
            PB1 = get_power(power_string(3));
            PB2 = get_power(power_string(4));

            RA1 = get_P(W, H11, HAA, H12, PB1, PA2, PB2, PA1, PB2, PA2, HA1C, HB2C, HA2C) - get_Q(W, HAA, H12,HB1C, PA2, PB2 ,PB1, PA1, PB2, PA2,HA1C, HB2C, HA2C);
            RA2 = get_P(W, H22, HAA, H21, PB2, PA1, PB1, PA1, PB1, PA2, HA1C, HB1C, HA2C) - get_Q(W, HAA, H21,HB2C, PA1, PB1 ,PB2, PA1, PB1, PA2,HA1C, HB1C, HA2C);
            RB1 = get_P(W, H11, HBB, H21, PA1, PB2, PA2, PB1, PB2, PA2, HB1C, HB2C, HA2C) - get_Q(W, HBB, H21,HA1C, PB2, PA2 ,PA1, PB1, PB2, PA2,HB1C, HB2C, HA2C);
            RB2 = get_P(W, H22, HBB, H12, PA2, PB1, PA1, PA1, PB1, PB2, HA1C, HB1C, HB2C) - get_Q(W, HBB, H12,HA2C, PB1, PA1 ,PA2, PA1, PB1, PB2,HA1C, HB1C, HB2C);

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

        for j = 1:inner_iter

            cvx_begin quiet
            cvx_solver mosek
            variable X(L+1,L+1) complex semidefinite %symmetric
            variable t

            minimize t

            subject to
                -log(PB1*real(trace((H11*H11')*X))+PA2*real(trace((HAA*HAA')*X))+PB2*real(trace((H12*H12')*X))+sigma_ab+sigma_loop) -log(PA1*real(trace((HA1C*HA1C')*X))    +PB2*real(trace((HB2C*HB2C')*X))+PA2*real(trace((HA2C*HA2C')*X))+sigma_c) - real(trace(get_grad_S(W, PA2,HAA, PB2, H12, PB1, HB1C, PA1, PB2, PA2, HA1C,    HB2C, HA2C)*(X-W))) <= t;

                -log(PB2*real(trace((H22*H22')*X))+PA1*real(trace((HAA*HAA')*X))+PB1*real(trace((H21*H21')*X))+sigma_ab+sigma_loop) -log(PA1*real(trace((HA1C*HA1C')*X))    +PB1*real(trace((HB1C*HB1C')*X))+PA2*real(trace((HA2C*HA2C')*X))+sigma_c) - real(trace(get_grad_S(W,  PA1,HAA, PB1, H21, PB2, HB2C, PA1, PB1, PA2, HA1C,   HB1C, HA2C)*(X-W))) <= t;

                -log(PA1*real(trace((H11*H11')*X))+PB2*real(trace((HBB*HBB')*X))+PA2*real(trace((H21*H21')*X))+sigma_ab+sigma_loop) -log(PB1*real(trace((HB1C*HB1C')*X))    +PB2*real(trace((HB2C*HB2C')*X))+PA2*real(trace((HA2C*HA2C')*X))+sigma_c) - real(trace(get_grad_S(W,  PB2,HBB, PA2, H21, PA1, HA1C, PB1, PB2, PA2, HB1C,   HB2C, HA2C)*(X-W))) <= t;

                -log(PA2*real(trace((H22*H22')*X))+PB1*real(trace((HBB*HBB')*X))+PA1*real(trace((H12*H12')*X))+sigma_ab+sigma_loop) -log(PB1*real(trace((HB1C*HB1C')*X))    +PB2*real(trace((HB2C*HB2C')*X))+PA1*real(trace((HA1C*HA1C')*X))+sigma_c) - real(trace(get_grad_S(W,  PB1,HBB, PA1, H12, PA2, HA2C, PA1, PB1, PB2, HA1C,   HB1C, HB2C)*(X-W))) <= t;

                diag(X) == 1;
                norm((X-W),1)<=2;

            cvx_end

            W = X;

            [U,D,U_] = svd(X);
            w0 = U(:,1)/U(L+1,1);
            w0 = w0./abs(w0);
            Wx = w0*w0';
    
            RA1 = get_secrecy_rate(Wx, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C);
            RA2 = get_secrecy_rate(Wx, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C);
            RB1 = get_secrecy_rate(Wx, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C);
            RB2 = get_secrecy_rate(Wx, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C);

            % Minimum rate at the end of each CVX iteration    
            min_rate = min([RA1,RA2,RB1,RB2])

            % Storing them in list

            A1_list = [A1_list, RA1];
            B1_list = [B1_list, RA2];
            A2_list = [A2_list, RB1];
            B2_list = [B2_list, RB2];

            min_list = [min_list, min_rate];

        end

    end

    A1_master_list = cat(1,A1_master_list,A1_list);
    B1_master_list = cat(1,B1_master_list,B1_list);
    A2_master_list = cat(1,A2_master_list,A2_list);
    B2_master_list = cat(1,B2_master_list,B2_list);
    min_master_list = cat(1,min_master_list,min_list);

end

close(f)

figure(1)
plot(iter_list,mean(min_master_list));
xlabel('Number of iterations')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')


% Create plots
figure(2)
t = tiledlayout(2,2);
nexttile
plot(iter_list,mean(A1_master_list));
xlabel('Number of iterations')
ylabel('A1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(iter_list,mean(A2_master_list));
xlabel('Number of iterations')
ylabel('A2 secrecy Rate(bits/sec/Hz)')
nexttile
plot(iter_list,mean(B1_master_list));
xlabel('Number of iterations')
ylabel('B1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(iter_list,mean(B2_master_list));
xlabel('Number of iterations')
ylabel('B2 secrecy Rate(bits/sec/Hz)')


function H = get_H(dti,dir,dtr,wave_l,L,Lo,pl,rng_val)

    gti = get_rayleigh_channel(dti, pl, L, 5*rng_val);
    hti = sqrt(Lo)*gti;
    gir = get_rayleigh_channel(dir, pl, L, 5*rng_val+2);
    hir = sqrt(Lo)*gir;
    gtr = get_rayleigh_direct_channel(dtr, pl, 5*rng_val+4);
    htr = sqrt(Lo)*gtr;
    H = cat(1,hti.*hir, htr);
end

function grad_S = get_grad_S(W, pinf1, inf1, pinf2, inf2, psec, sec, peves1, peves2, peves3, eves1, eves2, eves3)
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    sigma_c = 10^(-14.5);

    term1 = (pinf1*(inf1*inf1') + pinf2*(inf2*inf2'))/(get_A(W, pinf1, inf1, pinf2, inf2)+sigma_ab+sigma_loop);
    term2 = (psec*(sec*sec') + peves1*(eves1*eves1') + peves2*(eves2*eves2')+ peves3*(eves3*eves3'))/(get_B(W, psec, sec, peves1, peves2, peves3, eves1, eves2, eves3)+sigma_c);

    grad_S = (-1/log(2))*(term1 + term2)';
end

function channel_term = get_channel_term(power, W, H)
    H_tilda = H*H';

    channel_term = power*(real(trace(H_tilda*W)));
end

function A = get_A(W, pinf1, inf1, pinf2, inf2)

    A = get_channel_term(pinf1, W, inf1)+get_channel_term(pinf2, W, inf2);
end

function B = get_B(W, psec, sec, peves1, peves2, peves3, eves1, eves2, eves3)

    B = get_channel_term(psec, W, sec) + get_channel_term(peves1, W, eves1) + get_channel_term(peves2, W, eves2) + get_channel_term(peves3, W, eves3);
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

function c_rate = get_secrecy_rate(W, leg, inf1, inf2, eves, pleg, pinf1, pinf2, peves, peves1, peves2, peves3, eves1, eves2, eves3)

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    eves_inf = log2(1 + (get_channel_term(peves, W, eves))/(sigma_c + get_channel_term(peves1, W, eves1)+get_channel_term(peves2, W, eves2)+get_channel_term(peves3, W, eves3)));
    leg_inf = log2(1 + (get_channel_term(pleg, W, leg))/(get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop));
    c_rate = leg_inf - eves_inf;
%     c_rate = max(0,leg_inf - eves_inf);

end

function P = get_P(W, leg, inf1, inf2, pleg, pinf1, pinf2, peves1, peves2, peves3, eves1, eves2, eves3)

    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    sigma_c = 10^(-14.5);
    p1 = log(get_channel_term(pleg, W, leg)+get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop);
    p2 = log(get_channel_term(peves1, W, eves1) + get_channel_term(peves2, W, eves2) + get_channel_term(peves3, W, eves3) + sigma_c);
    P = p1 + p2;
    % P = log(get_channel_term(pleg, W, leg)+get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop);

end

function Q = get_Q(W,inf1, inf2, eves,pinf1, pinf2, peves, peves1, peves2, peves3, eves1, eves2, eves3)
    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    q1 = log(get_channel_term(pinf1, W, inf1)+get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop);
    q2 = log(get_channel_term(peves, W, eves) + get_channel_term(peves1, W, eves1) + get_channel_term(peves2, W, eves2) + get_channel_term(peves3, W, eves3)+sigma_c);
    Q = q1 + q2;
    % Q = log(get_channel_term(pinf1, W, inf1)+get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop)+log(get_channel_term(peves, W, eves) + sigma_c);

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