clear;
% Wavelength

lambda = 0.06;

% Random arrays for Gaussian Randomization
N = 10000;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
pl = 2;



L_list = 5:5:30;

A1_master_list = [];
B1_master_list = [];
min_master_list = [];

num_seeds = 10;
f = waitbar(0,'Please wait...');

for seed = 1:num_seeds

    if seed == 3
        continue
    end
    waitbar(seed/num_seeds,f,'Simulating...');

    A1_list = [];
    B1_list = [];
    min_list = [];

    % Number of elements of the IRS
    for L = 5:5:30

        % Rayleigh channels initialized

        % H11 = get_H(34.33,35.66,9.13,lambda,L,Lo,pl,3*seed);
        % HA1C = get_H(34.33,26.47,9.36,lambda,L,Lo,pl,3*seed+1);
        % HB1C = get_H(35.66,26.47,9.60,lambda,L,Lo,pl,3*seed+2);


        H11 = get_H(23.19,25.06,46.19,lambda,L,Lo,pl,3*seed);
        HA1C = get_H(23.19,6.76,22.50,lambda,L,Lo,pl,3*seed+1);
        HB1C = get_H(25.06,6.76,23.70,lambda,L,Lo,pl,3*seed+2);


        % IRS reflection matrix (Initialized randomly)
        wi = exp(1i*(2*rand(L,1)-1)*pi);
        w = [wi ; 1].';
        w_opt = w;

        PA1 = 0.04;
        PB1 = 0.04;
        % (W, leg, eves, pleg, peves, peves1, eves1)

        RA1 = get_secrecy_rate(w'*w, H11, HB1C, PB1, PB1, PA1, HA1C);
        RB1 = get_secrecy_rate(w'*w, H11, HA1C, PA1, PA1, PB1, HB1C);

        % Noise and residual loop interference
        sigma_ab = 10^(-14.5);
        sigma_loop = 10^(-14);
        sigma_c = 10^(-14.5);

        min_rate = 0;

        for j = 0:1

            max_min = 0;
            max_min_ind = 0;
            for i = 0:3
                power_string = num2str(dec2bin(i,2));
                PA1 = get_power(power_string(1));
                PB1 = get_power(power_string(2));
    % get_first(W, leg, pleg, peves1, eves1)
    % get_second(W, eves, peves, peves1, eves1)
                RA1 = get_first(w'*w, H11, PB1, PA1, HA1C) - get_second(w'*w, HB1C, PB1, PA1, HA1C);
                RB1 = get_first(w'*w, H11, PA1, PB1, HB1C) - get_second(w'*w, HA1C, PA1, PB1, HB1C);

                min_value = min([RA1,RB1]);

                if (min_value > max_min)

                    max_min = min_value;
                    max_min_ind = i;
                end

            end

            power_string = num2str(dec2bin(max_min_ind,2));
            PA1 = get_power(power_string(1));
            PB1 = get_power(power_string(2));

            t_old = 0;
            t_new = 1;

            for j = 0:100
    %%%%%%%%%%%%

                cvx_begin quiet
                cvx_solver Mosek
                variable w_opt(1,L+1) complex
                variable t
                maximize t
                subject to

                %A1
                (real(get_P(w'*w, PB1, H11)-get_Q(w'*w, PB1, HB1C, PA1, HA1C)) + ...
                    2*PB1*real(w_opt*H11*H11'*w')/(sigma_ab+sigma_loop+0.0000001) ...
                    -(get_channel_term(PB1, w'*w, H11)/((sigma_ab+sigma_loop+0.0000001)* ...
                    ((sigma_ab+sigma_loop+0.0000001)+get_channel_term(PB1, w'*w, H11))))*(sigma_ab+sigma_loop+PB1*pow_abs(w_opt * H11, 2)) ...
                    -(get_channel_term(PB1, w'*w, H11)/(sigma_ab+sigma_loop+0.0000001))...
                    - get_coef(w'*w, HB1C, HA1C, PB1, PA1)*PB1*quad_over_lin(w_opt * HB1C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w')))) >= t;

                %B1
                (real(get_P(w'*w, PA1, H11)-get_Q(w'*w, PA1, HA1C, PB1, HB1C)) + ...
                    2*PA1*real(w_opt*H11*H11'*w')/(sigma_ab+sigma_loop+0.0000001) ...
                    -(get_channel_term(PA1, w'*w, H11)/((sigma_ab+sigma_loop+0.0000001)* ...
                    ((sigma_ab+sigma_loop+0.0000001)+get_channel_term(PA1, w'*w, H11))))*(sigma_ab+sigma_loop+PA1*pow_abs(w_opt * H11, 2)) ...
                    -(get_channel_term(PA1, w'*w, H11)/(sigma_ab+sigma_loop+0.0000001))...
                    - get_coef(w'*w, HA1C, HB1C, PA1, PB1)*PA1*quad_over_lin(w_opt * HA1C,sigma_c+PB1*real((2*w_opt*HB1C-w*HB1C)*(HB1C'*w')))) >= t;

                % %A1
                % (log(1+ PB1*pow_abs(w_opt * H11, 2)/(sigma_ab+sigma_loop))-get_Q(w'*w, PB1, HB1C, PA1, HA1C) + ...
                %     - get_coef(w'*w, HB1C, HA1C, PB1, PA1)*PB1*quad_over_lin(w_opt * HB1C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w')))) >= t;

                % %B1
                % (log(1+PA1*pow_abs(w_opt * H11, 2)/(sigma_ab+sigma_loop))-get_Q(w'*w, PA1, HA1C, PB1, HB1C) + ...
                %     - get_coef(w'*w, HA1C, HB1C, PA1, PB1)*PA1*quad_over_lin(w_opt * HA1C,sigma_c+PB1*real((2*w_opt*HB1C-w*HB1C)*(HB1C'*w')))) >= t;

                for k = 1 : L
                    abs(w_opt(1,k)) <= 1;
                end
                w_opt(1,L+1)== 1;
                norm((w-w_opt),1)<=1;
                % t <= 1;
                % Add Taylor constraint here **** Important **** , update: Added

                cvx_end

                t_new = t
                if abs(t_new-t_old)<0.001
                    j
                    t_old = 0;
                    t_new = 1;
                    break
                end
                t_old = t_new;

    % P = get_P(w, pleg, leg)
    % Q = get_Q(w, peves, eves, peves1, eves1)

    %%%%%%%%%%%%
                w = w_opt; % Might need to take the other way around

                %% testing 
                % a = get_P(w'*w, PB1, H11)
                % a1 = -get_Q(w'*w, PB1, HB1C, PA1, HA1C)
                % b = 2*PB1*real(w_opt*H11*H11'*w')/(sigma_ab+sigma_loop+0.0000001)
                % c = -(get_channel_term(PB1, w'*w, H11)/((sigma_ab+sigma_loop+0.0000001)* ...
                % ((sigma_ab+sigma_loop+0.0000001)+get_channel_term(PB1, w'*w, H11))))*(sigma_ab+sigma_loop+PB1*pow_abs(w_opt * H11, 2))
                % d = -(get_channel_term(PB1, w'*w, H11)/(sigma_ab+sigma_loop+0.0000001))
                % e = - get_coef(w'*w, HB1C, HA1C, PB1, PA1)*PB1*quad_over_lin(w_opt * HB1C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w')))


                %% end testing




            end

            RA1 = get_secrecy_rate(w'*w, H11, HB1C, PB1, PB1, PA1, HA1C);
            RB1 = get_secrecy_rate(w'*w, H11, HA1C, PA1, PA1, PB1, HB1C);

            min_rate = min([RA1,RB1]);

        end


        RA1 = get_secrecy_rate(w'*w, H11, HB1C, PB1, PB1, PA1, HA1C)
        RB1 = get_secrecy_rate(w'*w, H11, HA1C, PA1, PA1, PB1, HB1C)


        A1_list = [A1_list, RA1];
        B1_list = [B1_list, RB1];

        min_rate = min([RA1,RB1])
        min_list = [min_list, min_rate];

    end

    A1_master_list = cat(1,A1_master_list,A1_list);
    B1_master_list = cat(1,B1_master_list,B1_list);
    min_master_list = cat(1,A1_master_list,A1_list);

end

close(f)

figure(1)
plot(L_list,min_list);
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')


% Create plots
figure(2)
t = tiledlayout(1,2);
nexttile
plot(L_list,A1_list);
xlabel('Number of IRS elements')
ylabel('A1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,B1_list);
xlabel('Number of IRS elements')
ylabel('B1 secrecy Rate(bits/sec/Hz)')

function H = get_H(dti,dir,dtr,wave_l,L,Lo,pl,rng_val)

    gti = get_rayleigh_channel(dti, pl, L, 5*rng_val);
    hti = sqrt(Lo)*gti;
    gir = get_rayleigh_channel(dir, pl, L, 5*rng_val+2);
    hir = sqrt(Lo)*gir;
    gtr = get_rayleigh_direct_channel(dtr, pl, 5*rng_val+4);
    htr = sqrt(Lo)*gtr;
    H = cat(1,hti.*hir, htr);
end

function channel_term = get_channel_term(power, W, H)
    H_tilda = H*H';

    channel_term = power*(real(trace(H_tilda*W)));
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

function c_rate = get_secrecy_rate(W, leg, eves, pleg, peves, peves1, eves1)

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    eves_inf = log2(1 + (get_channel_term(peves, W, eves))/(sigma_c + get_channel_term(peves1, W, eves1)));
    leg_inf = log2(1 + (get_channel_term(pleg, W, leg))/(sigma_ab + sigma_loop));
    c_rate = leg_inf - eves_inf;
%     c_rate = max(0,leg_inf - eves_inf);

end

function P = get_P(W, pleg, leg)

    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    P = log2(1 + (get_channel_term(pleg, W, leg))/(sigma_ab+sigma_loop));

end

function Q = get_Q(W, peves, eves, peves1, eves1)

    sigma_c = 10^(-14.5);
    Q = log2(1 + (get_channel_term(peves, W, eves))/(get_channel_term(peves1, W, eves1)+sigma_c));

end

function first = get_first(W, leg, pleg, peves1, eves1)

    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    sigma_c = 10^(-14.5);
    p1 = log(get_channel_term(pleg, W, leg)+ sigma_ab + sigma_loop);
    p2 = log(get_channel_term(peves1, W, eves1) + sigma_c);
    first = p1 + p2;
    % P = log(get_channel_term(pleg, W, leg)+get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop);

end

function second = get_second(W, eves, peves, peves1, eves1)
    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    q1 = log(sigma_ab + sigma_loop);
    q2 = log(get_channel_term(peves, W, eves) + get_channel_term(peves1, W, eves1)+sigma_c);
    second = q1 + q2;
    % Q = log(get_channel_term(pinf1, W, inf1)+get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop)+log(get_channel_term(peves, W, eves) + sigma_c);

end

function coef = get_coef(W,eves,eves1,peves, peves1)
    sigma_c = 10^(-14.5);
    value = (1 + (get_channel_term(peves,W,eves))/(sigma_c + get_channel_term(peves1,W,eves1)));
    coef = real(1/value);

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
