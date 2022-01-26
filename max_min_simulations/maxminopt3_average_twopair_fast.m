clear;
% Wavelength

lambda = 0.06;

% Random arrays for Gaussian Randomization
N = 10000;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
pl = 2;

maximum_L = 20;

L_list = 5:5:maximum_L;

A1_master_list = [];
B1_master_list = [];
A2_master_list = [];
B2_master_list = [];
min_master_list = [];

num_seeds = 8;
f = waitbar(0,'Please wait...');

for seed = 4:num_seeds
    if seed==2
        continue
    end
    waitbar(seed/num_seeds,f,'Simulating...');

    A1_list = [];
    B1_list = [];
    A2_list = [];
    B2_list = [];
    min_list = [];

    % Number of elements of the IRS
    for L = 5:5:maximum_L

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
        w_opt = w;

        PA1 = 0.04;
        PB1 = 0.04;
        PA2 = 0.04;
        PB2 = 0.04;

        RA1 = get_secrecy_rate(w'*w, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C)
        RA2 = get_secrecy_rate(w'*w, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C)
        RB1 = get_secrecy_rate(w'*w, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C)
        RB2 = get_secrecy_rate(w'*w, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C)

        % Noise and residual loop interference
        sigma_ab = 10^(-14.5);
        sigma_loop = 10^(-14);
        sigma_c = 10^(-14.5);

        min_rate = 0;


        
        for j = 0:1

            max_min = 0;
            max_min_ind = 0;
            for i = 0:15
                power_string = num2str(dec2bin(i,4));
                PA1 = get_power(power_string(1));
                PA2 = get_power(power_string(2));
                PB1 = get_power(power_string(3));
                PB2 = get_power(power_string(4));

                RA1 = get_P(w'*w, H11, HAA, H12, PB1, PA2, PB2, PA1, PB2, PA2, HA1C, HB2C, HA2C) - get_Q(w'*w, HAA, H12,HB1C, PA2, PB2 ,PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C);
                RA2 = get_P(w'*w, H22, HAA, H21, PB2, PA1, PB1, PA1, PB1, PA2, HA1C, HB1C, HA2C) - get_Q(w'*w, HAA, H21,HB2C, PA1, PB1 ,PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C);
                RB1 = get_P(w'*w, H11, HBB, H21, PA1, PB2, PA2, PB1, PB2, PA2, HB1C, HB2C, HA2C) - get_Q(w'*w, HBB, H21,HA1C, PB2, PA2 ,PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C);
                RB2 = get_P(w'*w, H22, HBB, H12, PA2, PB1, PA1, PA1, PB1, PB2, HA1C, HB1C, HB2C) - get_Q(w'*w, HBB, H12,HA2C, PB1, PA1 ,PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C);

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

            t_old = 0;
            t_new = 1;

            for j = 0:30

                cvx_begin
                cvx_solver Mosek
                variable w_opt(1,L+1) complex
                variable t
                maximize t
                subject to

                %A1
                (get_P_Q(w'*w, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C) + ...
                    2*PB1*real(w_opt*H11*H11'*w')/(sigma_ab+sigma_loop+get_A(w'*w, PA2, HAA, PB2, H12)) ...
                    -(get_channel_term(PB1, w'*w, H11)/((sigma_ab+sigma_loop+get_A(w'*w, PA2, HAA, PB2, H12))* ...
                    ((sigma_ab+sigma_loop+get_A(w'*w, PA2, HAA, PB2, H12))+get_channel_term(PB1, w'*w, H11))))*(sigma_ab+sigma_loop+PB1*pow_abs(w_opt * H11, 2)+PA2*pow_abs(w_opt * HAA, 2)+PB2*pow_abs(w_opt * H12, 2)) ...
                    -(get_channel_term(PB1, w'*w, H11)/(sigma_ab+sigma_loop+get_A(w'*w, PA2, HAA, PB2, H12)))...
                    - (1/(1+get_zbar(w'*w, HB1C, HA1C, HB2C, HA2C, PB1, PA1, PB2, PA2)))*(PB1*quad_over_lin(w_opt * HB1C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w')))-get_zbar(w'*w, HB1C, HA1C, HB2C, HA2C, PB1, PA1, PB2, PA2))) >= t;
                    % - (1/(1+get_zbar(w'*w, HB1C, HA1C, HB2C, HA2C, PB1, PA1, PB2, PA2)))*(PB1*quad_over_lin(w_opt * HB1C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w'))+PB2*real((2*w_opt*HB2C-w*HB2C)*(HB2C'*w'))+PA2*real((2*w_opt*HA2C-w*HA2C)*(HA2C'*w')))-get_zbar(w'*w, HB1C, HA1C, HB2C, HA2C, PB1, PA1, PB2, PA2))) >= t;

                %A2
                (get_P_Q(w'*w, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C) + ...
                    2*PB2*real(w_opt*H22*H22'*w')/(sigma_ab+sigma_loop+get_A(w'*w, PA1, HAA, PB1, H21)) ...
                    -(get_channel_term(PB2, w'*w, H22)/((sigma_ab+sigma_loop+get_A(w'*w, PA1, HAA, PB1, H21))* ...
                    ((sigma_ab+sigma_loop+get_A(w'*w, PA1, HAA, PB1, H21))+get_channel_term(PB2, w'*w, H22))))*(sigma_ab+sigma_loop+PB2*pow_abs(w_opt * H22, 2)+PA1*pow_abs(w_opt * HAA, 2)+PB1*pow_abs(w_opt * H21, 2)) ...
                    -(get_channel_term(PB2, w'*w, H22)/(sigma_ab+sigma_loop+get_A(w'*w, PA1, HAA, PB1, H21)))...
                    - (1/(1+get_zbar(w'*w, HB2C, HA1C, HB1C, HA2C, PB2, PA1, PB1, PA2)))*(PB2*quad_over_lin(w_opt * HB2C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w')))-get_zbar(w'*w, HB2C, HA1C, HB1C, HA2C, PB2, PA1, PB1, PA2))) >= t;
                    % - (1/(1+get_zbar(w'*w, HB2C, HA1C, HB1C, HA2C, PB2, PA1, PB1, PA2)))*(PB2*quad_over_lin(w_opt * HB2C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w'))+PB1*real((2*w_opt*HB1C-w*HB1C)*(HB1C'*w'))+PA2*real((2*w_opt*HA2C-w*HA2C)*(HA2C'*w')))-get_zbar(w'*w, HB2C, HA1C, HB1C, HA2C, PB2, PA1, PB1, PA2))) >= t;

                %B1
                (get_P_Q(w'*w, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C) + ...
                    2*PA1*real(w_opt*H11*H11'*w')/(sigma_ab+sigma_loop+get_A(w'*w, PA2, H21, PB2, HBB)) ...
                    -(get_channel_term(PA1, w'*w, H11)/((sigma_ab+sigma_loop+get_A(w'*w, PA2, H21, PB2, HBB))* ...
                    ((sigma_ab+sigma_loop+get_A(w'*w, PA2, H21, PB2, HBB))+get_channel_term(PA1, w'*w, H11))))*(sigma_ab+sigma_loop+PA1*pow_abs(w_opt * H11, 2)+PA2*pow_abs(w_opt * H21, 2)+PB2*pow_abs(w_opt * HBB, 2)) ...
                    -(get_channel_term(PA1, w'*w, H11)/(sigma_ab+sigma_loop+get_A(w'*w, PA2, H21, PB2, HBB)))...
                    - (1/(1+get_zbar(w'*w, HA1C, HB1C, HB2C, HA2C, PA1, PB1, PB2, PA2)))*(PA1*quad_over_lin(w_opt * HA1C,sigma_c+PB1*real((2*w_opt*HB1C-w*HB1C)*(HB1C'*w')))-get_zbar(w'*w, HA1C, HB1C, HB2C, HA2C, PA1, PB1, PB2, PA2))) >= t;
                    % - (1/(1+get_zbar(w'*w, HA1C, HB1C, HB2C, HA2C, PA1, PB1, PB2, PA2)))*(PA1*quad_over_lin(w_opt * HA1C,sigma_c+PB1*real((2*w_opt*HB1C-w*HB1C)*(HB1C'*w'))+PB2*real((2*w_opt*HB2C-w*HB2C)*(HB2C'*w'))+PA2*real((2*w_opt*HA2C-w*HA2C)*(HA2C'*w')))-get_zbar(w'*w, HA1C, HB1C, HB2C, HA2C, PA1, PB1, PB2, PA2))) >= t;

                %B2
                (get_P_Q(w'*w, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C) + ...
                    2*PA2*real(w_opt*H22*H22'*w')/(sigma_ab+sigma_loop+get_A(w'*w, PB1, HBB, PA1, H12)) ...
                    -(get_channel_term(PA2, w'*w, H22)/((sigma_ab+sigma_loop+get_A(w'*w, PB1, HBB, PA1, H12))* ...
                    ((sigma_ab+sigma_loop+get_A(w'*w, PB1, HBB, PA1, H12))+get_channel_term(PA2, w'*w, H22))))*(sigma_ab+sigma_loop+PA2*pow_abs(w_opt * H22, 2)+PB1*pow_abs(w_opt * HBB, 2)+PA1*pow_abs(w_opt * H12, 2)) ...
                    -(get_channel_term(PA2, w'*w, H22)/(sigma_ab+sigma_loop+get_A(w'*w, PB1, HBB, PA1, H12)))...
                    - (1/(1+get_zbar(w'*w, HA2C, HA1C, HB2C, HB1C, PA2, PA1, PB2, PB1)))*(PA2*quad_over_lin(w_opt * HA2C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w')))-get_zbar(w'*w, HA2C, HA1C, HB2C, HB1C, PA2, PA1, PB2, PB1))) >= t;
                    % - (1/(1+get_zbar(w'*w, HA2C, HA1C, HB2C, HB1C, PA2, PA1, PB2, PB1)))*(PA2*quad_over_lin(w_opt * HA2C,sigma_c+PA1*real((2*w_opt*HA1C-w*HA1C)*(HA1C'*w'))+PB2*real((2*w_opt*HB2C-w*HB2C)*(HB2C'*w'))+PB1*real((2*w_opt*HB1C-w*HB1C)*(HB1C'*w')))-get_zbar(w'*w, HA2C, HA1C, HB2C, HB1C, PA2, PA1, PB2, PB1))) >= t;

%% Above 4 parts are changed for testing purposes. uncomment last three lines of each and comments last two lines to retrieve the original lines.

                for k = 1 : L+1
                    abs(w_opt(1,k)) <= 1;
                end
                % w_opt(1,L+1)== 1;
                norm((w_opt-w),1)<=0.1;

                cvx_end

                t_new = t
                if abs(t_new-t_old)/t_old<0.01
                    j
                    t_old = 0;
                    t_new = 1;
                    break
                end
                t_old = t_new;

    %%%%%%%%%%%%
                w = w_opt; % Might need to take the other way around
                RA1 = get_secrecy_rate(w'*w, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C);
                RA2 = get_secrecy_rate(w'*w, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C);
                RB1 = get_secrecy_rate(w'*w, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C);
                RB2 = get_secrecy_rate(w'*w, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C);
    
    
                min_rate = min([RA1,RA2,RB1,RB2])

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
    %         fprintf(string(t))
    %         fprintf('\r')

            RA1 = get_secrecy_rate(w'*w, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C);
            RA2 = get_secrecy_rate(w'*w, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C);
            RB1 = get_secrecy_rate(w'*w, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C);
            RB2 = get_secrecy_rate(w'*w, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C);


            min_rate = min([RA1,RA2,RB1,RB2]);

        end

        RA1 = get_secrecy_rate(w'*w, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1, PA1, PB2, PA2, HA1C, HB2C, HA2C);
        RA2 = get_secrecy_rate(w'*w, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2, PA1, PB1, PA2, HA1C, HB1C, HA2C);
        RB1 = get_secrecy_rate(w'*w, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1, PB1, PB2, PA2, HB1C, HB2C, HA2C);
        RB2 = get_secrecy_rate(w'*w, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2, PA1, PB1, PB2, HA1C, HB1C, HB2C);

        A1_list = [A1_list, RA1];
        B1_list = [B1_list, RA2];
        A2_list = [A2_list, RB1];
        B2_list = [B2_list, RB2];

        min_rate = min([RA1,RA2,RB1,RB2]);
        min_list = [min_list, min_rate];

    end

    A1_master_list = cat(1,A1_master_list,A1_list);
    B1_master_list = cat(1,B1_master_list,B1_list);
    A2_master_list = cat(1,A2_master_list,A2_list);
    B2_master_list = cat(1,B2_master_list,B2_list);
    min_master_list = cat(1,A1_master_list,A1_list);

end

close(f)

figure(1)
plot(L_list,mean(min_master_list));
xlabel('Number of IRS elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')


% Create plots
figure(2)
t = tiledlayout(2,2);
nexttile
plot(L_list,mean(A1_master_list));
xlabel('Number of IRS elements')
ylabel('A1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,mean(A2_master_list));
xlabel('Number of IRS elements')
ylabel('A2 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,mean(B1_master_list));
xlabel('Number of IRS elements')
ylabel('B1 secrecy Rate(bits/sec/Hz)')
nexttile
plot(L_list,mean(B2_master_list));
xlabel('Number of IRS elements')
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

function c_rate = get_P_Q(W, leg, inf1, inf2, eves, pleg, pinf1, pinf2, peves, peves1, peves2, peves3, eves1, eves2, eves3)

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    eves_inf = log(1 + (get_channel_term(peves, W, eves))/(sigma_c + get_channel_term(peves1, W, eves1)+get_channel_term(peves2, W, eves2)+get_channel_term(peves3, W, eves3)));
    leg_inf = log(1 + (get_channel_term(pleg, W, leg))/(get_channel_term(pinf1, W, inf1) + get_channel_term(pinf2, W, inf2) + sigma_ab + sigma_loop));
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

function value = get_zbar(W,eves,eves1,eves2,eves3,peves, peves1, peves2, peves3)
    sigma_c = 10^(-14.5);
    value = (get_channel_term(peves,W,eves))/(sigma_c + get_channel_term(peves1,W,eves1) + get_channel_term(peves2,W,eves2)+ get_channel_term(peves3,W,eves3));

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
