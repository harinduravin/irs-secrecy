clear;
close all;
load('L5channels.mat','Hab','Hac','Hbc','Ha2b2','Ha2c','Hb2c','Ha2a','Ha2b','Hb2a','Hb2b','L','N','sigma_ab','sigma_c','sigma_loop');

Pmax = dbm2watt(15);
Pmin  = dbm2watt(0);

I1 = 3;
R  = 100;
convergence_data = zeros(R,I1);
store_data = zeros(R,9);
rand_data = zeros(R,1);

for r = 1:R

    fprintf(['Attempt: ',num2str(r),'\n'])

    % Initialize IRS elements
    w = exp(1i*rand(L,1)*(pi));
    wrand = w;
    w_sdp = [wrand; 1];

    Pa = Pmin;
    Pb = Pmin;
    Pa2 = Pmin;
    Pb2 = Pmin;
    
    output_rand = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w_sdp(1:end-1),sigma_ab,sigma_loop,sigma_c);
    rand_data(r) = output_rand(1);
    % output_rand
    tic
    for i1 = 1:I1
        [Pa, Pb, Pa2, Pb2] = power_optimization_maxmin_iter(Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,w_sdp,25);
        fprintf('|')
        w_sdp = irs_optimization_maxmin_iter(L, Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,w_sdp,8);
        fprintf('|')
        output_convergence = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w_sdp(1:end-1),sigma_ab,sigma_loop,sigma_c);
        convergence_data(r,i1) = output_convergence(1);
        % output_convergence
    end
    store_data(r,:) = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w_sdp(1:end-1),sigma_ab,sigma_loop,sigma_c);
    toc
    fprintf(['minvalue: ',num2str(store_data(r,1))])
    fprintf(['\n','randvalue: ',num2str(rand_data(r))])
    fprintf('\n')
end

save('data1')

function funcoutput = min_secrecy_rate(Hab,Hac,Hbc,Ha2b2,Ha2c,Hb2c,Ha2a,Ha2b,Hb2a,Hb2b,Pa,Pb,Pa2,Pb2,w,sigma_ab,sigma_loop,sigma_c)
    wx = [w;1];

    I_b_A = log2(1 + (Pb*abs(wx'*Hab)^2)/(sigma_ab+sigma_loop+Pb2*abs(wx'*Hb2a)^2+Pa2*abs(wx'*Ha2a)^2));
    I_a_B = log2(1 + (Pa*abs(wx'*Hab)^2)/(sigma_ab+sigma_loop+Pb2*abs(wx'*Hb2b)^2+Pa2*abs(wx'*Ha2b)^2));
    I_ab_C = log2(1 + (Pa*abs(wx'*Hac)^2 + Pb*abs(wx'*Hbc)^2)/(sigma_c+Pb2*abs(wx'*Hb2c)^2+Pa2*abs(wx'*Ha2c)^2));
    p1rate = max(I_b_A+I_a_B-I_ab_C,0);

    I2_b_A = log2(1 + (Pb2*abs(wx'*Ha2b2)^2)/(sigma_ab+sigma_loop+Pb*abs(wx'*Ha2b)^2+Pa*abs(wx'*Ha2a)^2));
    I2_a_B = log2(1 + (Pa2*abs(wx'*Ha2b2)^2)/(sigma_ab+sigma_loop+Pb*abs(wx'*Hb2b)^2+Pa*abs(wx'*Hb2a)^2));
    I2_ab_C = log2(1 + (Pa2*abs(wx'*Ha2c)^2 + Pb2*abs(wx'*Hb2c)^2)/(sigma_c+Pb*abs(wx'*Hbc)^2+Pa*abs(wx'*Hac)^2));
    p2rate = max(I2_b_A+I2_a_B-I2_ab_C,0);

    minrate = min(p1rate, p2rate);
    funcoutput = [minrate,p1rate,p2rate,I_b_A,I_a_B,I_ab_C,I2_b_A,I2_a_B,I2_ab_C];
end

function power = dbm2watt(dbm_value)
    % Converting power values from dBm to Watt for calculations
    power = (10^(dbm_value/10))*(10^(-3));
end

