clear;
close all;
load('L40channels.mat','Hab','Hac','Hbc','Ha2b2','Ha2c','Hb2c','Ha2a','Ha2b','Hb2a','Hb2b','L','N','sigma_ab','sigma_c','sigma_loop');

Pmax = dbm2watt(15);
Pmin  = dbm2watt(0);

% Initialize IRS elements
w = exp(1i*rand(L,1)*(pi));
wrand = w;
w_sdp = [wrand; 1];
r = 16;

Pa = Pmin;
Pb = Pmin;
Pa2 = Pmin;
Pb2 = Pmin;

[minrate_rand,~,~] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w_sdp(1:end-1),sigma_ab,sigma_loop,sigma_c);

I1 = 12;

minrate_rand

tic
for i1 = 1:I1
    [Pa, Pb, Pa2, Pb2] = power_optimization_maxmin_iter(Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,w_sdp,5);
    w_sdp = irs_optimization_maxmin_iter(L, Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,w_sdp,3);
    [minrate_both,~,~] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w_sdp(1:end-1),sigma_ab,sigma_loop,sigma_c)
end
[minrate_both,~,~] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w_sdp(1:end-1),sigma_ab,sigma_loop,sigma_c);
toc
minrate_both
exhtrials = 22;
Pavals = linspace(Pmin, Pmax, exhtrials);
Pbvals = linspace(Pmin, Pmax, exhtrials);
Pa2vals = linspace(Pmin, Pmax, exhtrials);
Pb2vals = linspace(Pmin, Pmax, exhtrials);
F = exhtrials;

optvalueexh = 0;
optPa = 0;
optPb = 0;
optPa2 = 0;
optPb2 = 0;
optw = 0;
minvalues = zeros(exhtrials,exhtrials,exhtrials,exhtrials,exhtrials);
for i = 1:length(Pavals)
    i
    for j = 1:length(Pbvals)
        for k = 1:length(Pa2vals)
            for l = 1:length(Pb2vals)
                for a = 1:F
                    w(1) = exp(1i*(a/F)*(2*pi));
                    [minrate,~,~] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pavals(i),Pbvals(j),Pa2vals(k),Pb2vals(l),w,sigma_ab,sigma_loop,sigma_c);
                    if minrate > optvalueexh
                        optPa =  Pavals(i);
                        optPb =  Pbvals(j);
                        optPa2 =  Pa2vals(k);
                        optPb2 =  Pb2vals(l);
                        optw = exp(1i*(a/F)*(2*pi)); 
                        optw;
                    end
                    optvalueexh = max(optvalueexh,minrate);
                    minvalues(i,j,k,l,a) = minrate;

                end
                
            end
        end
    end
end

optvalueexh
minrate_both

F = 1000;
optimumvalueexh = 0;
for a = 1:F
    w(1) = exp(1i*(a/F)*(2*pi));
    [minrate,~,~] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa,Pb,Pa2,Pb2,w,sigma_ab,sigma_loop,sigma_c);
    optimumvalueexh = max(optimumvalueexh,minrate);
end
optimumvalueexh
minrate_sdp

[Rp1, Rp2] = power_optimization_maxmin_iter_test(Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,[w;1]);

tic
[Pa_iter, Pb_iter, Pa2_iter, Pb2_iter] = power_optimization_maxmin_iter(Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,[w;1]);
toc

[minrate_iter,p1rate_iter,p2rate_iter] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa_iter,Pb_iter,Pa2_iter,Pb2_iter,w,sigma_ab,sigma_loop,sigma_c);


function [minrate,p1rate,p2rate] = min_secrecy_rate(Hab,Hac,Hbc,Ha2b2,Ha2c,Hb2c,Ha2a,Ha2b,Hb2a,Hb2b,Pa,Pb,Pa2,Pb2,w,sigma_ab,sigma_loop,sigma_c)
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
end

function power = dbm2watt(dbm_value)
    % Converting power values from dBm to Watt for calculations
    power = (10^(dbm_value/10))*(10^(-3));
end

% w = zeros(L,1);

% allmaxsumrates = zeros(N,1);
% F = 10000;
% I1 = 4;
% I2 = 4;

% Patrial = Pmin;
% Pbtrial = Pmin;
% 
% [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,w);
% 
% w_hat = [w;1];
% for zz = 1:7
%     zz
%     y = sqrt(sigma_ab + sigma_loop + Pbtrial*abs(w_hat'*Hab(:,r))^2)/(Patrial*abs(w_hat'*Hac(:,r))^2+Pbtrial*abs(w_hat'*Hbc(:,r))^2+sigma_c);
%     [Patrial, Pbtrial] = power_optimization_iter(y,Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,[w;1]);
% end
% Patrial
% Pbtrial

% function [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab,Hac,Hbc,sigma_ab,sigma_loop,sigma_c,w)
%     % Using bisection on grad to obtain the max point
%     pa_low = Pmin;
%     pa_high = Pmax;
%     while abs(pa_low - pa_high) > 1e-10
%         pa = (pa_high+pa_low)/2;
% 
%         gradvalue = calgrad(Hab,Hac,Hbc,w,sigma_ab,sigma_loop,sigma_c,Pmax,pa);
% 
%         if gradvalue < 0
%             pa_high = pa;
%         else
%             pa_low = pa; 
%         end
%     end 
%     Pasol = pa;
%     Pbsol = Pmax - Pasol;
% 
%     optimalvalrate = sum_secrecy_rate(Hab,Hac,Hbc,Pasol,Pbsol,w,sigma_ab,sigma_loop,sigma_c);
% 
%     corner1 = sum_secrecy_rate(Hab,Hac,Hbc,Pmin,Pmin,w,sigma_ab,sigma_loop,sigma_c);
%     corner2 = sum_secrecy_rate(Hab,Hac,Hbc,Pmin,Pmax,w,sigma_ab,sigma_loop,sigma_c);
%     corner3 = sum_secrecy_rate(Hab,Hac,Hbc,Pmax,Pmin,w,sigma_ab,sigma_loop,sigma_c);
% 
%     % Optimizing between the corner points and the maximum point on Pa + Pb =
%     % Pmax line
%     [optimum,index]= max([corner1 corner2 corner3 optimalvalrate]);
% 
%     switch index
%         case 1
%             % disp(strcat('Pa = ',num2str(Pamin) ,' Pb = ',num2str(Pbmin) ))
%             Paopt = Pmin;
%             Pbopt = Pmin;
%         case 2
%             % disp(strcat('Pa = ',num2str(Pamin) ,' Pb = ',num2str(Pbmax) ))
%             Paopt = Pmin;
%             Pbopt = Pmax;
%         case 3
%             % disp(strcat('Pa = ',num2str(Pamax) ,' Pb = ',num2str(Pbmin) ))
%             Paopt = Pmax;
%             Pbopt = Pmin;
%         otherwise
%             % disp(strcat('Pa = ',num2str(Pasol) ,' Pb = ',num2str(Pbsol) ))
%             Paopt = Pasol;
%             Pbopt = Pbsol;
%     end
% end

% function grad = calgrad(Hab,Hac,Hbc,w,sigma_ab,sigma_loop,sigma_c,Pmax,pa)
% 
%     wx = [w;1];
% 
%     x = real(trace(Hab*Hab'*wx*wx'));
%     y = real(trace(Hac*Hac'*wx*wx'));
%     z = real(trace(Hbc*Hbc'*wx*wx'));
% 
%     % Both gradients are correct
%     grad = x*(sigma_ab + sigma_loop + (Pmax - pa)*x)*(sigma_c + Pmax*z + (y-z)*pa)...
%         - x*(sigma_ab + sigma_loop+pa*x)*(sigma_c + Pmax*z + (y-z)*pa)...
%         - (y-z)*(sigma_ab + sigma_loop+pa*x)*(sigma_ab + sigma_loop + (Pmax - pa)*x);
% end

% function sumrate = sum_secrecy_rate(Hab,Hac,Hbc,Pa,Pb,w,sigma_ab,sigma_loop,sigma_c)
%     wx = [w;1];
% 
%     I_b_A = log2(1 + (Pb*abs(wx'*Hab)^2)/(sigma_ab+sigma_loop));
%     I_a_B = log2(1 + (Pa*abs(wx'*Hab)^2)/(sigma_ab+sigma_loop));
%     I_ab_C = log2(1 + (Pa*abs(wx'*Hac)^2 + Pb*abs(wx'*Hbc)^2)/(sigma_c));
%     sumrate = max(I_b_A+I_a_B-I_ab_C,0);
% end

