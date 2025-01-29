clear;
close all;
load('L20channels.mat','Hab','Hac','Hbc','Ha2b2','Ha2c','Hb2c','Ha2a','Ha2b','Hb2a','Hb2b','L','N','sigma_ab','sigma_c','sigma_loop');

Pmax = dbm2watt(15);
Pmin  = dbm2watt(0);

% Initialize IRS elements
w = exp(1i*rand(L,1)*(pi));
wrand = w;

num_iterations = 40;
exhtrials = 50;
numrealizations = 1;

Pavals = linspace(Pmin, Pmax, exhtrials);
Pbvals = linspace(Pmin, Pmax, exhtrials);
Pa2vals = linspace(Pmin, Pmax, exhtrials);
Pb2vals = linspace(Pmin, Pmax, exhtrials);

optvalueexh = zeros(1,numrealizations);
optPa = zeros(1,numrealizations);
optPb = zeros(1,numrealizations);
optPa2 = zeros(1,numrealizations);
optPb2 = zeros(1,numrealizations);
minvalues = zeros(exhtrials,exhtrials,exhtrials,exhtrials);
obj_evolutionset = zeros(num_iterations,numrealizations);
Pa_evolutionset = zeros(num_iterations,numrealizations);
Pb_evolutionset = zeros(num_iterations,numrealizations);
Pa2_evolutionset = zeros(num_iterations,numrealizations);
Pb2_evolutionset = zeros(num_iterations,numrealizations);

% r = 10;

for r = 1:numrealizations

    Pa = Pmax;
    Pb = Pmax;
    Pa2 = Pmax;
    Pb2 = Pmax;

    tic
    [Pa_iter, Pb_iter, Pa2_iter, Pb2_iter, Pa_evolution, Pb_evolution, Pa2_evolution, Pb2_evolution, obj_evolution] = power_optimization_maxmin_iter_plot(Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),sigma_ab,sigma_loop,sigma_c,[w;1],num_iterations);
    toc
    obj_evolutionset(:,r) = obj_evolution;
    Pa_evolutionset(:,r) = Pa_evolution;
    Pb_evolutionset(:,r) = Pb_evolution;
    Pa2_evolutionset(:,r) = Pa2_evolution;
    Pb2_evolutionset(:,r) = Pb2_evolution;

    [minrate_iter,p1rate_iter,p2rate_iter] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pa_iter,Pb_iter,Pa2_iter,Pb2_iter,w,sigma_ab,sigma_loop,sigma_c);

    tic
    for i = 1:length(Pavals)
        i
        for j = 1:length(Pbvals)
            for k = 1:length(Pa2vals)
                for l = 1:length(Pb2vals)
                    [minrate,p1rate,p2rate] = min_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Ha2b2(:,r),Ha2c(:,r),Hb2c(:,r),Ha2a(:,r),Ha2b(:,r),Hb2a(:,r),Hb2b(:,r),Pavals(i),Pbvals(j),Pa2vals(k),Pb2vals(l),w,sigma_ab,sigma_loop,sigma_c);
                    optvalueexh(r) = max(optvalueexh(r),minrate);
                    if minrate == optvalueexh(r)
                        optPa(r) =  Pavals(i);
                        optPb(r) =  Pbvals(j);
                        optPa2(r) =  Pa2vals(k);
                        optPb2(r) =  Pb2vals(l);
                    end
                    minvalues(i,j,k,l) = minrate;
                end
            end
        end
    end
    toc
    optvalueexh(r)

end

% Create the first figure for Paopt and Pbopt
figure;
yline(optPa(r)/optPb2(r), 'r--', 'optPa/optPb2', 'LabelHorizontalAlignment', 'left');
hold on;
yline(optPb(r)/optPb2(r), 'g--', 'optPb/optPb2', 'LabelHorizontalAlignment', 'left');
yline(optPa2(r)/optPb2(r), 'g--', 'optPa2/optPb2', 'LabelHorizontalAlignment', 'left');
% yline(optPb2/optPb2, 'g--', 'Pb2opt', 'LabelHorizontalAlignment', 'left');
plot(1:num_iterations, Pa_evolutionset(:,r)./Pb2_evolutionset(:,r), 'r-o');
plot(1:num_iterations, Pb_evolutionset(:,r)./Pb2_evolutionset(:,r), 'g-o');
plot(1:num_iterations, Pa2_evolutionset(:,r)./Pb2_evolutionset(:,r), 'm-o');
% plot(1:num_iterations, Pb2_evolution, 'b-o');
% ylim([Pmin, Pmax]);
xlim([0, num_iterations+1]);
xlabel('X-axis');
ylabel('Power Values');
title('Optimization Results: Paopt and Pbopt');
% legend('Paopt', 'Pbopt', 'Location', 'best');
grid on;
hold off;



% fractional_increases = zeros(num_iterations-1,numrealizations);
% diffwithexh = zeros(num_iterations,numrealizations);
% for r = 1:numrealizations
%     fractional_increases(:,r) = diff(obj_evolutionset(:,r))./obj_evolutionset(1:end-1,r);
%     diffwithexh(:,r) = abs(obj_evolutionset(:,r) - optvalueexh(r))/optvalueexh(r);
% end
% 
% save('converge_power')
% 
% figure;
% semilogy(1:num_iterations-1, mean(fractional_increases,2))
% hold on;
% semilogy(1:num_iterations, mean(diffwithexh,2));

% Create the second figure for optimum
figure;
yline(mean(optvalueexh(r)), 'b--', 'Optimum', 'LabelHorizontalAlignment', 'left');hold on;
plot(1:num_iterations, mean(obj_evolutionset(:,r),2), 'b-o');
ylim([min(obj_evolutionset,[],"all") - 1, max(obj_evolutionset,[],"all") + 1]);
xlim([0, num_iterations+1]);
xlabel('X-axis');
ylabel('Power Value');
title('Optimization Result: Optimum');
% legend('Optimum', 'Location', 'best');
grid on;
hold off;

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

