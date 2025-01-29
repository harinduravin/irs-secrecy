clear;
close all;
load('L20channels.mat','Hab','Hac','Hbc','L','N','sigma_ab','sigma_c','sigma_loop');

Pmax = dbm2watt(15);
Pmin  = dbm2watt(0);

F = 10000;
I1 = 4;
I2 = 4;
% Initialize IRS elements
w = exp(1i*rand(L,1)*(pi));
wrand = w;
% w = zeros(L,1);

allmaxsumrates = zeros(N,1);
numrealizations = 2;
num_iterations = 20;

% Initialize arrays to store Patrial and Pbtrial values
Patrial_evolution = zeros(num_iterations, numrealizations);
Pbtrial_evolution = zeros(num_iterations, numrealizations);
obj_evolution = zeros(num_iterations, numrealizations);

Paopt = zeros(1,numrealizations);
Pbopt = zeros(1,numrealizations);
optimum = zeros(1,numrealizations);

for r = 1:numrealizations
    Patrial = Pmin;
    Pbtrial = Pmin;
    
    [Paopt(r),Pbopt(r),optimum(r)] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,w);
    

    
    w_hat = [w;1];
    for zz = 1:num_iterations
        zz
        Patrial_evolution(zz,r) = Patrial;
        Pbtrial_evolution(zz,r) = Pbtrial;
    
        obj_evolution(zz,r) = sum_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Patrial,Pbtrial,w,sigma_ab,sigma_loop,sigma_c);
    
        y = sqrt(sigma_ab + sigma_loop + Pbtrial*abs(w_hat'*Hab(:,r))^2)/(Patrial*abs(w_hat'*Hac(:,r))^2+Pbtrial*abs(w_hat'*Hbc(:,r))^2+sigma_c);
        [Patrial, Pbtrial] = power_optimization_iter(y,Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,[w;1]);
    end
end

Patrial
Pbtrial

% Create the first figure for Paopt and Pbopt
figure;
yline(Paopt(1), 'r--', 'Paopt', 'LabelHorizontalAlignment', 'left');
hold on;
yline(Pbopt(1), 'g--', 'Pbopt', 'LabelHorizontalAlignment', 'left');
plot(1:num_iterations, Patrial_evolution(:,1), 'r-o');
plot(1:num_iterations, Pbtrial_evolution(:,1), 'g-o');
ylim([Pmin, Pmax]);
xlim([0, num_iterations+1]);
xlabel('X-axis');
ylabel('Power Values');
title('Optimization Results: Paopt and Pbopt');
legend('Paopt', 'Pbopt', 'Location', 'best');
grid on;
hold off;

% save('poweropt_sum_data','Paopt','Pbopt','num_iterations','Patrial_evolution','Pbtrial_evolution','Pmin','Pmax')

% Create the second figure for optimum
figure;
yline(optimum(1), 'b--', 'Optimum', 'LabelHorizontalAlignment', 'left');hold on;
plot(1:num_iterations, obj_evolution(:,1), 'b-o');
ylim([min(obj_evolution(:,1)) - 1, max(obj_evolution(:,1)) + 1]);
xlim([0, num_iterations+1]);
xlabel('X-axis');
ylabel('Power Value');
title('Optimization Result: Optimum');
legend('Optimum', 'Location', 'best');
grid on;
hold off;

fractional_increases = zeros(num_iterations-1,numrealizations);
diffwithexh = zeros(num_iterations,numrealizations);
for r = 1:numrealizations
    fractional_increases(:,r) = diff(obj_evolution(:,r))./obj_evolution(1:end-1,r);
    diffwithexh(:,r) = abs(obj_evolution(:,r) - optimum(r))/optimum(r);
end

figure;
semilogy(1:num_iterations-1, mean(fractional_increases,2))
hold on;
semilogy(1:num_iterations, mean(diffwithexh,2));
grid on;

% save('convergence_sum_power_2')

function [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab,Hac,Hbc,sigma_ab,sigma_loop,sigma_c,w)
    % Using bisection on grad to obtain the max point
    pa_low = Pmin;
    pa_high = Pmax;
    while abs(pa_low - pa_high) > 1e-10
        pa = (pa_high+pa_low)/2;
    
        gradvalue = calgrad(Hab,Hac,Hbc,w,sigma_ab,sigma_loop,sigma_c,Pmax,pa);
    
        if gradvalue < 0
            pa_high = pa;
        else
            pa_low = pa; 
        end
    end 
    Pasol = pa;
    Pbsol = Pmax - Pasol;
    
    optimalvalrate = sum_secrecy_rate(Hab,Hac,Hbc,Pasol,Pbsol,w,sigma_ab,sigma_loop,sigma_c);
    
    corner1 = sum_secrecy_rate(Hab,Hac,Hbc,Pmin,Pmin,w,sigma_ab,sigma_loop,sigma_c);
    corner2 = sum_secrecy_rate(Hab,Hac,Hbc,Pmin,Pmax-Pmin,w,sigma_ab,sigma_loop,sigma_c);
    corner3 = sum_secrecy_rate(Hab,Hac,Hbc,Pmax-Pmin,Pmin,w,sigma_ab,sigma_loop,sigma_c);
    
    % Optimizing between the corner points and the maximum point on Pa + Pb =
    % Pmax line
    [optimum,index]= max([corner1 corner2 corner3 optimalvalrate]);
    
    switch index
        case 1
            % disp(strcat('Pa = ',num2str(Pamin) ,' Pb = ',num2str(Pbmin) ))
            Paopt = Pmin;
            Pbopt = Pmin;
        case 2
            % disp(strcat('Pa = ',num2str(Pamin) ,' Pb = ',num2str(Pbmax) ))
            Paopt = Pmin;
            Pbopt = Pmax;
        case 3
            % disp(strcat('Pa = ',num2str(Pamax) ,' Pb = ',num2str(Pbmin) ))
            Paopt = Pmax;
            Pbopt = Pmin;
        otherwise
            % disp(strcat('Pa = ',num2str(Pasol) ,' Pb = ',num2str(Pbsol) ))
            Paopt = Pasol;
            Pbopt = Pbsol;
    end
end

function grad = calgrad(Hab,Hac,Hbc,w,sigma_ab,sigma_loop,sigma_c,Pmax,pa)

    wx = [w;1];

    x = real(trace(Hab*Hab'*wx*wx'));
    y = real(trace(Hac*Hac'*wx*wx'));
    z = real(trace(Hbc*Hbc'*wx*wx'));

    % Both gradients are correct
    % grad = x/(sigma_ab + sigma_loop+pa*x) - x/(sigma_ab + sigma_loop + (Pmax - pa)*x) - (y-z)/(sigma_c + Pmax*z + (y-z)*pa);
    grad = x*(sigma_ab + sigma_loop + (Pmax - pa)*x)*(sigma_c + Pmax*z + (y-z)*pa)...
        - x*(sigma_ab + sigma_loop+pa*x)*(sigma_c + Pmax*z + (y-z)*pa)...
        - (y-z)*(sigma_ab + sigma_loop+pa*x)*(sigma_ab + sigma_loop + (Pmax - pa)*x);

end

function sumrate = sum_secrecy_rate(Hab,Hac,Hbc,Pa,Pb,w,sigma_ab,sigma_loop,sigma_c)
    wx = [w;1];

    I_b_A = log2(1 + (Pb*abs(wx'*Hab)^2)/(sigma_ab+sigma_loop));
    I_a_B = log2(1 + (Pa*abs(wx'*Hab)^2)/(sigma_ab+sigma_loop));
    I_ab_C = log2(1 + (Pa*abs(wx'*Hac)^2 + Pb*abs(wx'*Hbc)^2)/(sigma_c));
    sumrate = max(I_b_A+I_a_B-I_ab_C,0);
end

function power = dbm2watt(dbm_value)
    % Converting power values from dBm to Watt for calculations

    power = (10^(dbm_value/10))*(10^(-3));
end