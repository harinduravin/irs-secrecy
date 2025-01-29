clear;
close all;
load('L4channels.mat','Hab','Hac','Hbc','L','N','sigma_ab','sigma_c','sigma_loop');

Pmax = dbm2watt(15);
Pmin  = dbm2watt(0);

F = 10000;
I1 = 5;
I2 = 5;
% Initialize IRS elements
w = exp(1i*rand(L,1)*(pi));
wrand = w;
% w = zeros(L,1);

allmaxsumrates = zeros(N,1);

% for r = 1:N
%     r
%     maxsumrate = 0;
%     for a = 1:F
%         % w(1) = exp(1i*(a/F)*(2*pi));
%         w(1) = 0;        
%         [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,w);
%         maxsumrate = max(maxsumrate,optimum);
%     end
%     allmaxsumrates(r) = maxsumrate;
% end

% For C1 values in script and paper became very close when zeros were
% ignored in averaging. In C2 no-irs case resulted in 0.64 where 0.85 is
% expected. Anyway conclusion is to use Lo twice for IRS reflection. Not to
% be stuck with channel generation (Possibilities of error in original paper).
% Check whether the SDP method results in
% exactly similar value as exhaustion.
r = 2;
for i1 = 1:I1
    [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,w);
    for i2 = 1:I2
        opt_w = irs_optimization(L,Paopt,Pbopt,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,[w;1]);
        w = opt_w(1:end-1);
    end
end
optimumvaluesdp = sum_secrecy_rate(Hab(:,r),Hac(:,r),Hbc(:,r),Paopt,Pbopt,opt_w(1:end-1),sigma_ab,sigma_loop,sigma_c);

% w_hat = [w;1];
% y = sqrt(sigma_ab + sigma_loop + Pb*abs(w_hat'*Hab)^2)/(Pa*abs(w_hat'*Hac)^2+Pb*abs(w_hat'*Hbc)^2+sigma_c);
% [Patrial, Pbtrial] = power_optimization_iter(y,Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,[w;1]);

optimumvalueexh = 0;
for a = 1:F
    w(1) = exp(1i*(a/F)*(2*pi));
    % w(1) = 0;        
    [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,w);
    optimumvalueexh = max(optimumvalueexh,optimum);
end
optimumvalueexh
optimumvaluesdp
[~,~,optimumvaluerand] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,wrand)
% 
% for r = 1:N
%      [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab(:,r),Hac(:,r),Hbc(:,r),sigma_ab,sigma_loop,sigma_c,w);
%      allmaxsumrates(r) = optimum;
% end
% 
% mean(nonzeros(allmaxsumrates))
% histogram(allmaxsumrates)

% [Paopt,Pbopt,optimum] = optimizepower(Pmin,Pmax,Hab(:,2),Hac(:,2),Hbc(:,2),sigma_ab,sigma_loop,sigma_c,w);

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
    corner2 = sum_secrecy_rate(Hab,Hac,Hbc,Pmin,Pmax,w,sigma_ab,sigma_loop,sigma_c);
    corner3 = sum_secrecy_rate(Hab,Hac,Hbc,Pmax,Pmin,w,sigma_ab,sigma_loop,sigma_c);
    
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