clear;
close all;
load('L4channels.mat','Hab','Hac','Hbc','L','N','sigma_ab','sigma_c','sigma_loop');
Pmax = dbm2watt(15);
Pamax = dbm2watt(15);
Pbmax = dbm2watt(15);

Pmin  = dbm2watt(0);
Pamin = dbm2watt(0);
Pbmin = dbm2watt(0);

ticks = 40;
Pavalues = linspace(Pamin,Pamax,ticks);

Pbvalues = linspace(Pbmin,Pbmax,ticks);

% Initialize IRS elements
w = exp(1i*rand(L,1)*(pi));
% w = zeros(L,1);

% Using bisection on grad to obtain the max point
pa_low = Pmin;
pa_high = Pmax;
while abs(pa_low - pa_high) > 1e-10
    pa = (pa_high+pa_low)/2;

    gradvalue = calgrad(Hab(:,2),Hac(:,2),Hbc(:,2),w,sigma_ab,sigma_loop,sigma_c,Pmax,pa);

    if gradvalue < 0
        pa_high = pa;
    else
        pa_low = pa; 
    end
end 
Pasol = pa;
Pbsol = Pmax - Pasol;

optimalval = get_G(Hab(:,2),Hac(:,2),Hbc(:,2),Pasol,Pbsol,w,sigma_ab,sigma_loop,sigma_c,Pmax);
optimalvalrate = sum_secrecy_rate(Hab(:,2),Hac(:,2),Hbc(:,2),Pasol,Pbsol,w,sigma_ab,sigma_loop,sigma_c);

crosssec = zeros(1,ticks);
crosssecrate = zeros(1,ticks);
gradplot = zeros(1,ticks);

for q = 1:ticks
    crosssec(q) = get_G(Hab(:,2),Hac(:,2),Hbc(:,2),Pavalues(q),Pmax -Pavalues(q),w,sigma_ab,sigma_loop,sigma_c,Pmax);
    crosssecrate(q) = sum_secrecy_rate(Hab(:,2),Hac(:,2),Hbc(:,2),Pavalues(q),Pmax -Pavalues(q),w,sigma_ab,sigma_loop,sigma_c);
    gradplot(q) = calgrad(Hab(:,2),Hac(:,2),Hbc(:,2),w,sigma_ab,sigma_loop,sigma_c,Pmax,Pavalues(q));
end

figure
plot(Pavalues,crosssec); hold on;
scatter(Pasol,optimalval)
title('Value of G function along Pmax = Pa + Pb line')
xlabel('Value of Pa')
ylabel('Function G')
grid on;

figure;
plot(Pavalues,crosssecrate); hold on;
scatter(Pasol,optimalvalrate);
title('Sum Secrecy rate along Pmax = Pa + Pb line')
xlabel('Value of Pa')
ylabel('Sum Secrecy rate')
grid on;

figure;
plot(Pavalues(1:end-3),gradplot(1:end-3)); hold on;
yline(0)
title('Gradient of G function along Pmax = Pa + Pb line')
xlabel('Value of Pa')
ylabel('Sum Secrecy rate')
grid on;

gfunc = zeros(ticks);
for a = 1:ticks
    for b = 1:ticks
        gfunc(a,b) = sum_secrecy_rate(Hab(:,2),Hac(:,2),Hbc(:,2),Pavalues(a),Pbvalues(b),w,sigma_ab,sigma_loop,sigma_c);
    end
end

figure;
surf(Pavalues,Pbvalues,real(gfunc.')); 
hold on;
plot3(Pasol,Pbsol,optimalvalrate,'.r','MarkerSize', 30)
title('Sum secrecy rate for possible Pa and Pb values')
xlabel('Value of Pa')
ylabel('Value of Pb')
zlabel('Sum Secrecy rate')


corner1 = sum_secrecy_rate(Hab(:,2),Hac(:,2),Hbc(:,2),Pamin,Pbmin,w,sigma_ab,sigma_loop,sigma_c);
corner2 = sum_secrecy_rate(Hab(:,2),Hac(:,2),Hbc(:,2),Pamin,Pbmax-Pbmin,w,sigma_ab,sigma_loop,sigma_c);
corner3 = sum_secrecy_rate(Hab(:,2),Hac(:,2),Hbc(:,2),Pamax-Pamin,Pbmin,w,sigma_ab,sigma_loop,sigma_c);

% Optimizing between the corner points and the maximum point on Pa + Pb =
% Pmax line
[optimum,index]= max([corner1 corner2 corner3 optimalvalrate]);

switch index
    case 1
        disp(strcat('Pa = ',num2str(Pamin) ,' Pb = ',num2str(Pbmin) ))
    case 2
        disp(strcat('Pa = ',num2str(Pamin) ,' Pb = ',num2str(Pbmax) ))
    case 3
        disp(strcat('Pa = ',num2str(Pamax) ,' Pb = ',num2str(Pbmin) ))
    otherwise
        disp(strcat('Pa = ',num2str(Pasol) ,' Pb = ',num2str(Pbsol) ))
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

function gvalue = get_G(Hab,Hac,Hbc,Pa,Pb,w,sigma_ab,sigma_loop,sigma_c,Pmax)

    wx = [w;1];
    x = real(trace(Hab*Hab'*wx*wx'));
    y = real(trace(Hac*Hac'*wx*wx'));
    z = real(trace(Hbc*Hbc'*wx*wx'));

    % gvalue = log2(sigma_ab + sigma_loop +Pa*trace(Hab*Hab'*wx*wx')) + log2(sigma_ab + sigma_loop + Pb*trace(Hab*Hab'*wx*wx')) - log2(sigma_c + Pb*trace(Hbc*Hbc'*wx*wx') + Pa*trace(Hac*Hac'*wx*wx'));
    gvalue = log2(sigma_ab + sigma_loop +Pa*x) + log2(sigma_ab + sigma_loop + (Pmax-Pa)*x)- log2(sigma_c + (Pmax-Pa)*z + Pa*y);

end

function maxval = maxline(Hab,Hac,Hbc,w,sigma_ab,sigma_loop,sigma_c,Pmax)

    wx = [w;1];

    x = real(trace(Hab*Hab'*wx*wx'));
    y = real(trace(Hac*Hac'*wx*wx'));
    z = real(trace(Hbc*Hbc'*wx*wx'));

    syms pa

    eqn = (x/(sigma_ab + sigma_loop+pa*x) - x/(sigma_ab + sigma_loop + (Pmax - pa)*x) - (y-z)/(sigma_c + Pmax*z + (y-z)*pa)) == 0;
    maxval = solve(eqn, pa);

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