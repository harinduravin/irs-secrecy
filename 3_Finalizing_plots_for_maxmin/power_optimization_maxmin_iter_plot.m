function [Pa, Pb, Pa2, Pb2, Pa_evolution, Pb_evolution, Pa2_evolution, Pb2_evolution, obj_evolution] = power_optimization_maxmin_iter(Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab,Hac,Hbc,Ha2b2,Ha2c,Hb2c,Ha2a,Ha2b,Hb2a,Hb2b,sigma_ab,sigma_loop,sigma_c,w_hat,iter)
    Pa_evolution = zeros(iter, 1);
    Pb_evolution = zeros(iter, 1);
    Pa2_evolution = zeros(iter, 1);
    Pb2_evolution = zeros(iter, 1);
    obj_evolution = zeros(iter, 1);
    % fprintf('\n');
    for r = 1:iter
        Pa_evolution(r) = Pa;
        Pb_evolution(r) = Pb;
        Pa2_evolution(r) = Pa2;
        Pb2_evolution(r) = Pb2;
        obj_evolution(r) = min_secrecy_rate(Hab,Hac,Hbc,Ha2b2,Ha2c,Hb2c,Ha2a,Ha2b,Hb2a,Hb2b,Pa,Pb,Pa2,Pb2,w_hat,sigma_ab,sigma_loop,sigma_c);

        fprintf([num2str(r)]);
        y1 = sqrt(Pb*abs(w_hat'*Hab)^2)/(Pa2*abs(w_hat'*Ha2a)^2 + Pb2*abs(w_hat'*Hb2a)^2 + sigma_ab + sigma_loop);
        y2 = sqrt(Pa*abs(w_hat'*Hab)^2)/(Pa2*abs(w_hat'*Ha2b)^2 + Pb2*abs(w_hat'*Hb2b)^2 + sigma_ab + sigma_loop);
        y3 = sqrt(Pa2*abs(w_hat'*Ha2c)^2+Pb2*abs(w_hat'*Hb2c)^2+sigma_c)/(Pa*abs(w_hat'*Hac)^2 + Pb*abs(w_hat'*Hbc)^2 + Pa2*abs(w_hat'*Ha2c)^2+Pb2*abs(w_hat'*Hb2c)^2+sigma_c);
        y4 = sqrt(Pb2*abs(w_hat'*Ha2b2)^2)/(Pa*abs(w_hat'*Ha2a)^2 + Pb*abs(w_hat'*Ha2b)^2 + sigma_ab + sigma_loop);
        y5 = sqrt(Pa2*abs(w_hat'*Ha2b2)^2)/(Pa*abs(w_hat'*Hb2a)^2 + Pb*abs(w_hat'*Hb2b)^2 + sigma_ab + sigma_loop);
        y6 = sqrt(Pa*abs(w_hat'*Hac)^2+Pb*abs(w_hat'*Hbc)^2+sigma_c)/(Pa2*abs(w_hat'*Ha2c)^2 + Pb2*abs(w_hat'*Hb2c)^2 + Pa*abs(w_hat'*Hac)^2+Pb*abs(w_hat'*Hbc)^2+sigma_c);     
    
        cvx_begin quiet
            variable Pa
            variable Pb
            variable Pa2
            variable Pb2
            variable t
    
    
            maximize t
    
            subject to
    
            log(10000+2*y1*sqrt(100000000*Pb*abs(w_hat'*Hab)^2)...
            -(10000*Pa2*abs(w_hat'*Ha2a)^2 + 10000*Pb2*abs(w_hat'*Hb2a)^2 + 10000*sigma_ab + 10000*sigma_loop)*y1^2)...
            + log(10000+2*y2*sqrt(100000000*Pa*abs(w_hat'*Hab)^2)...
            -(10000*Pa2*abs(w_hat'*Ha2b)^2 + 10000*Pb2*abs(w_hat'*Hb2b)^2 + 10000*sigma_ab + 10000*sigma_loop)*y2^2)...
            + log(2*y3*sqrt(100000000*Pa2*abs(w_hat'*Ha2c)^2+100000000*Pb2*abs(w_hat'*Hb2c)^2+100000000*sigma_c)...
            -(10000*Pa*abs(w_hat'*Hac)^2 + 10000*Pb*abs(w_hat'*Hbc)^2 + 10000*Pa2*abs(w_hat'*Ha2c)^2+ 10000*Pb2*abs(w_hat'*Hb2c)^2+10000*sigma_c)*y3^2) >= t
    
            log(10000+2*y4*sqrt(100000000*Pb2*abs(w_hat'*Ha2b2)^2)...
            -(10000*Pa*abs(w_hat'*Ha2a)^2 + 10000*Pb*abs(w_hat'*Ha2b)^2 + 10000*sigma_ab + 10000*sigma_loop)*y4^2)...
            + log(10000+2*y5*sqrt(100000000*Pa2*abs(w_hat'*Ha2b2)^2)...
            -(10000*Pa*abs(w_hat'*Hb2a)^2 + 10000*Pb*abs(w_hat'*Hb2b)^2 + 10000*sigma_ab + 10000*sigma_loop)*y5^2)...
            + log(2*y6*sqrt(100000000*Pa*abs(w_hat'*Hac)^2+ 100000000*Pb*abs(w_hat'*Hbc)^2+sigma_c)...
            -(10000*Pa2*abs(w_hat'*Ha2c)^2 + 10000*Pb2*abs(w_hat'*Hb2c)^2 + 10000*Pa*abs(w_hat'*Hac)^2+ 10000*Pb*abs(w_hat'*Hbc)^2 + 10000*sigma_c)*y6^2) >= t
    
            Pa >= Pmin
            Pb >= Pmin 
            Pa2 >= Pmin
            Pb2 >= Pmin 
            Pb <= Pmax
            Pa <= Pmax
            Pb2 <= Pmax
            Pa2 <= Pmax
        cvx_end

    end

    % Rp1 = Rp1/log(2);
    % Rp2 = Rp2/log(2);    
end

function [minrate,p1rate,p2rate] = min_secrecy_rate(Hab,Hac,Hbc,Ha2b2,Ha2c,Hb2c,Ha2a,Ha2b,Hb2a,Hb2b,Pa,Pb,Pa2,Pb2,wx,sigma_ab,sigma_loop,sigma_c)
    % wx = [w;1];

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
