function [Pa Pb] = power_optimization_iter(y,Pmin,Pmax,Hab,Hac,Hbc,sigma_ab,sigma_loop,sigma_c,w_hat)
    % optimal beamformer SDP solution
    % y = 1;
    % Hab = Hab * 1000;
    % Hbc = Hbc * 1000;
    % Hac = Hac * 1000;
    cvx_begin quiet
        variable Pa
        variable Pb
        maximize (log(10000*sigma_ab + 10000*sigma_loop + 10000*Pa*abs(w_hat'*Hab)^2)...
        + log(2*y*sqrt(100000000*sigma_ab + 100000000*sigma_loop + 100000000*Pb*abs(w_hat'*Hab)^2)...
        -(10000*Pa*abs(w_hat'*Hac)^2+10000*Pb*abs(w_hat'*Hbc)^2+10000*sigma_c)*y^2))
        subject to
        Pa >= Pmin
        Pb >= Pmin 
        Pa + Pb <= Pmax
    cvx_end

    % cvx_status

end
