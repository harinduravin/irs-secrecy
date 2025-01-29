function opt_w = irs_optimization_maxmin_iter(L, Pmin, Pmax, Pa,Pb,Pa2,Pb2,Hab,Hac,Hbc,Ha2b2,Ha2c,Hb2c,Ha2a,Ha2b,Hb2a,Hb2b,sigma_ab,sigma_loop,sigma_c,w_hat,iter)

    wx = w_hat;
    Wt = wx*wx';
    % fprintf('\n');

    for r = 1:iter
        fprintf([num2str(r)]);
        
        cvx_begin quiet
            variable W(L+1,L+1) complex semidefinite
            variable t

            maximize t
            % User pair 1
            Q11 = log(10000*sigma_ab+10000*sigma_loop+Pb2*real(trace(10000*Hb2a*Hb2a'*W))+Pa2*real(trace(10000*Ha2a*Ha2a'*W)) + Pb*real(trace(10000*Hab*Hab'*W)))-log(10000);
            Q12 = log(10000*sigma_ab+10000*sigma_loop+Pb2*real(trace(10000*Hb2b*Hb2b'*W))+Pa2*real(trace(10000*Ha2b*Ha2b'*W)) + Pa*real(trace(10000*Hab*Hab'*W)))-log(10000);
            S13 = log(sigma_c+Pb2*real(trace(Hb2c*Hb2c'*Wt))+Pa2*real(trace(Ha2c*Ha2c'*Wt)) + Pa*real(trace(Hac*Hac'*Wt)) + Pb*real(trace(Hbc*Hbc'*Wt)));
            S11 = log(sigma_ab+sigma_loop+Pb2*real(trace(Hb2a*Hb2a'*Wt))+Pa2*real(trace(Ha2a*Ha2a'*Wt)));
            S12 = log(sigma_ab+sigma_loop+Pb2*real(trace(Hb2b*Hb2b'*Wt))+Pa2*real(trace(Ha2b*Ha2b'*Wt)));
            Q13 = log(10000*sigma_c+Pb2*real(trace(10000*Hb2c*Hb2c'*W))+Pa2*real(trace(10000*Ha2c*Ha2c'*W)))-log(10000);

            deltaS13 = (Pb2*Hb2c*Hb2c' + Pa2*Ha2c*Ha2c' + Pa*Hac*Hac' + Pb*Hbc*Hbc')/(sigma_c+Pb2*real(trace(Hb2c*Hb2c'*Wt))+Pa2*real(trace(Ha2c*Ha2c'*Wt)) + Pa*real(trace(Hac*Hac'*Wt)) + Pb*real(trace(Hbc*Hbc'*Wt)));
            deltaS11 = (Pb2*Hb2a*Hb2a' + Pa2*Ha2a*Ha2a')/(sigma_ab+sigma_loop+Pb2*real(trace(Hb2a*Hb2a'*Wt))+Pa2*real(trace(Ha2a*Ha2a'*Wt)));
            deltaS12 = (Pb2*Hb2b*Hb2b' + Pa2*Ha2b*Ha2b')/(sigma_ab+sigma_loop+Pb2*real(trace(Hb2b*Hb2b'*Wt))+Pa2*real(trace(Ha2b*Ha2b'*Wt)));
            deltaS1 = (deltaS11+deltaS12+deltaS13)/log(2);

            % User pair 2
            Q21 = log(10000*sigma_ab+10000*sigma_loop+Pb*real(trace(10000*Ha2b*Ha2b'*W))+Pa*real(trace(10000*Ha2a*Ha2a'*W)) + Pb2*real(trace(10000*Ha2b2*Ha2b2'*W)))-log(10000);
            Q22 = log(10000*sigma_ab+10000*sigma_loop+Pb*real(trace(10000*Hb2b*Hb2b'*W))+Pa*real(trace(10000*Hb2a*Hb2a'*W)) + Pa2*real(trace(10000*Ha2b2*Ha2b2'*W)))-log(10000);
            S23 = log(sigma_c+Pb*real(trace(Hbc*Hbc'*Wt))+Pa*real(trace(Hac*Hac'*Wt)) + Pa2*real(trace(Ha2c*Ha2c'*Wt)) + Pb2*real(trace(Hb2c*Hb2c'*Wt)));
            S21 = log(sigma_ab+sigma_loop+Pb*real(trace(Ha2b*Ha2b'*Wt))+Pa*real(trace(Ha2a*Ha2a'*Wt)));
            S22 = log(sigma_ab+sigma_loop+Pb*real(trace(Hb2b*Hb2b'*Wt))+Pa*real(trace(Hb2a*Hb2a'*Wt)));
            Q23 = log(10000*sigma_c+Pb*real(trace(10000*Hbc*Hbc'*W))+Pa*real(trace(10000*Hac*Hac'*W)))-log(10000);

            deltaS23 = (Pb*Hbc*Hbc' + Pa*Hac*Hac' + Pa2*Ha2c*Ha2c' + Pb2*Hb2c*Hb2c')/(sigma_c+Pb*real(trace(Hbc*Hbc'*Wt))+Pa*real(trace(Hac*Hac'*Wt)) + Pa2*real(trace(Ha2c*Ha2c'*Wt)) + Pb2*real(trace(Hb2c*Hb2c'*Wt)));
            deltaS21 = (Pb*Ha2b*Ha2b' + Pa*Ha2a*Ha2a')/(sigma_ab+sigma_loop+Pb*real(trace(Ha2b*Ha2b'*Wt))+Pa*real(trace(Ha2a*Ha2a'*Wt)));
            deltaS22 = (Pb*Hb2b*Hb2b' + Pa*Hb2a*Hb2a')/(sigma_ab+sigma_loop+Pb*real(trace(Hb2b*Hb2b'*Wt))+Pa*real(trace(Hb2a*Hb2a'*Wt)));
            deltaS2 = (deltaS21+deltaS22+deltaS23)/log(2);

            subject to
            Q11+Q12+Q13-S11-S12-S13 - real(trace(deltaS1*(W-Wt))) >= t
            Q21+Q22+Q23-S21-S22-S23 - real(trace(deltaS2*(W-Wt))) >= t

            diag(W) == 1;
            % norm(W-Wt) <= 0.2
        cvx_end
        Wt = W;
    end 
    [U,~,~] = svd(W);
    w_extracted = U(:,1);
    opt_w = exp(1i*angle(w_extracted/w_extracted(end)));
end

    % I_b_A = log2(1 + (Pb*real(trace(Hab*Hab'*wx*wx')))/(sigma_ab+sigma_loop+Pb2*real(trace(Hb2a*Hb2a'*wx*wx'))+Pa2*real(trace(Ha2a*Ha2a'*wx*wx'))));
    % I_a_B = log2(1 + (Pa*real(trace(Hab*Hab'*wx*wx')))/(sigma_ab+sigma_loop+Pb2*real(trace(Hb2b*Hb2b'*wx*wx'))+Pa2*real(trace(Ha2b*Ha2b'*wx*wx'))));
    % I_ab_C = log2(1 + (Pa*real(trace(Hac*Hac'*wx*wx')) + Pb*real(trace(Hbc*Hbc'*wx*wx')))/(sigma_c+Pb2*real(trace(Hb2c*Hb2c'*wx*wx'))+Pa2*real(trace(Ha2c*Ha2c'*wx*wx'))));

    % Q11 = log(sigma_ab+sigma_loop+Pb2*real(trace(Hb2a*Hb2a'*wx*wx'))+Pa2*real(trace(Ha2a*Ha2a'*wx*wx')) + Pb*real(trace(Hab*Hab'*wx*wx')));
    % Q12 = log(sigma_ab+sigma_loop+Pb2*real(trace(Hb2b*Hb2b'*wx*wx'))+Pa2*real(trace(Ha2b*Ha2b'*wx*wx')) + Pa*real(trace(Hab*Hab'*wx*wx')));
    % S13 = log(sigma_c+Pb2*real(trace(Hb2c*Hb2c'*wx*wx'))+Pa2*real(trace(Ha2c*Ha2c'*wx*wx')) + Pa*real(trace(Hac*Hac'*wx*wx')) + Pb*real(trace(Hbc*Hbc'*wx*wx')));
    % S11 = log(sigma_ab+sigma_loop+Pb2*real(trace(Hb2a*Hb2a'*wx*wx'))+Pa2*real(trace(Ha2a*Ha2a'*wx*wx')));
    % S12 = log(sigma_ab+sigma_loop+Pb2*real(trace(Hb2b*Hb2b'*wx*wx'))+Pa2*real(trace(Ha2b*Ha2b'*wx*wx')));
    % Q13 = log(sigma_c+Pb2*real(trace(Hb2c*Hb2c'*wx*wx'))+Pa2*real(trace(Ha2c*Ha2c'*wx*wx')));
    % p1rate = max(Q11+Q12+Q13-S11-S12-S13,0)/log(2);
    
    % I2_b_A = log2(1 + (Pb2*real(trace(Ha2b2*Ha2b2'*wx*wx')))/(sigma_ab+sigma_loop+Pb*real(trace(Ha2b*Ha2b'*wx*wx'))+Pa*real(trace(Ha2a*Ha2a'*wx*wx'))));
    % I2_a_B = log2(1 + (Pa2*real(trace(Ha2b2*Ha2b2'*wx*wx')))/(sigma_ab+sigma_loop+Pb*real(trace(Hb2b*Hb2b'*wx*wx'))+Pa*real(trace(Hb2a*Hb2a'*wx*wx'))));
    % I2_ab_C = log2(1 + (Pa2*real(trace(Ha2c*Ha2c'*wx*wx')) + Pb2*real(trace(Hb2c*Hb2c'*wx*wx')))/(sigma_c+Pb*real(trace(Hbc*Hbc'*wx*wx'))+Pa*real(trace(Hac*Hac'*wx*wx'))));

    % Q21 = log(sigma_ab+sigma_loop+Pb*real(trace(Ha2b*Ha2b'*wx*wx'))+Pa*real(trace(Ha2a*Ha2a'*wx*wx')) + Pb2*real(trace(Ha2b2*Ha2b2'*wx*wx')));
    % Q22 = log(sigma_ab+sigma_loop+Pb*real(trace(Hb2b*Hb2b'*wx*wx'))+Pa*real(trace(Hb2a*Hb2a'*wx*wx')) + Pa2*real(trace(Ha2b2*Ha2b2'*wx*wx')));
    % S23 = log(sigma_c+Pb*real(trace(Hbc*Hbc'*wx*wx'))+Pa*real(trace(Hac*Hac'*wx*wx')) + Pa2*real(trace(Ha2c*Ha2c'*wx*wx')) + Pb2*real(trace(Hb2c*Hb2c'*wx*wx')));
    % S21 = log(sigma_ab+sigma_loop+Pb*real(trace(Ha2b*Ha2b'*wx*wx'))+Pa*real(trace(Ha2a*Ha2a'*wx*wx')));
    % S22 = log(sigma_ab+sigma_loop+Pb*real(trace(Hb2b*Hb2b'*wx*wx'))+Pa*real(trace(Hb2a*Hb2a'*wx*wx')));
    % Q23 = log(sigma_c+Pb*real(trace(Hbc*Hbc'*wx*wx'))+Pa*real(trace(Hac*Hac'*wx*wx')));
    % p2rate = max(Q21+Q22+Q23-S21-S22-S23,0)/log(2);
    % 
    % minrate = min(p1rate, p2rate)

        % deltaS13 = (Pb2*Hb2c + Pa2*Ha2c + Pa*Hac+Pb*Hbc)/(sigma_c+Pb2*real(trace(Hb2c*Hb2c'*Wt))+Pa2*real(trace(Ha2c*Ha2c'*Wt)) + Pa*real(trace(Hac*Hac'*Wt)) + Pb*real(trace(Hbc*Hbc'*Wt)));
    % deltaS11 = (Pb2*Hb2a + Pa2*Ha2a)/(sigma_ab+sigma_loop+Pb2*real(trace(Hb2a*Hb2a'*Wt))+Pa2*real(trace(Ha2a*Ha2a'*Wt)));
    % deltaS12 = (Pb2*Hb2b + Pa2*Ha2b)/(sigma_ab+sigma_loop+Pb2*real(trace(Hb2b*Hb2b'*Wt))+Pa2*real(trace(Ha2b*Ha2b'*Wt)));
    % deltaS1 = (deltaS11+deltaS12+deltaS13)/log(2);