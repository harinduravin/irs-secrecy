function opt_w = irs_optimization(L,Pa,Pb,Hab,Hac,Hbc,sigma_ab,sigma_loop,sigma_c,w_hat)
    % optimal beamformer SDP solution

    W_hat = w_hat*w_hat';
    Hab = Hab*100;
    Hbc = Hbc;
    Hac = Hac;

    cvx_begin quiet
        variable W(L+1,L+1) complex semidefinite
        maximize (-2*log(0.0001)+(log(10000*sigma_ab + 10000*sigma_loop + real(Pa*trace(Hab*Hab'*W)))...
        + log(10000*sigma_ab + 10000*sigma_loop +Pb*real(trace(Hab*Hab'*W))))/log(2)...
        - (Pa*real(trace(Hac*Hac'*(W-W_hat))) + Pb*real(trace(Hbc*Hbc'*(W-W_hat))))...
        /(sigma_c + Pa*real(trace(Hac*Hac'*W_hat)) + Pb*real(trace(Hbc*Hbc'*W_hat))))
        subject to
        diag(W) == 1;
        norm(W-W_hat) <= 1
    cvx_end

    cvx_status

        % - log(sigma_c + 0.0001*Pa*real(trace(Hac*Hac'*W_hat)) + 0.0001*Pb*real(trace(Hbc*Hbc'*W_hat)))

    [U,~,~] = svd(W);
    w_extracted = U(:,1);
    opt_w = exp(1i*angle(w_extracted/w_extracted(end)));
end

% function  removed = remove_small(X)
%     indices = abs(X)<=(10^(-7));
%     X(indices) = 0 + 0*1i;
%     removed = X;
% end