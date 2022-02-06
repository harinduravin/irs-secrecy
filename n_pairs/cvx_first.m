max_min = 0;
max_min_ind = 0;

for i = 0:15

    power_string = num2str(dec2bin(i,4));
    PA1 = get_power(power_string(1));
    PA2 = get_power(power_string(2));
    PB1 = get_power(power_string(3));
    PB2 = get_power(power_string(4));

    RA1 = get_P('A1') - get_Q('A1');
    RA2 = get_P('A2') - get_Q('A2');
    RB1 = get_P('B1') - get_Q('B1');
    RB2 = get_P('B2') - get_Q('B2');

    min_value = min([RA1,RA2,RB1,RB2]);

    if (min_value > max_min)

        max_min = min_value;
        max_min_ind = i;
    end

end

power_string = num2str(dec2bin(max_min_ind,4));
PA1 = get_power(power_string(1));
PA2 = get_power(power_string(2));
PB1 = get_power(power_string(3));
PB2 = get_power(power_string(4));

for j = 1:inner_iter

    cvx_begin quiet
    cvx_solver mosek
    variable X(L+1,L+1) complex semidefinite %symmetric
    variable t

    minimize t

    subject to

        -log(real(trace(get_leg_inf('A1')*X))+sigma_ab+sigma_loop) -log(real(trace(get_eve_inf('A1')*X)) +sigma_c) - real(trace(get_grad_S('A1')*(X-W))) <= t;

        -log(real(trace(get_leg_inf('B1')*X))+sigma_ab+sigma_loop) -log(real(trace(get_eve_inf('B1')*X)) +sigma_c) - real(trace(get_grad_S('B1')*(X-W))) <= t;

        -log(real(trace(get_leg_inf('A2')*X))+sigma_ab+sigma_loop) -log(real(trace(get_eve_inf('A2')*X)) +sigma_c) - real(trace(get_grad_S('A2')*(X-W))) <= t;

        -log(real(trace(get_leg_inf('B2')*X))+sigma_ab+sigma_loop) -log(real(trace(get_eve_inf('B2')*X)) +sigma_c) - real(trace(get_grad_S('B2')*(X-W))) <= t;

        diag(X) == 1;
        norm((X-W),1)<=2;

    cvx_end
end
% function grad = get_grad_S(node_name)

% end

function leg_inf = get_leg_inf(node_name)

    node_name(1)
    node_name(2)


end

% function eve_inf = get_eve_inf('ss')


% end

% permute(pagemtimes(permute(r,[1 3 2]),t),[1 3 2])