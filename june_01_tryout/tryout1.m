% Initializing number of pairs
global n;
n = 1;

% Position of all the nodes
coords = [-15.83 12.03;-15.96 0.03;-0.25 11.97;
0.08 -0.27;];

% Preparing the Euclidean distance matrix [A1, B1, A2, B2, ...., IRS, C]
global dist;
dist = sqDistance(coords, coords);

L_list = 5:5:25;

% Residual interference and noise values
sigma_c = 10^(-14.5);
sigma_ab = 10^(-14.5);
sigma_loop = 10^(-14);

% Maximum and minimum transmit powers of users in dB
max_power = 15;
min_power = 0;

% Preparing all the users in a list
global user_list;
user_list = [];

for j = 1:n

    string1 = string(strcat('A',num2str(j)));
    user_list = [user_list string1];
    string2 = string(strcat('B',num2str(j)));
    user_list = [user_list string2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_master_list = [];

num_seeds = 100 - 1;
f = waitbar(0,'Please wait...');


global seed;
for seed = 25:num_seeds
    seed_flag = 0;

    waitbar(seed/num_seeds,f,'Simulating...');
    min_list = [];

    % Number of elements in the IRS
    global L;
    for L = 5:5:25

        % IRS reflection matrix (Initialized randomly)
        wi = exp(1i*(2*rand(L,1)-1)*pi);
        w = [wi ; 1].';

        global W; % Check this again
        W = w'*w;

        % Generating all the channels before algorithm loop starts
        %
        global leg_inf_with_opp_stacks;
        global eves_inf_with_self_stacks;
        global leg_inf_stacks;
        global eves_stack;


        % Revise this part
        leg_inf_stacks = get_stacks("leg_inf");
        leg_inf_with_opp_stacks = get_stacks("leg_inf_with_opp");
        eves_inf_with_self_stacks = get_stacks("eves_inf_with_self");
        eves_stack = eves_inf_all();

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        outer_iter = 20;

        for j = 1:outer_iter
            %
            % Power optimization - Successive convex approximation
            %
            power_string = num2str(dec2bin(0,2*n));
            P_hat = (get_power_values(power_string,max_power,min_power))';

            % Remove comment when removing power optimization (Fixed average power)
            % power_string_high = num2str(dec2bin(2^(2*n)-1,2*n));
            % P_hat_high = (get_power_values(power_string_high,max_power,min_power))';
            % P_hat = P_hat_high;

            %%%%%%%%%%%%%%%%%%%%%%%%%


            cvx_begin quiet
            cvx_solver Mosek
            variable P(2*n,1)
            variable z

            minimize z

            subject to

                for u = 1:length(user_list)


                    % for n = 1
                    -log(sigma_ab + sigma_loop + P(p_(user_list(u),1,1))*H_(user_list(u),1,1)+ P(p_(user_list(u),1,2))*H_(user_list(u),1,2) + P(p_(user_list(u),1,3))*H_(user_list(u),1,3))...
                    -log(sigma_c + P(p_(user_list(u),2,1))*H_(user_list(u),2,1)+ P(p_(user_list(u),2,2))*H_(user_list(u),2,2) + P(p_(user_list(u),2,3))*H_(user_list(u),2,3))...
                    -get_S(user_list(u),P_hat)...
                    -(get_grad_P(user_list(u),P_hat,W))'*(P-P_hat) <= z;

                    % for n = 2
                    % -log(sigma_ab + sigma_loop + P(p_(user_list(u),1,1))*H_(user_list(u),1,1)+ P(p_(user_list(u),1,2))*H_(user_list(u),1,2) + P(p_(user_list(u),1,3))*H_(user_list(u),1,3))...
                    % -log(sigma_c + P(p_(user_list(u),2,1))*H_(user_list(u),2,1)+ P(p_(user_list(u),2,2))*H_(user_list(u),2,2) + P(p_(user_list(u),2,3))*H_(user_list(u),2,3))...
                    % -get_S(user_list(u),P_hat)...
                    % -(get_grad_P(user_list(u),P_hat,W))'*(P-P_hat) <= z;

                    % for n = 3
                    % -log(sigma_ab + sigma_loop + P(p_(user_list(u),1,1))*H_(user_list(u),1,1)+ P(p_(user_list(u),1,2))*H_(user_list(u),1,2) + P(p_(user_list(u),1,3))*H_(user_list(u),1,3) + P(p_(user_list(u),1,4))*H_(user_list(u),1,4) + P(p_(user_list(u),1,5))*H_(user_list(u),1,5))...
                    % -log(sigma_c + P(p_(user_list(u),2,1))*H_(user_list(u),2,1)+ P(p_(user_list(u),2,2))*H_(user_list(u),2,2) + P(p_(user_list(u),2,3))*H_(user_list(u),2,3) + P(p_(user_list(u),2,4))*H_(user_list(u),2,4)+ P(p_(user_list(u),2,5))*H_(user_list(u),2,5))...
                    % -get_S(user_list(u),P_hat)...
                    % - (get_grad_P(user_list(u),P_hat,W))'*(P-P_hat) <= z;

                end

                for i = 1:length(user_list)

                    P(i) <= dbm2watt(max_power); 
                    P(i) >= dbm2watt(min_power);

                end
                norm((P-P_hat),1)<=1;
            cvx_end

            P_hat = P;

            %%%%%%%%%%%%%%%%%%%%%%

            % Reflection phase optimization
            % Successive Convex Approximation (SCA)
            %

            try
                cvx_begin quiet
                % cvx_precision high
                cvx_solver Mosek
                variable X(L+1,L+1) complex semidefinite %symmetric
                variable t

                minimize t

                subject to

                    for u = 1:length(user_list)

                        % -log(real(trace(get_cvx_leg_inf(user_list(u),P_hat)*X))+sigma_ab+sigma_loop) -log(real(trace(get_cvx_eve_inf(user_list(u),P_hat)*X)) +sigma_c) - real(trace(get_grad_S(user_list(u),P_hat,W)*(X-W))) <= t;
                        -log(real(trace(remove_small(get_cvx_leg_inf(user_list(u),P_hat))*X))+sigma_ab+sigma_loop) -log(real(trace(remove_small(get_cvx_eve_inf(user_list(u),P_hat))*X)) +sigma_c) - real(trace(remove_small(get_grad_S(user_list(u),P_hat,W))*(X-remove_small(W)))) <= t;
                    end

                    diag(X) == 1;
                    % norm((X-W),1)<=2*exp(-outer_iter/2);
                    norm((X-W),1)<=2;
                cvx_end
            catch

                fprintf("Error occured")
                cvx_status
                seed
                seed_flag = 1;
                break;
            end

            if cvx_status == 'Failed'
                cvx_status
                seed
                L
                continue;
            else
                W = X;
            end

            if seed_flag == 1
                break;
            end
            [min_value, rate_list] = get_min_rate(P_hat);
        end % end of outer iteration
        if seed_flag == 1
            break;
        end

        [min_value, rate_list] = get_min_rate(P_hat);
        min_list = [min_list, min_value];

    end % end of L iteration
    min_master_list = cat(1,min_master_list,min_list);
end % end of random seed iteration

close(f)

figure(1)
plot(L_list,mean(min_master_list));
ylim([2 4])
xlabel('Number of elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function D = sqDistance(X, Y)
    % Obtaining the Euclidean distance matrix from given x,y coordinates.
    %

    D = sqrt(bsxfun(@plus,dot(X,X,2)',dot(Y,Y,2))-2*(X*Y'));
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leg_inf = get_leg_inf(node_index,indices_list)
    % Stacks channels from given set of users to a given node.
    % These 2D channels are stacked along the 3rd dimension for future use.

    stack = [];

    for p = 1:length(indices_list)
        stack = cat(3, stack,get_leg_channel(node_index,indices_list(p)));
    end

    leg_inf = stack;

end

function eves_inf = get_eves_inf(indices_list)
    % Stacks channels from given set of users to the eavesdropper.
    % These 2D channels are stacked along the 3rd dimension for future use.

    stack = [];

    for p = 1:length(indices_list)
        stack = cat(3, stack,get_eves_channel(indices_list(p)));
    end

    eves_inf = stack;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leg_channel = get_leg_channel(rx_index,tx_index)
    % Provides the channel between given transmitter and receiver node
    % This depends on,
    % 1) Pathloss factor 
    % 2) Number of elements in the IRS
    % 3) Distances between RX, IRS, TX

    global dist;
    global n;
    global L;
    global seed;

    Lo = 10^(-3);
    pl = 2;

    irs_index = 2*n+1;

    % Obtaining distances from Euclidean distance matrix
    dti = dist(tx_index, irs_index);
    dir = dist(irs_index, rx_index);
    dtr = dist(rx_index,tx_index);

    % Calculation of the random seed from rx_index and tx_index
    total_channels = 2*n^2 + n;
    % total_leg_channels = 2*n^2 - n;

    start_val_leg = seed*total_channels;
    u = zeros(2*n,2*n);
    p = start_val_leg;
    for i = 1:2*n
        for j = 1:2*n
            if j>i
                p = p+1;
                u(j,i) = p;
            end

        end
    end
    nmatrix = u + u';

    input_seed = nmatrix(rx_index,tx_index);

    H = get_H(dti,dir,dtr,L,Lo,pl,input_seed); % TODO Change this line to complete the randomness
    leg_channel = H*H';
end

function eves_channel = get_eves_channel(tx_index)
    % Provides the channel between given transmitter and eavesdropper node
    % This depends on,
    % 1) Pathloss factor 
    % 2) Number of elements in the IRS
    % 3) Distances between Eavesdropper, IRS, TX

    global dist;
    global n;
    global L;
    global seed;

    Lo = 10^(-3);
    pl = 2;

    irs_index = 2*n+1;
    rx_index = 2*n+2;

    % Obtaining distances from Euclidean distance matrix
    dti = dist(tx_index, irs_index);
    dir = dist(irs_index, rx_index);
    dtr = dist(rx_index,tx_index);

    % Calculation of the random seed from rx_index and tx_index
    total_channels = 2*n^2 + n;
    total_leg_channels = 2*n^2 - n;
    start_val_eve = seed*total_channels+total_leg_channels;

    input_seed = tx_index + start_val_eve;

    H = get_H(dti,dir,dtr,L,Lo,pl,input_seed); % TODO Change this line to complete the randomness

    eves_channel = H*H';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices_list = convert2indices(node_list)
    % [A1, B1, A2, B2, ...]  ---- This is converted to 1, 2, 3, 4,... indices.
    % This helps accessing channels from compact arrays.
    %

    initial_list = [];
    for p = 1:length(node_list)

        node = char(node_list(p));

        node_type = node(1);
        pair_num = str2num(node(2));

        if node_type == 'A'
            index = 2*pair_num -1;
        else 
            index = 2*pair_num;
        end

        initial_list = [initial_list index];
    end
    indices_list = initial_list;

end

function index_val = get_index(user)
    % Similar to convert2indices function. Instead of a list of user names, a single 
    % name is converted to its index.
    %

    node = char(user);
    node_type = node(1);
    pair_num = str2num(node(2));

    if node_type == 'A'
        index = 2*pair_num -1;
    else 
        index = 2*pair_num;
    end

    index_val = index;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = get_H(dti,dir,dtr,L,Lo,pl,rng_val)
    % Generates H. This matrix stores both direct channel and the channel through
    % the IRS. Two channels are concatenated to form a single column matrix.
    %
    
    Lo2 = 0.001;
    Lo3 = 0.001;

    gti = get_ricean_channel(dti, 2, L, 5*rng_val);
    hti = sqrt(Lo2)*gti;
    gir = get_ricean_channel(dir, 2, L, 5*rng_val+2);
    hir = sqrt(Lo2)*gir;
    gtr = get_ricean_direct_channel(dtr, 2.4, 5*rng_val+4);
    htr = sqrt(Lo3)*gtr;
    H = cat(1,hti.*hir, htr);
end

function rayleigh_channel = get_rayleigh_channel(dist, pl, L,rng_val_)
    % Generates rayleigh fading channel through the IRS
    %

    rng(rng_val_);
    real_part = randn(L,1);
    rng(rng_val_+1);
    img_part = randn(L,1);
    rayleigh_channel = sqrt(1/(2*dist^pl))*(real_part+1i*img_part);

end

function rayleigh_direct_channel = get_rayleigh_direct_channel(dist, pl, rng_val_)
    % Generates rayleigh fading channel for the direct path
    %

    rng(rng_val_)
    rayleigh_direct_channel = sqrt(1/(2*dist^pl))*(randn(1)+1i*randn(1));
end

function ricean_channel = get_ricean_channel(dist, pl, L,rng_val_)
    % Generates rayleigh fading channel through the IRS
    %
    kappa = 5;
    g1 = sqrt(kappa/(1+kappa));
    g2 = sqrt(1/(2*(1+kappa)));

    rng(rng_val_);
    real_part = randn(L,1);
    rng(rng_val_+1);
    img_part = randn(L,1);
    cg = real_part+1i*img_part;
    ricean_channel = remove_small(sqrt(1/(dist^pl))*(g1+g2*cg));

end

function ricean_direct_channel = get_ricean_direct_channel(dist, pl, rng_val_)
    % Generates ricean fading channel for the direct path
    %
    kappa = 10;
    g1 = sqrt(kappa/(1+kappa));
    g2 = sqrt(1/(2*(1+kappa)));

    rng(rng_val_)
    cg = (randn(1)+1i*randn(1));
    ricean_direct_channel = remove_small(sqrt(1/(dist^pl))*(g1+g2*cg));
end


function  removed = remove_small(X)
    indices = abs(X)<=(10^(-10));
    X(indices) = 0 + 0*1i;
    removed = X;
end