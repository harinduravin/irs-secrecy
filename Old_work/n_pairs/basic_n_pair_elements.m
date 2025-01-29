% Initializing number of pairs
global n;
n = 5;

% Position of all the nodes
% n = 5
coords = [-24.48 16.97;-23.93 11.52;-22.36 -7.87;
-16.70 -10.30;-3.39 -2.34;1.74 -2.41;
16.44 -6.14;20.99 -2.66;23.06 18.20;
21.79 12.58;-0.83 18.29;2.15 13.24;
];

% n = 2
% coords = [-24.48 16.97;-23.93 11.52;-22.36 -7.87;
% -16.70 -10.30;-0.83 18.29;2.15 13.24;
% ];

% n = 3
% coords = [-24.48 16.97;-23.93 11.52;-22.36 -7.87;
% -16.70 -10.30;-3.39 -2.34;1.74 -2.41;
% -0.83 18.29;2.15 13.24;];

% n = 4
% coords = [-24.48 16.97;-23.93 11.52;-22.36 -7.87;
% -16.70 -10.30;-3.39 -2.34;1.74 -2.41;
% 16.44 -6.14;20.99 -2.66;-0.83 18.29;2.15 13.24;
% ];

% Preparing the Euclidean distance matrix
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

num_seeds = 105 - 1;
f = waitbar(0,'Please wait...');


global seed;
for seed = 0:num_seeds
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

        leg_inf_stacks = get_stacks("leg_inf");
        leg_inf_with_opp_stacks = get_stacks("leg_inf_with_opp");
        eves_inf_with_self_stacks = get_stacks("eves_inf_with_self");
        eves_stack = eves_inf_all();

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        outer_iter = 2;

        for j = 1:outer_iter
            % Power optimization
            % All the vertices in the (2N)-D box are checked to find the max-min power
            % setting
            %
            max_min = 0;
            max_min_ind = 0;

            for i = 0:(2^(2*n)-1)

                power_string = num2str(dec2bin(i,2*n));
                P_array = get_power_values(power_string,max_power,min_power);
                min_value = get_min_rate(P_array);

                if (min_value > max_min)

                    max_min = min_value;
                    max_min_ind = i;
                end
            end

            power_string = num2str(dec2bin(max_min_ind,2*n));
            if j > 1
                if power_string_old == power_string
                    continue;
                end
            end
            P_array = get_power_values(power_string,max_power,min_power);
            power_string_old = power_string;


            % Reflection phase optimization
            % inner_iter times Successive Convex Approximation (SCA)
            %
            inner_iter = 15;

            for j = 1:inner_iter

                try

                    cvx_begin quiet
                    cvx_solver Mosek
                    variable X(L+1,L+1) complex semidefinite %symmetric
                    variable t

                    minimize t

                    subject to

                        for u = 1:length(user_list)

                            -log(real(trace(get_cvx_leg_inf(user_list(u),P_array)*X))+sigma_ab+sigma_loop) -log(real(trace(get_cvx_eve_inf(user_list(u),P_array)    *X)) +sigma_c) - real(trace(get_grad_S(user_list(u),P_array,W)*(X-W))) <= t;
                        end

                        diag(X) == 1;
                        % norm((X-W),1)<=2*exp(-outer_iter/2);
                        norm((X-W),1)<=2;
                    cvx_end
                catch

                    fprintf("Error occured")
                    seed
                    seed_flag = 1;
                    break;
                end

                W = X;


            end % end of inner iteration
            if seed_flag == 1
                break;
            end
        end % end of outer iteration
        if seed_flag == 1
            break;
        end

        min_value = get_min_rate(P_array);
        min_list = [min_list, min_value];

    end % end of L iteration
    min_master_list = cat(1,min_master_list,min_list);
end % end of random seed iteration

close(f)

figure(1)
plot(L_list,mean(min_master_list));
xlabel('Number of elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cvx_leg_inf = get_cvx_leg_inf(user, power_list)
    % Calculating legitimate user interference terms
    %

    user_id = get_index(user);
    global leg_inf_with_opp_stacks;
    leg_power_with_opp = filter_power_values(char(user), power_list,"leg_inf_with_opp");
    cvx_leg_inf = dot_product(leg_power_with_opp', leg_inf_with_opp_stacks(:,:,:,user_id));

end

function cvx_eve_inf = get_cvx_eve_inf(user, power_list)
    % Calculating evesdropper interference terms
    %

    user_id = get_index(user);
    global eves_inf_with_self_stacks;
    eves_power_with_self = filter_power_values(char(user), power_list,"eves_inf_with_self");
    cvx_eve_inf = dot_product(eves_power_with_self',eves_inf_with_self_stacks(:,:,:,user_id));
    
end

function grad_S = get_grad_S(user, power_list, W)
    % Obtaining grad_S for Taylor approximation lower bound in the CVX implementation.
    % The return matrix depends on the specified user.
    %

    user_id = get_index(user);

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);

    global leg_inf_stacks;
    global eves_stack;

    eves_power_all = filter_power_values(char(user), power_list,"eves_inf_all");
    leg_power = filter_power_values(char(user), power_list,"leg_inf");

    leg_inf_term1 = dot_product(leg_power', leg_inf_stacks(:,:,:,user_id));
    eves_inf_term2 = dot_product(eves_power_all',eves_stack);

    term_1 = leg_inf_term1/(trace(leg_inf_term1*W) + sigma_ab + sigma_loop);
    term_2 = eves_inf_term2/(trace(eves_inf_term2*W) + sigma_c);

    grad_S = -(term_1+term_2)/log(2);
    
end

function stacks = get_stacks(stack_type)
    % Collecting stacks of channels into 4-D array for later use
    % according to stack_type specified
    %

    global user_list;

    temp_stack = [];

    for i = 1:length(user_list)
        if stack_type == "leg_inf"
            temp_stack = cat(4,temp_stack,leg_inf(char(user_list(i))));
        end
        if stack_type == "leg_inf_with_opp"
            temp_stack = cat(4,temp_stack,leg_inf_with_opp(char(user_list(i))));
        end
        if stack_type == "eves_inf_with_self"
            temp_stack = cat(4,temp_stack,eves_inf_with_self(char(user_list(i))));            
        end
    end

    stacks = temp_stack;
end


function power_array = get_power_values(string_val, dbmax, dbmin)
    % This function generates a list of power values from the representative string
    % Eg-: '1001001'  --> [pmax pmin pmin pmax pmin pmin pmax]
    % 1 is replaced by pmax, 0 is replaced by pmin. pmin and pmax are specified in dBm.

    % Converting from dBm to watt
    pmax = dbm2watt(dbmax);
    pmin = dbm2watt(dbmin);

    power_list = [];

    % Generating the list
    for r = 1:length(string_val)
        if string_val(r) == '0' 
            power_list = [power_list pmin];
        else
            power_list = [power_list pmax];
        end
    end
    
    power_array = power_list;

end

function min_rate = get_min_rate(power_list)
    % Obtains a list of data rates for all users for a given power setting.
    % After obtaining the list, the minimum of all rates is obtained.
    %

    global user_list;
    rate_list = [];

    % Looping through all users
    for i = 1:length(user_list)
        rate_list = [rate_list get_rate(user_list(i), power_list)];
    end

    rate_list;

    % Returning the minimum rate
    min_rate = min(rate_list);

end

function rate = get_rate(user, power_list)
    % Returns the datarate achieved by a given user, at a given power setting.
    % This requires retrieving the channels stacks stored after generation.
    % The dot_product function is used to obtain the interference terms.
    %

    user_id = get_index(user);

    global W;
    global leg_inf_with_opp_stacks;
    global eves_inf_with_self_stacks;
    global leg_inf_stacks;
    global eves_stack;

    % Choose filter modes from ["leg_inf" "eves_inf_with_self" "leg_inf_with_opp" "eves_inf_all"]
    % Retrieving suitable power arrays for dot_product function
    eves_power_with_self = filter_power_values(char(user), power_list,"eves_inf_with_self");
    leg_power_with_opp = filter_power_values(char(user), power_list,"leg_inf_with_opp");
    eves_power_all = filter_power_values(char(user), power_list,"eves_inf_all");
    leg_power = filter_power_values(char(user), power_list,"leg_inf");

    leg_inf_term1 = dot_product(leg_power', leg_inf_stacks(:,:,:,user_id));
    leg_inf_term2 = dot_product(leg_power_with_opp', leg_inf_with_opp_stacks(:,:,:,user_id));
    eves_inf_term1 = dot_product(eves_power_with_self',eves_inf_with_self_stacks(:,:,:,user_id));
    eves_inf_term2 = dot_product(eves_power_all',eves_stack);

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);

    eves_inf = log2(1 + trace((eves_inf_term2 - eves_inf_term1)*W)/(trace(eves_inf_term1*W) + sigma_c));

    leg_inf = log2(1 + trace((leg_inf_term2 - leg_inf_term1)*W)/(trace(leg_inf_term1*W) + sigma_ab + sigma_loop));
    rate = max(0,real(leg_inf - eves_inf));

end


% TODO test this function properly.
function y = dot_product(power, channel_cube)
    % Here ordered lists of power values and stacks of 2D channels are multiplied
    % 3D and 1D matrices are multiplied together and added
    % Returns a 2D matrix (L+1 x L+1)
    
    y = permute(pagemtimes(permute(channel_cube,[1 3 2]),power),[1 3 2]);
end

function filtered_power = filter_power_values(user, power_list, filter_type)
    % Power_list contains all the power values used for transmission at all users.
    % This function can filter and provide ordered lists of these values in a required 
    % order.  This order should match with the channel stacks that are created before
    % looping the algorithm.
    %

    global n;
    node_type = user(1);
    pair_num = str2num(user(2));

    eligible_list = [];

    % Creating a list of eligible users as per the requirement
    for j = 1:n

        if filter_type ~= "eves_inf_all"
            if j == pair_num
                continue
            end
        end
    
        string1 = string(strcat('A',num2str(j)));
        eligible_list = [eligible_list string1];
        string2 = string(strcat('B',num2str(j)));
        eligible_list = [eligible_list string2];
    end

    if filter_type == "leg_inf_with_opp"

        if node_type == 'A'
            opp_node = 'B';
        else
            opp_node = 'A';
        end

        opposite = string(strcat(opp_node,user(2)));
        eligible_list = [eligible_list opposite];
    end

    if filter_type == "eves_inf_with_self"

        eligible_list = [eligible_list string(user)];
    end

    indices_list = convert2indices(eligible_list);

    % Adding the power values to the list
    temp_list  = [];
    for q = 1:length(eligible_list)

        temp_list = [temp_list power_list(indices_list(q))];
    end

    filtered_power = temp_list;
end

function power = dbm2watt(dbm_value)
    % Converting power values from dBm to Watt for calculations
    %

    power = (10^(dbm_value/10))*(10^(-3));
end

function D = sqDistance(X, Y)
    % Obtaining the Euclidean distance matrix from given x,y coordinates.
    %

    D = sqrt(bsxfun(@plus,dot(X,X,2)',dot(Y,Y,2))-2*(X*Y'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Below functions are used to collect stacks of channels used for
% calculating the interference terms. 
%

function eves_inf = eves_inf_all()
    % Stacks all the channels from legitimate users to the eavesdroppers.
    %

    global n;

    eligible_eves_list = [];

    % Looping through all n pairs
    for j = 1:n

        % Considering A nodes
        string1 = string(strcat('A',num2str(j)));
        eligible_eves_list = [eligible_eves_list string1];

        % Considering B nodes
        string2 = string(strcat('B',num2str(j)));
        eligible_eves_list = [eligible_eves_list string2];

    end

    indices_list = convert2indices(eligible_eves_list);
    eves_inf = get_eves_inf(indices_list);

end

function eves_inf = eves_inf_with_self(node_name)
    % Stacks all the channels from legitimate users to the eavesdroppers, 
    % except for the opposite user of a given node.
    % 

    global n;

    node_type = node_name(1);
    pair_num = str2num(node_name(2));

    eligible_eves_list = [];

    % Looping through all n pairs
    for j = 1:n

        % Skipping the current pair of the given node
        if j == pair_num
            continue
        end

        string1 = string(strcat('A',num2str(j)));
        eligible_eves_list = [eligible_eves_list string1];
        string2 = string(strcat('B',num2str(j)));
        eligible_eves_list = [eligible_eves_list string2];

    end

    % Adding back the current node
    eligible_eves_list = [eligible_eves_list string(node_name)];

    indices_list = convert2indices(eligible_eves_list);
    eves_inf = get_eves_inf(indices_list);

end

function leg_inf = leg_inf_with_opp(node_name)
    % Stacks all the channels from legitimate users to the given node, 
    % 

    global n;

    node_type = node_name(1);
    pair_num = str2num(node_name(2));

    eligible_leg_list = [];

    for j = 1:n

        if j == pair_num
            continue
        end

        string1 = string(strcat('A',num2str(j)));
        eligible_leg_list = [eligible_leg_list string1];
        string2 = string(strcat('B',num2str(j)));
        eligible_leg_list = [eligible_leg_list string2];

    end

    if node_type == 'A'
        opp_node = 'B';
    else
        opp_node = 'A';
    end

    if node_type == 'A'
        node_index = 2*pair_num -1;
    else 
        node_index = 2*pair_num;
    end

    % Adding the opposite node the list
    opposite = string(strcat(opp_node,node_name(2)));
    eligible_leg_list = [eligible_leg_list opposite];

    indices_list = convert2indices(eligible_leg_list);
    leg_inf = get_leg_inf(node_index,indices_list);

end

function leg_inf = leg_inf(node_name)
    % Stacks all the channels from legitimate users to the given node, 
    % except for the opposite node of the current pair.
    % 

    global n;

    node_type = node_name(1);
    pair_num = str2num(node_name(2));

    eligible_leg_list = [];

    for j = 1:n

        if j == pair_num
            continue
        end

        string1 = string(strcat('A',num2str(j)));
        eligible_leg_list = [eligible_leg_list string1];
        string2 = string(strcat('B',num2str(j)));
        eligible_leg_list = [eligible_leg_list string2];

    end

    if node_type == 'A'
        node_index = 2*pair_num -1;
    else 
        node_index = 2*pair_num;
    end

    indices_list = convert2indices(eligible_leg_list);
    leg_inf = get_leg_inf(node_index,indices_list);

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

    gti = get_rayleigh_channel(dti, pl, L, 5*rng_val);
    hti = sqrt(Lo)*gti;
    gir = get_rayleigh_channel(dir, pl, L, 5*rng_val+2);
    hir = sqrt(Lo)*gir;
    gtr = get_rayleigh_direct_channel(dtr, pl, 5*rng_val+4);
    htr = sqrt(Lo)*gtr;
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
