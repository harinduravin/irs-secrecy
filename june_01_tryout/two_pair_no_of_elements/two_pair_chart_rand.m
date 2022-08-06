
% Initializing number of pairs
global n;
n = 2;

% Position of all the nodes
% n = 2

coords = [50.13 76.03;31.56 40.60;-85.14 -38.49;
-88.72 -78.33;
0.00 97.50;0.00 -97.50];

% Preparing the Euclidean distance matrix
global dist;
dist = sqDistance(coords, coords);

L_list = 5:5:40;

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
all_master_list = [];

num_seeds = 40 - 1;
f = waitbar(0,'Please wait...');


global seed;
for seed = 0:num_seeds

    seed
    seed_flag = 0;

    waitbar(seed/num_seeds,f,'Simulating...');
    min_list = [];
    all_list = [];

    % Number of elements in the IRS
    global L;
    for L = 5:5:40
        fprintf('\n')
        fprintf(strcat("starting with ",string(L), " IRS eleements."))
        fprintf('\n')

        % IRS reflection matrix (Initialized randomly)
        wi = exp(1i*(2*rand(L,1)-1)*pi);
        w = [wi ; 1].';

        global W; % Check this again
        % W = remove_small(w'*w);
        W = w'*w;
%         temp_W = zeros(size(W));

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
        outer_iter = 10;

        upper_t = 25;
        lower_t = 0;
        t = upper_t;

        % upper_z = 25;
        % lower_z =0;
        % z = upper_z;

        power_string = num2str(dec2bin(0,2*n));
        P_hat = (get_power_values(power_string,max_power,min_power))';

        for j = 1:outer_iter

            %
            % Power optimization - Successive convex approximation
            %
            % fprintf('\n')
            % fprintf('l :')

            try

                cvx_begin quiet
                cvx_solver Mosek
                variable P(2*n,1)
                % variable z

                minimize z

                subject to

                    for u = 1:length(user_list)

                        % for n = 2
                        20-log(sigma_ab + sigma_loop + P(p_(user_list(u),1,1))*H_(user_list(u),1,1)+ P(p_(user_list(u),1,2))*H_(user_list(u),1,2) + P(p_(user_list(u),1,3))*H_(user_list(u),1,3))...
                        -log(sigma_c + P(p_(user_list(u),2,1))*H_(user_list(u),2,1)+ P(p_(user_list(u),2,2))*H_(user_list(u),2,2) + P(p_(user_list(u),2,3))*H_(user_list(u),2,3))...
                        -get_S(user_list(u),P_hat)...
                        -(get_grad_P(user_list(u),P_hat,W))'*(P-P_hat) <= z;

                    end

                    for i = 1:length(user_list)

                        P(i) <= dbm2watt(max_power); 
                        P(i) >= dbm2watt(min_power);

                    end
                    norm((P-P_hat),1)<=1;
                cvx_end

            catch
                fprintf("skipped1")
                cvx_status
                seed;
                seed_flag = 1;
                break;
            end

            P_hat = P;

            if seed_flag == 1
                break;
            end

            fprintf('\n')
            fprintf('P values :')
            for x = 1:(2*n)
                fprintf('|')
                fprintf(string(P_hat(x)))
            end

            % fprintf(' z: ')
            % fprintf(string(z))


            fprintf('\n')
            fprintf('m :')



%             % Reflection phase optimization
%             % Successive Convex Approximation (SCA)
%             %
% 
%             for m = 1:100
% 
% 
%                 fprintf(string(m))
%                 fprintf('|')
%                 try
%                     cvx_begin quiet
%                     % cvx_precision high
%                     cvx_solver Mosek
%                     variable X(L+1,L+1) complex semidefinite %symmetric
%     
%                     minimize 0
%     
%                     subject to
% 
%                         for u = 1:length(user_list)
%     
%                             % -log(real(trace(get_cvx_leg_inf(user_list(u),P_hat)*X))+sigma_ab+sigma_loop) -log(real(trace(get_cvx_eve_inf(user_list(u),P_hat)*X)) +sigma_c) - real(trace(get_grad_S(user_list(u),P_hat,W)*(X-W))) <= t;
%                             10-log(real(trace(remove_small(get_cvx_leg_inf(user_list(u),P_hat))*X))+sigma_ab+sigma_loop) -log(real(trace(remove_small(get_cvx_eve_inf(user_list(u),P_hat))*X)) +sigma_c)-get_S(user_list(u),P_hat) - real(trace(remove_small(get_grad_S(user_list(u),P_hat,W))*(X-remove_small(W)))) <= t;
%                         end
%     
%                         diag(X) == 1;
%                         % norm((X-W),1)<=2*exp(-outer_iter/2);
%                         % norm((X-W),1)<=8;
%                     cvx_end
% 
%                     cvx_status;
%     
%                 catch
%                     fprintf('skipped')
%                     cvx_status
%                     seed;
%                     seed_flag = 1;
%                     break;
%                 end
% 
%                 if string(cvx_status) == string('Infeasible')
%                     lower_t = t;
% 
%                 elseif string(cvx_status) == string('Failed')
%                     lower_t = t;
%                     fprintf('Failed')
%                     % seed_flag = 1;
%                     % break;
% 
% 
%                 else
%                     upper_t = t;
%                     if t < 0
%                         lower_t = (1.3)*t;
%                     else
%                         lower_t = (0.7)*t;
%                     end
% 
%                     temp_W = X;
% 
%                 end
% 
%                 t = (lower_t + upper_t)/2;
% 
%                 if (upper_t - lower_t) < 0.2
% %                     fprintf('reached IRS');
%                     t;
%                     break;
%                 end
% 
%             end
%             W = temp_W;
% 
%             fprintf(' t: ')
%             fprintf(string(t))
%             if seed_flag == 1
%                 break;
%             end


        end % end of outer iteration
        if seed_flag == 0
            [min_value, rate_list] = get_min_rate(P_hat);
        elseif seed_flag == 1
            min_value = nan;
            seed_flag = 0;
        end


        min_list = [min_list, min_value];
        all_list = cat(1,all_list, rate_list);

    end % end of L iteration
    min_master_list = cat(1,min_master_list,min_list);
    all_master_list = cat(3,all_master_list, all_list);
end % end of random seed iteration

mean_user_rates = mean(all_master_list,3);


close(f)

figure(1)
plot(L_list,mean(min_master_list,'omitnan'));
hold on;
plot(L_list,mean_user_rates(:,1)','DisplayName','A1');
plot(L_list,mean_user_rates(:,2)','DisplayName','B1');
plot(L_list,mean_user_rates(:,3)','DisplayName','A2');
plot(L_list,mean_user_rates(:,4)','DisplayName','B2');


ylim([-1 5])
xlabel('Number of elements')
ylabel('Minimum secrecy Rate(bits/sec/Hz)')
legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function S = get_S(user,P_hat)
    % Calculation of the S term in the constraint functions. Calculating this value is important
    % for the purpose of carrying out max-min optimization. All users' constraint functions should include constants
    % even though they do not affect the optimization variables.
    %

    global W;

    user = char(user);
    user_id = get_index(user);

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);


    global leg_inf_stacks;
    global eves_stack;

    eves_power_all = filter_power_values(char(user), P_hat,"eves_inf_all");
    leg_power = filter_power_values(char(user), P_hat,"leg_inf");

    leg_inf_term1 = dot_product(leg_power', leg_inf_stacks(:,:,:,user_id));
    eves_inf_term2 = dot_product(eves_power_all',eves_stack);

    term_1 = real(trace(leg_inf_term1*W) + sigma_ab + sigma_loop);
    term_2 = real(trace(eves_inf_term2*W) + sigma_c);

    S = -log(term_1) -log(term_2);

end


function index = p_(user,set_num,term_num)
    % This functions enables interactive use of power optimization variable in the cvx max-min type
    % optimization problem. This simplifies the cvx objective function.
    %

    user = char(user);

    global n;
    node_type = user(1);
    pair_num = str2num(user(2));

    eligible_list = [];

    % Creating a list of eligible users as per the requirement
    for j = 1:n

        if j == pair_num
            continue
        end

        string1 = string(strcat('A',num2str(j)));
        eligible_list = [eligible_list string1];
        string2 = string(strcat('B',num2str(j)));
        eligible_list = [eligible_list string2];
    end

    if set_num == 1 % filter_type == "leg_inf_with_opp"

        if node_type == 'A'
            opp_node = 'B';
        else
            opp_node = 'A';
        end

        opposite = string(strcat(opp_node,user(2)));
        eligible_list = [eligible_list opposite];
    end

    if set_num == 2 % filter_type == "eves_inf_with_self"

        eligible_list = [eligible_list string(user)];
    end

    indices_list = convert2indices(eligible_list);
    index = indices_list(term_num);

end

function channel = H_(user,set_num,term_num)
    % This function works similar to P_ function above. Instead of the power variable index,
    % this function returns the suitable channel for the expected term.
    %

    global W;

    user = char(user);

    user_id = get_index(user);
    global leg_inf_with_opp_stacks;
    global eves_inf_with_self_stacks;
    global n;
    node_type = user(1);
    pair_num = str2num(user(2));

    if set_num == 1 % filter_type == "leg_inf_with_opp"

        channel = real(trace(leg_inf_with_opp_stacks(:,:,term_num,user_id)*W));  % Check again

    end

    if set_num == 2 % filter_type == "eves_inf_with_self"

        channel = real(trace(eves_inf_with_self_stacks(:,:,term_num,user_id)*W));  % Check again

    end

end

function grad_P = get_grad_P(user, P_hat,W)
    % This function returns a column of gradient values of P vector. (grad P)
    %

    user = char(user);
    user_id = get_index(user);
    pair_num = str2num(user(2));

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);

    global user_list;
    global leg_inf_stacks;
    global eves_stack;

    eves_power_all = filter_power_values(char(user), P_hat,"eves_inf_all");
    leg_power = filter_power_values(char(user), P_hat,"leg_inf");

    leg_inf_term1 = dot_product(leg_power', leg_inf_stacks(:,:,:,user_id));
    eves_inf_term2 = dot_product(eves_power_all',eves_stack);

    term_1 = (trace(leg_inf_term1*W) + sigma_ab + sigma_loop);
    term_2 = (trace(eves_inf_term2*W) + sigma_c);

    grad_list = [];
    counter = 0;

    for i = 1:length(user_list) 
        value = trace(eves_stack(i)*W)/term_2;

        if ~((i == 2*pair_num) | (i == 2*pair_num-1))
            counter = counter + 1;
            value = value + trace(leg_inf_stacks(:,:,counter,user_id)*W)/term_1;
        end

        grad_list = [grad_list real(value)];
    end

    grad_P = grad_list';

end

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

    grad_S = -(term_1+term_2);
    
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

function [min_rate, rates] = get_min_rate(power_list)
    % Obtains a list of data rates for all users for a given power setting.
    % After obtaining the list, the minimum of all rates is obtained.
    %

    global user_list;
    rate_list = [];

    % Looping through all users
    for i = 1:length(user_list)
        rate_list = [rate_list get_rate(user_list(i), power_list)];
    end

    rates = rate_list;

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
    rate = real(leg_inf - eves_inf);

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
    % Using 3GPP Urban Micro pathloss models; A frequency of 750 MHz is assumed.


    % LOS channel : Used for Users to IRS channels
    % Lo2 = 10^(2*-2.55012);
    Lo2 = 0.001;
    % NLOS channel : Used for direct path
    % Lo3 = 10^(2*-1.94515);
    Lo3 = 0.001;

    gti = get_ricean_channel(dti, 2, L, 5*rng_val);
    gir = get_ricean_channel(dir, 2, L, 5*rng_val+2); %2.2
    gtr = get_ricean_direct_channel(dtr, 4, 5*rng_val+4); %3.67
    htr = sqrt(Lo3)*gtr;
    H = cat(1,sqrt(Lo2)*gti.*gir, htr);
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
    kappa =20;
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
    kappa = 20;
    g1 = sqrt(kappa/(1+kappa));
    g2 = sqrt(1/(2*(1+kappa)));

    rng(rng_val_)
    cg = (randn(1)+1i*randn(1));
    ricean_direct_channel = remove_small(sqrt(1/(dist^pl))*(g1+g2*cg));
end


function  removed = remove_small(X)
    indices = abs(X)<=(10^(-7));
    X(indices) = 0 + 0*1i;
    removed = X;
end