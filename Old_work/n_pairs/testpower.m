global n;
n =5;
max_min_ind = 7;

max_power = 15;
min_power = 0;
power_string = num2str(dec2bin(max_min_ind,2*n));
P_array = get_power_values(power_string,max_power,min_power);

filtered_array = filter_power_values('A1', P_array, "eves_inf_with_self") % "leg_inf" "eves_inf_with_self" "leg_inf_with_opp" "eves_inf_all"

function power_array = get_power_values(string_val, dbmax, dbmin)

    pmax = dbm2watt(dbmax);
    pmin = dbm2watt(dbmin);

    power_list = [];

    for r = 1:length(string_val)
        if string_val(r) == '0' 
            power_list = [power_list pmin];
        else
            power_list = [power_list pmax];
        end
    end
    
    power_array = power_list;

end

function power = dbm2watt(dbm_value)

    power = (10^(dbm_value/10))*(10^(-3));
end


function filtered_power = filter_power_values(user, power_list, filter_type)

    global n;
    node_type = user(1);
    pair_num = str2num(user(2));

    eligible_list = [];

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

    temp_list  = [];
    for q = 1:length(eligible_list)

        temp_list = [temp_list power_list(indices_list(q))];

    end

    filtered_power = temp_list;

end

function indices_list = convert2indices(node_list)

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

% if stack_type == "leg_inf"
%     temp_stack = cat(4,temp_stack,leg_inf(char(user_list(i))));
% end
% if stack_type == "leg_inf_with_opp"
%     temp_stack = cat(4,temp_stack,leg_inf_with_opp(char(user_list(i))));
% end
% if stack_type == "eves_inf_with_self"
%     temp_stack = cat(4,temp_stack,eves_inf_with_self(char(user_list(i))));            
% end


% idx=[]; %get the index of all values classified in category 1.
% for i = 1: length(P)
%   if P(i)==1, idx = [idx,i]; end
% end

% for j = 1:n

%     if j == pair_num
%         continue
%     end

%     string1 = string(strcat('A',num2str(j)));
%     eligible_leg_list = [eligible_leg_list string1];
%     string2 = string(strcat('B',num2str(j)));
%     eligible_leg_list = [eligible_leg_list string2];

% end