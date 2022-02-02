rng('default')
X = rand(3,2);

coords = [-2.24 6.18;4.59 7.21;-23.78 5.05;
-18.59 8.85;-12.37 -4.52;-14.86 -9.01;
23.38 10.60;17.39 13.63;20.16 -3.91;
21.84 -10.31;-0.09 -16.47;5.98 -6.85;
];

global dist;
dist = sqDistance(coords, coords)
global n;
n = 5;

% leg_H_stack = leg_inf_with_opp('B4');
eves_H_stack_self = eves_inf_with_self('A3');

function D = sqDistance(X, Y)
    D = sqrt(bsxfun(@plus,dot(X,X,2)',dot(Y,Y,2))-2*(X*Y'));
end

function eves_inf = eves_inf_with_self(node_name)

    global n;

    node_type = node_name(1);
    pair_num = str2num(node_name(2));

    eligible_eves_list = [];

    for j = 1:n

        if j == pair_num
            continue
        end

        string1 = string(strcat('A',num2str(j)));
        eligible_eves_list = [eligible_eves_list string1];
        string2 = string(strcat('B',num2str(j)));
        eligible_eves_list = [eligible_eves_list string2];

    end


    eligible_eves_list = [eligible_eves_list string(node_name)];

    eligible_eves_list

    indices_list = convert2indices(eligible_eves_list)
    eves_inf = get_eves_inf(indices_list);

end

function leg_inf = leg_inf_with_opp(node_name)

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

    opposite = string(strcat(opp_node,node_name(2)));
    eligible_leg_list = [eligible_leg_list opposite];

    eligible_leg_list

    indices_list = convert2indices(eligible_leg_list)
    leg_inf = get_leg_inf(node_index,indices_list);

end

function leg_inf = leg_inf(node_name)

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

    eligible_leg_list

    indices_list = convert2indices(eligible_leg_list);
    leg_inf = get_leg_inf(node_index,indices_list);

end

function leg_inf = get_leg_inf(node_index,indices_list)

    stack = [];

    for p = 1:length(indices_list)
        stack = cat(3, stack,get_leg_channel(node_index,indices_list(p)));
    end

    leg_inf = stack;

end

function eves_inf = get_eves_inf(indices_list)

    stack = [];

    for p = 1:length(indices_list)
        stack = cat(3, stack,get_eves_channel(indices_list(p)));
    end

    eves_inf = stack;

end

function leg_channel = get_leg_channel(rx_index,tx_index)

    global dist;
    global n;
    L = 10;
    Lo = 10^(-3);
    pl = 2;

    irs_index = 2*n+1;
    dti = dist(tx_index, irs_index);
    dir = dist(irs_index, rx_index);
    dtr = dist(rx_index,tx_index);

    H = get_H(dti,dir,dtr,L,Lo,pl,1); % TODO Change this line to complete the randomness

    leg_channel = H*H';

end

function eves_channel = get_eves_channel(tx_index)

    global dist;
    global n;
    L = 10;
    Lo = 10^(-3);
    pl = 2;

    irs_index = 2*n+1;
    rx_index = 2*n+2;
    dti = dist(tx_index, irs_index);
    dir = dist(irs_index, rx_index);
    dtr = dist(rx_index,tx_index);

    H = get_H(dti,dir,dtr,L,Lo,pl,1); % TODO Change this line to complete the randomness

    eves_channel = H*H';

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

function H = get_H(dti,dir,dtr,L,Lo,pl,rng_val)

    gti = get_rayleigh_channel(dti, pl, L, 5*rng_val);
    hti = sqrt(Lo)*gti;
    gir = get_rayleigh_channel(dir, pl, L, 5*rng_val+2);
    hir = sqrt(Lo)*gir;
    gtr = get_rayleigh_direct_channel(dtr, pl, 5*rng_val+4);
    htr = sqrt(Lo)*gtr;
    H = cat(1,hti.*hir, htr);
end

function rayleigh_channel = get_rayleigh_channel(dist, pl, L,rng_val_)

    rng(rng_val_);
    real_part = randn(L,1);
    rng(rng_val_+1);
    img_part = randn(L,1);
    rayleigh_channel = sqrt(1/(2*dist^pl))*(real_part+1i*img_part);
%     rayleigh_channel = sqrt(1/(2*dist^pl))*(randn(L,1)+1i*randn(L,1));
end

function rayleigh_direct_channel = get_rayleigh_direct_channel(dist, pl, rng_val_)
    rng(rng_val_)
    rayleigh_direct_channel = sqrt(1/(2*dist^pl))*(randn(1)+1i*randn(1));
end