leg_inf_with_opp('B1')

function leg_inf = leg_inf_with_opp(node_name)

    n = 5;

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

    opposite = string(strcat(opp_node,node_name(2)));
    eligible_leg_list = [eligible_leg_list opposite];

    leg_inf = get_leg_inf(node_name,eligible_leg_list);

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

    indices_list = convert2indices(node_list);
    leg_inf = get_leg_inf(node_index,indices_list);

end

function leg_inf = get_leg_inf(node_index,indices_list)

    sum = 0;

    for p = 1:length(indices_list)
        sum = sum + get_leg_term(node_index,indices_list(p))
    end

end

function leg_term = get_leg_term(rx_index,tx_index)

    global dist;
    global n;

    irs_index = 2*n+1;
    dti = dist(rx_index, irs_index);
    dir = dist(irs_index, tx_index);
    dtr = dist(rx_index,tx_index);

    leg_term = dti + dir + dtr;

end

function indices_list = convert2indices(node_list)

    initial_list = [];
    for p = 1:length(node_list)

        node_type = (node_list(p))(1);
        pair_num = str2num((node_list(p))(2));

        if node_type == 'A'
            index = 2*pair_num -1;
        else 
            index = 2*pair_num;
        end

        initial_list = [initial_list index];
    end
    indices_list = initial_list;

end

function H = get_H(dti,dir,dtr,wave_l,L,Lo,pl,rng_val)

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
