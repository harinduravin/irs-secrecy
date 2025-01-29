matr = get_nmatrix(5)

function nmatrix = get_nmatrix(n)
    u = zeros(2*n,2*n);
    p = 0;
    for i = 1:2*n
        for j = 1:2*n
            if j>i
                p = p+1;
                u(j,i) = p;
            end

        end
    end
    nmatrix = u + u';

end