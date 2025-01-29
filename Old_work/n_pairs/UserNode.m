global vari;
for vari = 1:5
    test_func(3)
end

function ret = test_func(q)
    global vari;
    ret = q*vari;
end