iter = 5;
large_array = [];
filename = 'test2.mat';
for i = 1:iter
    randarray = randn(1,5);
    large_array = cat(1,large_array,randarray);
    save(filename,'large_array')

end

large_array
