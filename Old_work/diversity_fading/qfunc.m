a = -4:0.1:4;
y = Q(a);
figure;
plot(a,y)
grid

% Q function implementation
function out=Q(x)
out=0.5.*erfc(x/sqrt(2));
end