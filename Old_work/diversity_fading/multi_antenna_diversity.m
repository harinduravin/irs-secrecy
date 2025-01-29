N = 100000;

%Number of antennas
M = 1;

% Average bit error rate for BPSK
ebnodB = 0:1:10;
ebno = 10.^(ebnodB./10);

figure;
Pbr = array_channel(M,N);
semilogy(ebnodB, Pbr)
hold on
M = 2;
Pbr = array_channel(M,N);
semilogy(ebnodB, Pbr)
M = 4;
Pbr = array_channel(M,N);
semilogy(ebnodB, Pbr)
legend('M = 1', 'M = 2', 'M = 4')

% Rayleigh fading channel tx-->rx
function Pb = array_channel(M,N)
ebnodB = 0:1:10;
ebno = 10.^(ebnodB./10);
g = sqrt(1/2)*(randn(M,N)+1i*randn(M,N));
g_norm_s = sum(abs(g).^2,1);
ebnor = ebno.'*(g_norm_s);
Pb = mean(Q(sqrt(2*ebnor)),2);
end

% Q function implementation
function out=Q(x)
out=0.5.*erfc(x/sqrt(2));
end
