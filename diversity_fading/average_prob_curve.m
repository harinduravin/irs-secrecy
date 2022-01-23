N = 100000;

%Number of antennas
M = 1;

% Rayleigh fading channel
g = sqrt(1/2)*(randn(1,N)+1i*randn(1,N));
meang = mean(abs(g).^2);
% Average bit error rate for BPSK
ebnodB = 0:1:10;
ebno = 10.^(ebnodB./10);
Pb = Q(sqrt(2*ebno));

ebnor = ebno.'*(abs(g).^2);
Pbr = mean(Q(sqrt(2*ebnor)),2);
% ebnordB = 10*log10(ebnor);
figure;
semilogy(ebnodB, Pb)
hold on
semilogy(ebnodB, Pbr)



% Q function implementation
function out=Q(x)
out=0.5.*erfc(x/sqrt(2));
end