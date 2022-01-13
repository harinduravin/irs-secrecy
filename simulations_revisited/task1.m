clear;
N = 10^4;
N2 = 100;

ES = logspace(-1,1,10); % SNR values from -10dB to +10dB

% Initialize variables for rayleigh BER 
BER_Rayleigh_E = zeros(length(ES));
BER_Rayleigh_T = zeros(length(ES));

% Initialize variables for noise BER 
BER_Noise_E = zeros(length(ES));
BER_Noise_T = zeros(length(ES));

f = waitbar(0,'Please wait...');
for e = 1:length(ES)
    Es = ES(e);
    BER_Rayleigh = 0;
    BER_Noise = 0;

    for i = 1:N2
        Data = randi(1,1,N); % N random data bits
        T_symbol = sqrt(Es).*(2.*Data-1); % Mapping to BPSK transmit symbols

        h = sqrt(1/2)*(randn(1,N)+1i*(randn(1,N))); % rayleigh fading channel
        n = sqrt(1/2)*(randn(1,N)+1i*(randn(1,N))); % complex noise

        y = T_symbol.*h + n;
        y_n = T_symbol + n;

        z = conj(h).*y; % Coherent Detection
        z_n = y_n;

        decision = z>0;
        BER_Rayleigh = BER_Rayleigh + mean(abs(decision-Data));

        decision_n = z_n>0;
        BER_Noise = BER_Noise + mean(abs(decision_n-Data));

        waitbar(e/length(ES),f,'Simulating...');

    end

    BER_Rayleigh_E(e) = BER_Rayleigh/N2;
    BER_Rayleigh_T(e) = (1/2)*(1-sqrt(Es/(Es+1)));

    BER_Noise_E(e) = BER_Noise/N2;
    BER_Noise_T(e) = Q(sqrt(2*Es));

end

close(f)
semilogy(10*log10(ES),BER_Rayleigh_E,'-o',10*log10(ES),BER_Noise_E,'-o',10*log10(ES),BER_Rayleigh_T,'--o',10*log10(ES),BER_Noise_T,'--o');
% hold on;
% semilogy(10*log10(ES),BER_Noise_E);

% Q function implementation
function out=Q(x)
out=0.5.*erfc(x/sqrt(2));
end