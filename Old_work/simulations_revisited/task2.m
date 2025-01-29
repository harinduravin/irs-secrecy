clear all;
N = 10^3;
N2 = 1000;

ES = logspace(-1,1,10); % SNR values from -10dB to +10dB

% Initialize variables for rayleigh BER 
BER_Rayleigh_E = zeros(1,length(ES));
BER_Rayleigh_T = zeros(1,length(ES));
BER_Rayleigh_SC = zeros(1,length(ES));
BER_Rayleigh_SC_T = zeros(1,length(ES));
BER_Rayleigh_EG = zeros(1,length(ES));
BER_Rayleigh_EG_T = zeros(1,length(ES));

f = waitbar(0,'Please wait...');
for e = 1:length(ES)
    Es = ES(e);
    BER_Rayleigh = 0;
    BER_SC = 0;
    BER_EG = 0;

    for i = 1:N2
        Data = randi([0 1],1,N); % N random data bits
        T_symbol = sqrt(Es/3).*(2.*Data-1); % Mapping to BPSK transmit symbols

        h1 = sqrt(1/2)*(randn(1,N)+1i*(randn(1,N))); % rayleigh fading channel 1
        h2 = sqrt(1/2)*(randn(1,N)+1i*(randn(1,N))); % rayleigh fading channel 2
        h3 = sqrt(1/2)*(randn(1,N)+1i*(randn(1,N))); % rayleigh fading channel 3

        n1 = sqrt(1/6)*(randn(1,N)+1i*(randn(1,N))); % complex noise (Adds up to No=1 in Eb/N0)
        n2 = sqrt(1/6)*(randn(1,N)+1i*(randn(1,N))); % complex noise
        n3 = sqrt(1/6)*(randn(1,N)+1i*(randn(1,N))); % complex noise

        % Signals received at 3 receiver antennas
        y1 = T_symbol.*h1 + n1;
        y2 = T_symbol.*h2 + n2;
        y3 = T_symbol.*h3 + n3;

        z = conj(h1).*y1 + conj(h2).*y2 + conj(h3).*y3; % 3-MRC receiver

        z_eg = (conj(h1)./abs(conj(h1))).*y1 + (conj(h2)./abs(conj(h2))).*y2 + (conj(h3)./abs(conj(h3))).*y3; % Equal gain combining

        all_channels = [h1.*conj(h1); h2.*conj(h2); h3.*conj(h3)];
        concat_signal = [conj(h1).*y1; conj(h2).*y2; conj(h3).*y3];
        [abs_sc, idx] = max(abs(all_channels),[],1);
        z_sc = zeros(1,N);
        for rows= 1:N
            z_sc(rows) = concat_signal(idx(rows),rows);
        end

        decision = z>0;
        BER_Rayleigh = BER_Rayleigh + mean(abs(decision-Data));

        decision_sc = z_sc>0;
        BER_SC = BER_SC + mean(abs(decision_sc-Data));      

        decision_eg = z_eg>0;
        BER_EG = BER_EG + mean(abs(decision_eg-Data)); 

        waitbar(e/length(ES),f,'Simulating...');

    end

    BER_Rayleigh_E(e) = BER_Rayleigh/N2;
    BER_Rayleigh_SC(e) = BER_SC/N2;
    BER_Rayleigh_EG(e) = BER_EG/N2;

    Gamma = sqrt(Es/(Es+1));
    BER_Rayleigh_T(e) = 0;
    BER_Rayleigh_SC_T(e) = 0;
    for k = 0:(3-1)
        BER_Rayleigh_T(e) = BER_Rayleigh_T(e) + (((1/2)*(1-Gamma))^3)*(((Gamma+1)/2)^k)*nchoosek(3-1+k,k);
    end

    for k = 0:3
        BER_Rayleigh_SC_T(e) = BER_Rayleigh_SC_T(e) + (1/2)*((-1)^k)*nchoosek(3,k)*(1+k/Es)^(-0.5);
    end

end

close(f)
semilogy(10*log10(ES),BER_Rayleigh_E,'-x',10*log10(ES),BER_Rayleigh_T,'--o',10*log10(ES),BER_Rayleigh_SC,'-^',10*log10(ES),BER_Rayleigh_SC_T,'--^',10*log10(ES),BER_Rayleigh_EG,'-<');
legend('3-MRC E','3-MRC T','3-SC E','3-SC T','3-EG E','Location','southwest')
title('3-branch MRC, Selection Combining and EGC - BER')
xlabel("SNR (dB)")
ylabel("BER")
grid on

% Q function implementation
function out=Q(x)
out=0.5.*erfc(x/sqrt(2));
end