clear;

% Wavelength

lambda = 0.06;

% Number of elements of the IRS
L = 200;

% Pathloss exponent, distance, pathloss at a distance of 1m
Lo = 10^(-3);
dti = 5.17;
dir = 20;
drt = 25;
pl = 1.5;

% Rician channels initialized

HAA = get_H(9.82,9.42,17.35,lambda,L,Lo,pl);

HBB = get_H(15.24,14.95,17.56,lambda,L,Lo,pl);

H11 = get_H(9.82,15.24,7.96,lambda,L,Lo,pl);

 H21 = get_H(9.42,15.24,19.42,lambda,L,Lo,pl);

H12 = get_H(9.82,14.95,19.08,lambda,L,Lo,pl);

H22 = get_H(9.42,14.95,8.28,lambda,L,Lo,pl);

HA1C = get_H(9.82,8.94,9.90,lambda,L,Lo,pl);

HA2C = get_H(9.42,8.94,9.91,lambda,L,Lo,pl);

HB1C = get_H(15.24,8.94,9.57,lambda,L,Lo,pl);

HB2C = get_H(14.95,8.94,9.22,lambda,L,Lo,pl);

% HAA = get_H(10,17.32,10,lambda,L,Lo,pl);
% HBB = get_H(10,17.32,10,lambda,L,Lo,pl);
% 
% H11 = get_H(10,10,17.32,lambda,L,Lo,pl);
% H21 = get_H(17.32,10,20,lambda,L,Lo,pl);
% H12 = get_H(10,17.32,20,lambda,L,Lo,pl);
% H22 = get_H(17.32,17.32,17.32,lambda,L,Lo,pl);
% 
% HA1C = get_H(10,20,17.32,lambda,L,Lo,pl);
% HA2C = get_H(17.32,20,10,lambda,L,Lo,pl);
% HB1C = get_H(10,20,17.32,lambda,L,Lo,pl);
% HB2C = get_H(17.32,20,10,lambda,L,Lo,pl);

% IRS reflection matrix (Initialized randomly)
wi = (2*rand(L,1)-1)*pi;
w = [wi ; 1].';

% Initializing power values
PA1 = 0.001;
PA2 = 0.001;
PB1 = 0.001;
PB2 = 0.001;

PA1_opt = PA1;
PA2_opt = PA2;
PB1_opt = PB1;
PB2_opt = PB2;



i = 20;
while i>0
    i = i-1;

cvx_begin quiet
cvx_solver mosek
  variable w_opt(1,L+1) complex
  variable t
  maximize t
  subject to
%A1
  (get_P(PB1, w, H11, PA2, HAA, PB2, H12)-get_Q(PB1, w, HB1C) + ...
      2*real((sqrt(PB1)*w_opt*H11)*(sqrt(PB1)*w*H11))/get_A(PA2, w, HAA, PB2, H12) ...
      -(get_channel_term(PB1, w, H11)/(get_A(PA2, w, HAA, PB2, H12)* ...
      (get_A(PA2, w, HAA, PB2, H12)+get_channel_term(PB1, w, H11))))*(PA2*pow_abs(w_opt * HAA, 2) ...
      + PB2*pow_abs(w_opt * H12, 2)+10^(-14.5)+10^(-14)+PB1*pow_abs(w_opt * H11, 2))-(get_channel_term(PB1, w, H11)/get_A(PA2, w, HAA, PB2, H12))...
      -(PB1*pow_abs(w_opt * HB1C, 2)-get_channel_term(PB1, w, HB1C))/(10^(-14.5)+get_channel_term(PB1, w, HB1C))) >= t;

%A2
  (get_P(PB2, w, H22, PA1, HAA, PB1, H21)-get_Q(PB2, w, HB2C) + ...
      2*real((sqrt(PB2)*w_opt*H22)*(sqrt(PB2)*w*H22))/get_A(PA1, w, HAA, PB1, H21) ...
      -(get_channel_term(PB2, w, H22)/(get_A(PA1, w, HAA, PB1, H21)* ...
      (get_A(PA1, w, HAA, PB1, H21)+get_channel_term(PB2, w, H22))))*(PA1*pow_abs(w_opt * HAA, 2) ...
      + PB1*pow_abs(w_opt * H21, 2)+10^(-14.5)+10^(-14)+PB2*pow_abs(w_opt * H22, 2))-(get_channel_term(PB2, w, H22)/get_A(PA1, w, HAA, PB1, H21))...
      -(PB2*pow_abs(w_opt * HB2C, 2)-get_channel_term(PB2, w, HB2C))/(10^(-14.5)+get_channel_term(PB2, w, HB2C))) >= t;

%B1
  (get_P(PA1, w, H11, PB2, HBB, PA2, H21)-get_Q(PA1, w, HA1C) + ...
      2*real((sqrt(PA1)*w_opt*H11)*(sqrt(PA1)*w*H11))/get_A(PB2, w, HBB, PA2, H21) ...
      -(get_channel_term(PA1, w, H11)/(get_A(PB2, w, HBB, PA2, H21)* ...
      (get_A(PB2, w, HBB, PA2, H21)+get_channel_term(PA1, w, H11))))*(PB2*pow_abs(w_opt * HBB, 2) ...
      + PA2*pow_abs(w_opt * H21, 2)+10^(-14.5)+10^(-14)+PA1*pow_abs(w_opt * H11, 2))-(get_channel_term(PA1, w, H11)/get_A(PB2, w, HBB, PA2, H21))...
      -(PA1*pow_abs(w_opt * HA1C, 2)-get_channel_term(PA1, w, HA1C))/(10^(-14.5)+get_channel_term(PA1, w, HA1C))) >= t;

%B2
  (get_P(PA2, w, H22, PB1, HBB, PA1, H12)-get_Q(PA2, w, HA2C) + ...
      2*real((sqrt(PA2)*w_opt*H22)*(sqrt(PA2)*w*H22))/get_A(PB1, w, HBB, PA1, H12) ...
      -(get_channel_term(PA2, w, H22)/(get_A(PB1, w, HBB, PA1, H12)* ...
      (get_A(PB1, w, HBB, PA1, H12)+get_channel_term(PA2, w, H22))))*(PB1*pow_abs(w_opt * HBB, 2) ...
      + PA1*pow_abs(w_opt * H12, 2)+10^(-14.5)+10^(-14)+PA2*pow_abs(w_opt * H22, 2))-(get_channel_term(PA2, w, H22)/get_A(PB1, w, HBB, PA1, H12))...
      -(PA2*pow_abs(w_opt * HA2C, 2)-get_channel_term(PA2, w, HA2C))/(10^(-14.5)+get_channel_term(PA2, w, HA2C))) >= t;

   
    for k = 1 : L
        abs(w_opt(1,k)) <= 1;
    end
  w_opt(1,L+1)== 1;
  t <= 1000;


cvx_end
t
w = w_opt;

cvx_begin quiet
cvx_solver mosek
  variable PA1_opt
  variable PA2_opt
  variable PB1_opt
  variable PB2_opt
  variable t
  maximize t

  subject to 
%A1
  (get_P(PB1, w, H11, PA2, HAA, PB2, H12)-get_Q(PB1, w, HB1C) + ...
      2*PB1_opt*real((w*H11)*(sqrt(PB1)*w*H11))/get_A(PA2, w, HAA, PB2, H12) ...
      -(get_channel_term(PB1, w, H11)/(get_A(PA2, w, HAA, PB2, H12)* ...
      (get_A(PA2, w, HAA, PB2, H12)+get_channel_term(PB1, w, H11))))*(PA2*pow_abs(w * HAA, 2) ...
      + PB2*pow_abs(w * H12, 2)+10^(-14.5)+10^(-14)+square(PB1_opt)*pow_abs(w * H11, 2))-(get_channel_term(PB1, w, H11)/get_A(PA2, w, HAA, PB2, H12))...
      -(square(PB1_opt)*pow_abs(w * HB1C, 2)-get_channel_term(PB1, w, HB1C))/(10^(-14.5)+get_channel_term(PB1, w, HB1C))) >= t;

%A2
  (get_P(PB2, w, H22, PA1, HAA, PB1, H21)-get_Q(PB2, w, HB2C) + ...
      2*PB2_opt*real((w*H22)*(sqrt(PB2)*w*H22))/get_A(PA1, w, HAA, PB1, H21) ...
      -(get_channel_term(PB2, w, H22)/(get_A(PA1, w, HAA, PB1, H21)* ...
      (get_A(PA1, w, HAA, PB1, H21)+get_channel_term(PB2, w, H22))))*(PA1*pow_abs(w * HAA, 2) ...
      + PB1*pow_abs(w * H21, 2)+10^(-14.5)+10^(-14)+square(PB2_opt)*pow_abs(w * H22, 2))-(get_channel_term(PB2, w, H22)/get_A(PA1, w, HAA, PB1, H21))...
      -(square(PB2_opt)*pow_abs(w * HB2C, 2)-get_channel_term(PB2, w, HB2C))/(10^(-14.5)+get_channel_term(PB2, w, HB2C))) >= t;

%B1
  (get_P(PA1, w, H11, PB2, HBB, PA2, H21)-get_Q(PA1, w, HA1C) + ...
      2*PA1_opt*real((w*H11)*(sqrt(PA1)*w*H11))/get_A(PB2, w, HBB, PA2, H21) ...
      -(get_channel_term(PA1, w, H11)/(get_A(PB2, w, HBB, PA2, H21)* ...
      (get_A(PB2, w, HBB, PA2, H21)+get_channel_term(PA1, w, H11))))*(PB2*pow_abs(w * HBB, 2) ...
      + PA2*pow_abs(w * H21, 2)+10^(-14.5)+10^(-14)+square(PA1_opt)*pow_abs(w * H11, 2))-(get_channel_term(PA1, w, H11)/get_A(PB2, w, HBB, PA2, H21))...
      -(square(PA1_opt)*pow_abs(w * HA1C, 2)-get_channel_term(PA1, w, HA1C))/(10^(-14.5)+get_channel_term(PA1, w, HA1C))) >= t;

%B2
  (get_P(PA2, w, H22, PB1, HBB, PA1, H12)-get_Q(PA2, w, HA2C) + ...
      2*PA2_opt*real((w*H22)*(sqrt(PA2)*w*H22))/get_A(PB1, w, HBB, PA1, H12) ...
      -(get_channel_term(PA2, w, H22)/(get_A(PB1, w, HBB, PA1, H12)* ...
      (get_A(PB1, w, HBB, PA1, H12)+get_channel_term(PA2, w, H22))))*(PB1*pow_abs(w * HBB, 2) ...
      + PA1*pow_abs(w * H12, 2)+10^(-14.5)+10^(-14)+square(PA2_opt)*pow_abs(w * H22, 2))-(get_channel_term(PA2, w, H22)/get_A(PB1, w, HBB, PA1, H12))...
      -(square(PA2_opt)*pow_abs(w * HA2C, 2)-get_channel_term(PA2, w, HA2C))/(10^(-14.5)+get_channel_term(PA2, w, HA2C))) >= t;

  0 <= PA1_opt;square(PA1_opt) <= 0.04;
  0 <= PA2_opt;square(PA2_opt) <= 0.04;
  0 <= PB1_opt;square(PB1_opt) <= 0.04;
  0 <= PB2_opt;square(PB2_opt) <= 0.04; 
  t <= 1000;


cvx_end

t

PA1 = PA1_opt^2
PA2 = PA2_opt^2
PB1 = PB1_opt^2
PB2 = PB2_opt^2

RA1 = get_secrecy_rate(w, H11, HAA, H12, HB1C, PB1, PA2, PB2, PB1)
RA2 = get_secrecy_rate(w, H22, HAA, H21, HB2C, PB2, PA1, PB1, PB2)
RB1 = get_secrecy_rate(w, H11, HBB, H21, HA1C, PA1, PB2, PA2, PA1)
RB2 = get_secrecy_rate(w, H22, HBB, H12, HA2C, PA2, PB1, PA1, PA2)

end





% function c_rate = get_secrecy_rate(w, eves, peves)
function c_rate = get_secrecy_rate(w, leg, inf1, inf2, eves, pleg, pinf1, pinf2, peves)

    sigma_c = 10^(-14.5);
    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    eves_inf = log2(1 + (get_channel_term(peves, w, eves))/(sigma_c));
    leg_inf = log2(1 + (get_channel_term(pleg, w, leg))/(get_channel_term(pinf1, w, inf1)+get_channel_term(pinf2, w, inf2)+sigma_ab+sigma_loop));
    c_rate = leg_inf - eves_inf;

end


function channel_term = get_channel_term(power, w, H)

    channel_term = power*(abs(w*H)^2);
end


function H = get_H(dti,dir,dtr,wave_l,L,Lo,pl)

    gti = get_rice_channel(dti, wave_l,L);
    hti = sqrt(Lo*dti^(-pl))*gti;
    gir = get_rice_channel(dir, wave_l,L);
    hir = sqrt(Lo*dir^(-pl))*gir;
    gtr = get_rice_direct_channel(dtr, wave_l);
    htr = sqrt(Lo*dtr^(-pl))*gtr;
    H = cat(1,hti.*hir, htr);
end


function rice_channel = get_rice_channel(dist, wave_l, L)
    rice_channel = sqrt(1/32)*(randn(L,1)+1i*randn(1)) + sqrt(30/32)*exp(1i*2*pi*(dist/wave_l-floor(dist/wave_l)));
end


function rice_direct_channel = get_rice_direct_channel(dist, wave_l)
    rice_direct_channel = sqrt(1/32)*(randn(1)+1i*randn(1)) + sqrt(30/32)*exp(1i*2*pi*(dist/wave_l-floor(dist/wave_l)));
end


function Q = get_Q(p, w, eves)

    sigma_c = 10^(-14.5);
    Q = log2(1 + (get_channel_term(p, w, eves))/sigma_c);
end


function P = get_P(pleg, w, leg, pinf1, inf1, pinf2, inf2)

    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    P = log2(1 + (get_channel_term(pleg, w, leg))/(get_channel_term(pinf1, w, inf1)+get_channel_term(pinf2, w, inf2)+sigma_ab+sigma_loop));
end


function A = get_A(pinf1, w, inf1, pinf2, inf2)

    sigma_ab = 10^(-14.5);
    sigma_loop = 10^(-14);
    A = get_channel_term(pinf1, w, inf1)+get_channel_term(pinf2, w, inf2)+sigma_ab+sigma_loop;
end
