%% data
clear;clc
wl = (1000:2000)*1e-9; % [1000 2000] nm
wlo = 1.55e-6; % 1550nm

no = 1.5;  ne = 1.7; 
nav = (no+ne)/2;
nglass = 1.5;
nH = 2.32;  nL = 1.45;

M = 9;

%% (i)
na = nav; nb = nglass;
r_start = ((nH-nb)/(nH+nb))*ones(1,length(wl));
rB = Bragg(wl, wlo, na, nH, nL, M, r_start);

RB = abs(rB).^2;
RB_dB = 10*log10(RB);
plot(wl*1e6,RB_dB), grid on
xlabel('{\lambda} [?m]')
ylabel('Power Reflection Coefficient (dB)')

%% (ii)

% glass -> LC: r1
na = nglass; nb = nav;
r_start = ((nH-nb)/(nH+nb))*ones(1,length(wl));
r1 = Bragg(wl, wlo, na, nH, nL, M, r_start);

% LC -> glass: r2
na = nav; nb = nglass;
r_start = ((nH-nb)/(nH+nb))*ones(1,length(wl));
r2 = Bragg(wl, wlo, na, nH, nL, M, r_start);

rmult = r1.*r2;
delta = 8*pi*(wlo./wl); %4*pi*d*nav/wl;
T1 = ((1-abs(r1).^2).^2)./abs(1-rmult.*exp(-1i*delta)).^2;

T1_dB = 10*log10(T1);
figure, plot(wl*1e6,T1_dB)
xlabel('{\lambda} [?m]')
ylabel('Power Transmission Coefficient (dB)')

%% (iii)

% the fisrt Bragg (na = nav, nb = nglass)
r_start = ((nH-nglass)/(nH+nglass))*ones(1,length(wl));
r_total = Bragg(wl, wlo, nav, nH, nL, M, r_start);
% the Liquid Crystal
r_mid = (nH-nav)/(nH+nav);
b_mid = 4*pi*(wlo./wl);
r_total = (r_mid+r_total.*exp(-2i*b_mid))./(1+r_mid*r_total.*exp(-2i*b_mid));
% the second Bragg mirror (na = nglass, nb: we use the initial value r_(M+1)=r_total)
r_total = Bragg(wl, wlo, nglass, nH, nL, M, r_total);

T2 = 1-abs(r_total).^2;
T2_dB = 10*log10(T2);

hold on
plot(wl*1e6,T2_dB,'*')
legend('Fabry-Perot','Multilayer')
xlim([1 2])