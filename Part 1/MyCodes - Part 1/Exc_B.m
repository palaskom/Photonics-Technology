%% (iv), (v)
clear;clc
wl = (1450:0.1:1700)*1e-9; % [1450 1700] nm
nglass = 1.5;

% data - Bragg
wlo = 1.55e-6; % 1550nm
M = 15;
nH = 2.32;
nL = 1.45;
% refractive index of GaAs
% nH = 3.37; 

% data - Liquid Crystal (LC)
d = 1.9e-6;  N = 20;

V = 1.5;%1.5:1.5:4.5;
wl_res = zeros(1,length(V)); % wavelength of resonance
bandwidth = zeros(1,length(V));

for i = 1:length(V)

    % First Bragg (-> LC): na = neff_N, nb = nglass
    na = nLC(V(i),d,d);  nb = nglass;
    r_start = ((nH-nb)/(nH+nb))*ones(1,length(wl));
    rB_1 = Bragg(wl, wlo, na, nH, nL, M, r_start);

    % LC: na = nH, 
    %     nb: use the previous total reflection coefficient -> rBragg
    na = nH;
    [rLC,neff] = LC(wl, na, d, V(i), N, rB_1);
    % neff -> layers

    % Second Bragg: 
    % na = nglass, 
    % nb: use the previous total reflection coefficient -> rLC
    na = nglass;
    rB_2 = Bragg(wl, wlo, na, nH, nL, M, rLC);

    % TOTAL
    rtot = rB_2;
    T_FPLC = 1-abs(rtot).^2;
    TdB = 10*log10(T_FPLC);
    hold on, plot(wl*1e6,TdB)

    % wavelength - resonance
    index = find(TdB==max(TdB));
    wl_res(i) = wl(index)*1e6; % [um]
    
    % Bandwidth
    f = false;  
    j = index;  
    Tlim = -3; % bandwidth -> Tlim-dB
    while (~f)
        j = j+1;
        if (TdB(j)<Tlim) 
            f = true;
            wl_up = wl(j);
        end
    end

    f = false;  
    j = index; 
    while (~f)
        j = j-1;
        if (TdB(j)<Tlim) 
            f = true;
            wl_down = wl(j);
        end
    end

    bandwidth(i) = (wl_up - wl_down)*1e9; %[nm]

end

xlabel('{\lambda} [?m]')
ylabel('Power Transmission Coefficient (dB)')
legend('1.5V','3V','4.5V')
grid on

% title('V=4.5')
% legend('N=10','N=20','N=30')
% ylim([-3 0])

%% (vi) Voltage Approximation @ 1530, 1560 nm
p = polyfit(wl_res,V,2);

% wavelength = 1.53 um
v1 = polyval(p,1.53);

% wavelength = 1.56 um
v2 = polyval(p,1.56);

