%% Parameters 
close all; clear;clc

wl = 1.55;

ns = 1.45;   ng = 3.2;
nLR = [ns ns];
nIL  = [ng ns ng]; % refr.index of intermediate layers

w = 0.25; % [um] x-width of "core" layers

% Discretize structure's cross-section plane (x-axis)
Nx = 5000;  % number of samples on x-axis
x = linspace( -3 , 3 , Nx ) ; % [um] x-axis vector
dx = mean(diff(x)); % [um] x-step


%% Lc-d diagram

d = 0.8 : -(0.8-0.02)/500 : 0.02;
% d = 0.42;
Lc_supmod = zeros(1,length(d)); % from (a)
Lc_theoretical = zeros(1,length(d)); % Lc = pi/(2k)

for i = 1:length(d)
    
    x_offset = (w+d(i)/2); % Offset required to focus beam on bottom waveguide

    N_iterations = 0;
    while N_iterations < 10
        [ neffsSM , EySM ] = MLSWG( 'TE' , wl , nLR , nIL , [w d(i) w] , x+x_offset );
        if length(neffsSM) > 1
            Lc_supmod(i) = 0.5*wl / abs(real(neffsSM(1))-real(neffsSM(2))); % [um] coupling length
            break;
        else
            N_iterations = N_iterations + 1;
            if N_iterations == 10
                disp( ' ## ERROR: MLSWG failed to find 2 supermodes.' );
                return;
            end        
        end
    end  
    
    % lone waveguide 1
    [neffsWG1 , EysWG1] = MLSWG( 'TE' , wl , nLR , ng , w , x+x_offset );
    % lone waveguide 2
    [neffsWG2 , EysWG2] = MLSWG( 'TE' , wl , nLR , ng , w , x-d(i)/2);
    % In case the structure supported more modes:
    neffWG = neffsWG1(1);
    EyWG_1 = EysWG1(1,:);
    EyWG_2 = EysWG2(1,:);
%     plot(x, EyWG_1)
%     hold on, plot(x, EyWG_2)
    
    P11 = sum( abs(EyWG_1).^2 )*dx ;    
    P21 = sum( EyWG_2(x>=-w-d(i)/2 & x<=-d(i)/2).* conj(EyWG_1(x>=-w-d(i)/2 & x<=-d(i)/2)) )*dx ;
     
    Lc_theoretical(i) = (P11/P21)*wl*neffWG/2/(ng^2-ns^2);

end

% plots
figure, plot(d, Lc_theoretical*1e-3, 'b')
hold on, plot(d(1:10:end), Lc_supmod(1:10:end)*1e-3, '*r')

grid on, title('w = 250nm')
xlabel('d [um]'), ylabel('L_c [mm]')
legend('CMT', 'SM')



