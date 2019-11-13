%% Parameters

clc; clear; close all; 

addpath( [pwd, '\BPM'] ); % Required path to access BPM-related m-files
addpath( [pwd, '\Solver'] ); % Required path to access solver m-files
addpath( [pwd, '\Misc'] ); % Some more misc functions

% =========================================================================
% Structure Parameters
% =========================================================================

% Fixed parameters
nsg = [1.45 3.2]; % [] refr.indices of substrate & guiding-layer 
w = 0.25; % [um] x-width of waveguide (input plane)
DSB = 4;    % [um] x-position of tap's output ports (input waveguide is at x=-dSB/2);
d = 0.42;  % [um] x-gap of tap-coupler waveguides

% Finite-Difference (FD) discretization of x-axis
Nx = 2000; % [.] number of x-samples -- Affects FD-BPM accuracy!
% Discretize structure's cross-section plane (x-axis)
x = linspace( -2 , 2 , Nx ) ; % [um] x-axis vector

%----------------- wavelength -------------------%
wl = 1.55; % [um] operating wavelength
%------------------------------------------------%

% Theory/Mode Analysis - Estimate coupling-length 
neffs_sm = MLSWG( 'TE' , 1.55 , nsg(1)*[1 1] , nsg([2 1 2]) , [w d w] );
Lc = 0.5*wl / abs(diff(real(neffs_sm))); % [um] coupling length
kappa = pi/2/Lc; % [rad/m] coupling coefficient

%% symmetric & antisymmatric modes for dn>0

%------------------ dn -------------------%
dn = 0.027; 
%-----------------------------------------%

[neffs_sm, Ey_dn] = MLSWG( 'TE' , wl , nsg(1)*[1 1] , [nsg(2) nsg(1) nsg(2)+dn] , [w d w],...
                           x+w+d/2 );
                     
% Ey_dn = Ey_dn / max(abs( Ey_dn(:) )); % normalize abs-value to 1

figure, plot(x, Ey_dn(1,:), 'b')
hold on, plot(x, Ey_dn(2,:), 'r')

%% (e)

%------------------ dn -------------------%
dn = 0:0.1/500:0.1 ;
%-----------------------------------------%
 
neff_1 = MLSWG( 'TE' , wl ,nsg(1)*[1 1] , nsg(2) , w );

Lb_cm = zeros(1,length(dn)); % [um] beating length -> cm: theory of coupled modes
Lb_sm = zeros(1,length(dn)); % [um] beating length -> sm: supermodes

diff_Lcm = 100*ones(1,length(dn));
diff_Lsm = 100*ones(1,length(dn));

for m = 1:length(dn)
    
    % SUPERMODES (SM)
    N_iterations = 0;
    while N_iterations < 5
        neff_2 = MLSWG( 'TE' , wl ,nsg(1)*[1 1] , nsg(2)+dn(m) , w );
        if length(neff_2) == 1
            Lb_cm(m) = wl / abs(neff_2-neff_1);
            break;
        else
            N_iterations = N_iterations + 1;
            if N_iterations == 5
                disp( ' ## ERROR: MLSWG failed to find 1 mode.' );
                return;
            end  
        end   
    end
    
    % THEORY OF COUPLED MODES (CMT)
    N_iterations = 0;
    while N_iterations < 10
        neffs_sm = MLSWG( 'TE' , wl , nsg(1)*[1 1] , [nsg(2) nsg(1) nsg(2)+dn(m)] , [w d w] );
        if length(neffs_sm) == 2
            Lb_sm(m) = wl / abs(diff(real(neffs_sm))); 
            break;       
        else % there must be 2 modes
            N_iterations = N_iterations + 1;
            if N_iterations == 10
                disp( ' ## ERROR: MLSWG failed to find 2 supermodes.' );
                return;
            end  
        end
    end
    
    diff_Lcm(m) = abs(Lc-0.5*sqrt(3)*Lb_cm(m));
    diff_Lsm(m) = abs(Lc-Lb_sm(m));
    if diff_Lcm(m) < 0.2
        dn_cm = dn(m); % dn -> supermodes
        Lc_cm = 0.5*sqrt(3)*Lb_cm(m);
    end
    if diff_Lsm(m) < 0.1
        dn_sm = dn(m); % dn -> supermodes
        Lc_sm = Lb_sm(m);
    end
    
end


min_Lcm = min(diff_Lcm);
min_Lsm = min(diff_Lsm);

% supermodes: Lb -> Lc
figure, plot(dn, diff_Lsm)

% theory of coupled modes: Lb -> 0.5*sqrt(3)*Lc
hold on, plot(dn(dn>0.007), diff_Lcm(dn>0.007))

xlabel('dn')
ylabel("|L_c - L'_b| [um]")
legend("SM: L'_b \equiv L_b \rightarrow  L_c", "CMT: L'_b \rightarrow  2L_c/\surd 3")
