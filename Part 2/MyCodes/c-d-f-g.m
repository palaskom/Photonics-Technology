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
LSB = 20;   % [um] z-length of S-Bendss
d = 0.42;  % [um] x-gap of tap-coupler waveguides

% Finite-Difference (FD) discretization of x-axis
Nx = 2000; % [.] number of x-samples -- Affects FD-BPM accuracy!

% PML (Perfectly Match Layers) for reflectionless absorption -- Be careful!
PML_Thk = 1*[ 1 1 ] ; % [um] PML Thickness for up/down layer, in x-dim
PML_Str = 1*[ 1 1 ] ; % [.] PML Absorption "strength"
PML_params = [ PML_Thk , PML_Str ];

%----------------- wavelength -------------------%
% wl_window = linspace(1.5, 1.6, 300);
wl_window = 1.55; % [um] operating wavelength
ER = zeros(1,length(wl_window));
%------------------------------------------------%

% Theory/Mode Analysis - Estimate coupling-length ( used in (c) )
neffs = MLSWG( 'TE' , 1.55 , nsg(1)*[1 1] , nsg([2 1 2]) , [w d w] );
Lc = 0.5*1.55 / abs(diff(real(neffs))); % [um] coupling length
kappa = pi/2/Lc; % [rad/m] coupling coefficient


%----------------- LDC --------------------%
% ldc = 43:0.05:49; % comment the plots!!
% ldc = 143:0.1:151; % comment the plots!!
ldc = 44; %144.5;
dP_Lc = zeros(1,length(ldc));
%------------------------------------------%


%----------------- dn --------------------%
% used in (d): LDC = 45 [um]
% dn = 0:0.001:0.1; % comment the plots!!

% used in (f)
% dn = 0.0125:0.0001:0.0175; % comment the plots!! 

dn = 0; %0.0147;
dT = dn*1e4/2;
dP_dn = zeros(1,length(dn));
%-----------------------------------------%


for m = 1:length(wl_window)

for i = 1:length(ldc)
    
    for j = 1:length(dn)
    
        LDC = ldc(i);
        wl = wl_window(m);

        % =========================================================================
        % Structure Layout (SL) -- TopView
        % =========================================================================

        % 1st Module: Input S-Bends waveguide
        XBR1 = [ -DSB/2 -(d+w)/2 w w 0 2 ]; 
        XBR2 = [ DSB/2 (d+w)/2 w w 0 2 ];
        SL1 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };

        % 2nd Module: Coupling region
        XBR1 = [ -(d+w)/2 -(d+w)/2 w w 0 2 ];
        XBR2 = [ (d+w)/2 (d+w)/2 w w  dn(j) 2 ];
        SL2 = { nsg , LDC , 2*round(LDC/wl) , XBR1 , XBR2 };

        % 3rd Module: Bend "tapped" waveguide away (decoupling)
        XBR1 = [ -(d+w)/2 -DSB/2 w w 0 2 ];
        XBR2 = [ (d+w)/2 DSB/2 w w 0 2 ];
        SL3 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };

        % Lets put everything together, in the ROWS of SL
        SL = [ SL1 ; SL2 ; SL3 ];

        % =========================================================================
        % Finite-Difference (FD) Discretization / TopView
        % =========================================================================

        % Discretize structure's cross-section plane (x-axis)
        x = linspace( -4 , 4 , Nx ) ; % [um] x-axis vector

        % Get index profile and z-axis discretized vector
        [ ns_zxProf, z ]=BPMFD2D_PreProcLayout( SL,x,PML_params );

        % Draw the FD-discretized SL
        figure;
        pcolor( z , x , ns_zxProf ); shading flat; colorbar;

        % Draw (on same axes) the side-walls of all waveguides in the SL
        LinCol = 'w' ; % Color of lines
        BPMFD2D_DrawLayout( gca, SL, x, PML_params, LinCol ); %axis equal tight;

        % Set axes cosmetics
        title( 'Structure TopView (refractive index)' );
        xlabel( 'z-axis (um)' ); 
        ylabel( 'x-axis (um)' );

        % =========================================================================
        % Generate BPM excitation (on input plane)
        % =========================================================================

        % input only in waveguide 1 (right)
        [ neffs , modeProfs ] = MLSWG( 'TE' , wl , nsg(1)*[1 1] , nsg(2) , w , x+(DSB+w)/2 );

        % In case the structure supported more modes:
        nref = neffs(1);
        Excitation = modeProfs(1,:);

        % =========================================================================
        % Beam Propagation Method (BPM)
        % =========================================================================

        % Call a single-routine to do the propagation, from the input of first
        % module of SL to the output of its last module. This returns the 2D matrix
        % "FA" that contains the field values in ALL the nodes of the FD-mesh.
        FA = BPMFD2D_DoProp( SL, x , PML_params , Excitation , nref , wl  );

        % =========================================================================
        % Post-processing
        % =========================================================================

        % Post-process
        FA = FA / max(abs( FA(:) )); % normalize abs-value to 1

        % Measure the device's Insertion Losses (IL), from simulation
        Ey_inp = FA(:,1); % [V/m] Ey field at input plane
        Ey_out = FA(:,end); % [V/m] Ey field at output plane
        dx = mean(diff(x)); % [um] x-step
        PT_inp = sum( abs( Ey_inp ).^2 )*dx;
        PT_out = sum( abs( Ey_out ).^2 )*dx;
%         fprintf( ' == Insertion Losses : %+4.2f dB \n' , 10*log10( PT_out/PT_inp ) );

        % Measure the tapped power ratio from simulation
        P1 = sum( abs( Ey_out( x < 0 ) ).^2 )*dx; 
        P2 = sum( abs( Ey_out( x > 0 ) ).^2 )*dx; 
%         fprintf( ' == Measured Tap : %+4.2f dB \n' , 10*log10( P2/(P1+P2) ) );
        % figure, plot(x,Ey_out)
        dP_Lc(i) = 10*log10(P2/P1); % used in (c)
        dP_dn(j) = 10*log10(P1/P2); % used in (d)
        
    end
    
%     plot(dn,dP_dn);
%     hold on
    
end

ER(m) = 10*log10(P1/P2); % used in (g)

end 

% =========================================================================
% Plots 1
% =========================================================================


%-------------------- (c) & (d)------------------------%

% % used in (c)
% figure, plot(ldc,dP_Lc)
% title('dn = 0')
% xlabel('L_{DC} [um]')
% ylabel('Rower Ratio (P_2/P_1) [dB]')
% hold on, plot(ldc, linspace(20, 20, length(ldc)), 'r')

% % used in (d)
% figure, plot(dn,dP_dn)
% hold on, plot(dn, linspace(20, 20, length(dn)), 'r')
% title('L_{DC} = 144.5 [um]')
% xlabel('dn')
% ylabel('Rower Ratio (P_1/P_2) [dB]')

% % find the dn_min that entails P1/P2 > 20
% index = find(dP_dn>20);
% dn_min = dn(min(index));

%-------------------- (f) & (g)------------------------%

% % used in (f)
% plot(wl_window, ER)
% hold on, plot(wl_window, linspace(20, 20, length(wl_window)), 'r')
% xlabel('? [?m]'), xlim([wl_window(1) wl_window(end)])
% ylabel('ER [dB]')
% title(['L_{DC} = ', num2str(ldc), '   dn = ', num2str(dn), '  -  Bar'])

% % used in (f)
% ER = -ER;
% plot(wl_window, ER)
% hold on, plot(wl_window, linspace(20, 20, length(wl_window)), 'r')
% xlabel('? [um]'), xlim([wl_window(1) wl_window(end)])
% ylabel('ER [dB]')
% title(['L_{DC} = ', num2str(ldc), ' um   dn = ', num2str(dn), '  -  Cross'])

% % Bandwidth (f)
% index = find(ER==max(ER));
% Tlim = 20;
% flag = false;  
% k = index;  
%  % bandwidth -> Tlim-dB
% while (~flag)
%     k = k+1;
%     if (ER(k)<Tlim) 
%         flag = true;
%         wl_up = wl_window(k);
%     end
% end
% 
% flag = false;  
% k = index; 
% while (~flag)
%     k = k-1;
%     if (ER(k)<Tlim) 
%         flag = true;
%         wl_down = wl_window(k);
%     end
% end
% 
% bandwidth = (wl_up - wl_down)*1e3; %[nm]

% =========================================================================
% Plots 2
% =========================================================================

% TopView -- Intensity
figure;
pcolor(z , x, 10*log10( abs(FA).^2 ) ); shading flat;
BPMFD2D_DrawLayout( gca, SL, x, PML_params, 'w' );
colorbar; 
caxis( [-30 0] );
title( 'BPM calculated : |Ey| (dB)' );
xlabel( 'z-axis (um)' ); 
ylabel( 'x-axis (um)' );

% Input and output profiles (along x-axis)
figure;
plot( x , 10*log10( abs(Ey_inp).^2 ) , 'b' ); hold on;
plot( x , 10*log10( abs(Ey_out).^2 ) , 'r' );
plot( x , -10+0*x , 'k:' );
legend( 'input' , 'output' );
xlabel( 'x-axis (um)' ); 
ylabel( '|Ey|' );
set( gca , 'YLim' , [-30 0] );

% Array the figures in the computer screen
fmfp;