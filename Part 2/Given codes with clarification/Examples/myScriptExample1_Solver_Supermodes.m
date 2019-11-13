%%
clear all; close all; clc;
path( pathdef );

% This ExampleScript studies the supermodes of a SOI-slab WG coupler using
% two different solvers (MLSWG & FIDODIES)
addpath( [pwd, '\Solver'] ); % Required path to access solver m-files
addpath( [pwd, '\Misc'] ); % Some more misc functions

% =========================================================================
% Parameters
% =========================================================================
wl  = 1.55; % [um] wavelength
nLR = [ 1.45 1.45 ]; % [.] refr.index of L/R semi-inf slabs (substrate/cladding)
ns  = [ 3.20 1.45 3.20 ]; % [.] refr.index of intermediate layers
wid = 0.30; % [um] x-width of "core" layers
gap = 0.50; % [um] x-gap (substrate) separating the cores
ts  = [ wid gap wid ]; % [um] thickness of intermediate layers 

%% =========================================================================
% MLSWG - Characteristic equation solver // NewtonRaphson
% =========================================================================

% Define an x-space for displaying the mode profiles
xPlot = linspace( -1 , sum(ts)+1 , 400 ); % [um] cross-section x-axis

% Solve char-eq using Newton-Raphson method and calculate modes on xPlot
[ neXE , modeProfsXE ] = MLSWG( 'TE' , wl , nLR , ns , ts , xPlot );

% Plot all the mode Profiles 
figure;
linCols = hsv( length(neXE) );
myLegends = cell( 1,length(neXE) );
for ii = 1 : length(neXE)
    plot( xPlot , real(modeProfsXE(ii,:)) , 'Color' , linCols(ii,:) ); hold on;
    myLegends{ii} = sprintf( ' TE%d | neff = %6.4f' , ii , real(neXE(ii)) );
end
legend( myLegends );
set( gca, 'XLim' , xPlot([1 end]) ); % set axis limits
title( 'MLSWG (Newton-Raphson) modes');

%% =========================================================================
% FIDODIES - Numerical FD-based eigenmode solver // SPTARN
% =========================================================================

% Now, get the same eigenmodes from numerical FD solver (FIDODIES). Here,
% we need to define a DENSE cross-section, since the accuracy of the solver
% depends on it. The x-axis (xPlot) defined previously was only used for 
% plotting, i.e. the mode-profile was calculated from closed-form 
% expressions on that particular x-vector. Here, the x-axis vector (xFD) is
% part of the solver.
xFD = linspace( -4 , sum(ts)+4 , 5000 ); % [um] x-axis

% We need to populate the nsFD matrix, with the refractive indices at each
% point of the xFD vector
nsFD = NaN*xFD; % initialize
nsFD( xFD <  0 ) = nLR(1); % left semi-inf region
nsFD( xFD >= sum(ts) ) = nLR(2); % right semi-inf region
for jj = 1 : length(ts) % intermediate layer #jj
    nsFD( xFD >= sum(ts(1:jj-1)) & xFD < sum(ts(1:jj)) ) = ns(jj); 
end

% Ready to run the numerical FD solver FIDODIES
[ neFD , modeProfsFD ] = FIDODIESv2( wl , xFD , nsFD );

% Plot all the modes found
figure;
linCols = hsv( length(neFD) );
myLegends = cell( 1,length(neFD) );
for ii = 1 : length(neFD)
    aux = real(modeProfsFD(:,ii)); aux = aux / max(abs(aux)); %normalize
    plot( xFD , aux , 'Color' , linCols(ii,:) ); hold on;
    myLegends{ii} = sprintf( ' TE%d | neff = %6.4f' , ii , real(neFD(ii)) );
end
legend( myLegends );
set( gca, 'XLim' , xPlot([1 end]) ); % set axis limits
title( 'FIDODIES (SPTARN) modes');

%% =========================================================================
% ERROR CHECKING
% =========================================================================

% In some cases, the Newton-Raphson solver in MLSWG might NOT find all the
% roots of the Characteristic Equation! FIDODIES (numerical FD solver)
% always finds all modes, but is less precise. You should rerun MLSWG
if length(neFD)~=length(neXE)
   fprintf( ' ## MLSWG missed %d mode(s)! -- Rerun script :)\n' , ...
       length(neFD)-length(neXE) );
end

% Finally, array the figures in the computer screen using a custom function
fmfp;

