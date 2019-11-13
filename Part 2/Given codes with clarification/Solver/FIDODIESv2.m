function [neff,V] = FIDODIESv2( wl , x , ns , sx )

% FIDODIES: FInite Difference One DImensional Eigenmode Solver / TE modes
%
% You provide a wavelength ("wl"), the 1D-space vector of the waveguide
% ("x") which should be equispaced (e.g. from linspace). At each of the
% nodes defined at "x", you assign a complex-valued refr.indx ("ns"
% vector) and a PML ("sx" vector, sx=1-1j*PML_Strength), for the nodes near
% the two edges of the 1D-window.
%
% Outputs is a list of the effective indices ("neff") of the eigenmodes as
% well as the corresponding field-profiles ("V") of the eigenmodes.
%
% Notes:
% ------
% ** Input "sx" is optional. It will be set to default if not inputted.
% ** This currently works only for TE polarization. TM can be unlocked (see
%    in-function parameters below), but it doesn't yield correct results.
% ** The denser the x-matrix (more points), the higher the FDM accuracy
% ** Set to use the SPTARN solver, which find ALL the eigenvalues within a
%    specific interval. This interval is appropriately chosen, assuming
%    regular photonic modes, that do not "oscillate" in the substrates.
% ** EIGS solver might be faster, but it only finds a given number of modes
%    near a specific value. We do not generally know how many modes are
%    supported, so it's not a good practice in this case.
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2011 July : Original version 
%  2015 Sept : Modified version

% Test params
if nargin == 0
    
    close all; clc;
    
    % Input params
    wl = 1.55; % wavelength
    h  = 0.25; % thickness
    g  = 0.42; % gap, (g=[] to disable)    
    ngui = 3.20; % gui-index
    nsub = 1.45; % sub-index
    ngap = 1.45; % gap/cladding-index
    
    % x-space grid
    x = 8*h * linspace(-1,+1,5001);    
    
    % PMLs tensors (sx)
    sx = []; % to use default values (see below)
    
    % Index-profile
    ns = nsub * ones(size(x)); % initialize -- all at n_substrate
    if isempty( g ) % Single-WG (step or graded profile)
        
        % Graded-index profile
        Delta = (ngui^2 - nsub^2) / ( 2 * ngui^2) ;
        gpr   = Inf; % Profile {Inf,1,2,...}->{Step,Triangle,Quadratic,...} 
        ii = abs(x)<h/2 ; 
        ns( ii )  = ngui * sqrt( 1 - 2*Delta*(abs(x(ii))/h*2).^gpr );
        ns( abs(x)>=h/2 ) = ngui * sqrt( 1 - 2*Delta );
        
    else % Coupler (step-index)
        
        ns( abs(x-(g+h)/2)<=h/2 ) = ngui; % Left 
        ns( abs(x+(g+h)/2)<=h/2 ) = ngui; % Right     
        ns( abs(x)<g/2 ) = ngap; % Gap
        
        ns( x>0 & abs(x-(g+h)/2)<=h/2 ) = ngui - 1*0.025; % de-sync branches
        
    end

end

% Default PML
if nargin < 4, sx = []; end % to use default values (see below)

%-------------------------------------------------------------------------
% More parameters
%-------------------------------------------------------------------------

% Method Related
ModePol     = 'TE' ; % {'TE','TM'} only TE works fine
Formulation = 'std' ; % {'std','hyb','alt'} use standard (std)
NumSolv     = 'sptarn'; % {'sptarn','eigs'} use sptarn
NEV         = 2; % Number of eigenvalues (if EIGS solver is used)

% Index-profile smoothing-window (odd). Set [] for no-smoothing
nsf = [];

% Tolerance in Delta-refr.indx (identifying Material Interfaces for plots)
dnTol = 0.05;

% PMLs tensors (sx) -- Default values
if isempty( sx )
    sx = ones(size(x)); % initialize -- no anisotropy/absorption
    sPML = 1; % strength of PML
    wPML = (max(x)-min(x))/10; % thickness of PML
    if sPML ~= 0
        iL = x<=(x( 1 )+wPML); % indices of Left PML
        iR = x>=(x(end)-wPML); % indices of Right PML
        sx( iL ) = sx( iL ) - 1j*sPML*( ( sum(iL==1):-1:1 )/sum(iL==1) ).^2 ; % Parabolic profile
        sx( iR ) = sx( iR ) - 1j*sPML*( ( 1:+1:sum(iR==1) )/sum(iR==1) ).^2 ; % Parabolic profile
    end
end


%-------------------------------------------------------------------------
% Misc operations
%-------------------------------------------------------------------------

% Locate material interfaces
xi = 0.5*( x(1:end-1) + x(2:end) ); % x-locs of interfaces between nodes
xInterfaces = xi( abs(diff(ns)) > dnTol ); % x-locs of material interfaces

% Constants
mu0 = 4*pi*1e-13; % [H/um] free-space magnetic permeability
ep0 = 8.854187817e-18; % [F/um] free-space electric permittivity
c0  =  1/sqrt(mu0*ep0); % [um/sec] free-space speed-of-light

% Define aux params
N = length(x); % [.] Grid-Points in x-vector
k0 = 2*pi/wl; % [1/um] free-space wave-number
wm = k0*c0; % [1/sec] omega
dx = mean(diff(x)); % x-step

% Smooth the refr-index profile (ns):
if ~isempty( nsf )
    fprintf(' ** Smoothing Refractice Index Profile n(x)\n' );
    if mod(length(nsf),2)~=1 , error(' ## Refr.Indx. Smoothing Filter should be of odd-length'); end
    ns = filter2( nsf/sum(nsf) , ns , 'same' );
    for k = 1 : (length(nsf)-1)/2
        ns(k)      =ns((length(nsf)-1)/2 + 1);
        ns(end-k+1)=ns(end-(length(nsf)-1)/2-1);
    end
end

% Error-check
if any( real(ns) < 1 )
    disp( ' ## Warning: Some ns<1 ! ' );
end

%-------------------------------------------------------------------------
% Construct Sparse-Equation System
%-------------------------------------------------------------------------

% Aux-matrices
sx_n = [ sx(2:end) , sx(end  ) ]; % PML sx / next-node
sx_p = [ sx( 1 ) , sx(1:end-1) ]; % PML sx / prev-node
es   = ns.^2; % dielectric-constant profile
es_n = [ es(2:end) , es(end  ) ]; % es / next-node
es_p = [ es( 1 ) , es(1:end-1) ]; % es / prev-node

% Formulation Dependent
if strcmp( ModePol , 'TE' )

    if strcmp( Formulation , 'std' ) %"Standard" working approach
        
        %Equation: d^2(Ey)/dx^2 + ( k0*n(x) )^2*Ey = beta^2*Ey
        % * PMLs: 1/s*diff( 1/s * diff(Ey) )/dx^2 + ...

        %FDM application: 1/u -> 2/(u_next+n_prev)
        T_next =     2 ./ sx_n ./ ( sx_n + sx );  % "Trans" coeff, "next-node"
        T_prev =     2 ./ sx_p ./ ( sx_p + sx );  % "Trans" coeff, "prev-node"
        R_next = 1 + 2 ./ sx ./ ( sx_n + sx ) ;  % "Refl" coeff, "next-node"
        R_prev = 1 + 2 ./ sx ./ ( sx_p + sx ) ;  % "Relf" coeff, "prev-node"

        %Sparse-Diagonals
        diag00 = (ns*k0).^2+(2-R_next-R_prev)/dx^2;
        diag1p = T_next/dx^2;
        diag1m = T_prev/dx^2;
        A = spdiags( [ diag1m(:) , diag00(:) , diag1p(:) ] , [-1 0 +1] , N , N );

    else %{Ey;Hx} dual-transverse component approach (works well)
        %U/R diagonals:
        diag1p = 2 ./ ( sx_n + sx ) /dx^2 /-wm/mu0 ; % "Trans" coeff, "next-node"
        diag1m = 2 ./ ( sx_p + sx ) /dx^2 /-wm/mu0; % "Trans" coeff, "prev-node"
        auxi = -diag1p-diag1m ;
        diag00 = -sx.*ns.^2.*wm*ep0 + auxi;

        %D/L diagonal:
        diag0m = -1./sx*wm*mu0;

        %Prepare Diagonals for sparse matrix
        ds(:,1) = [ zeros(N,1) ; diag1m(:) ];
        ds(:,2) = [ zeros(N,1) ; diag00(:) ];
        ds(:,3) = [ zeros(N,1) ; diag1p(:) ];
        ds(:,4) = [ diag0m(:) ; zeros(N,1) ];
        A = spdiags( ds , [ N+[-1 0 +1] -N ] , 2*N , 2*N );
    end

end
if strcmp( ModePol , 'TM' )

    if strcmp( Formulation , 'std' ) %"Standard" NON-working approach...
        
        %Equation: n^2*diff( 1/n^2*diff( Hy )  )/dx^2 + ( k0*n(x) )^2*Hy = beta^2*Hy
        % * PMLs: ns^2/s*diff( 1/s/n^2 * diff(Hy) )/dx^2 + ...
        
        %REVISED H-field formulation: No-good results. 
        % I have noted some sensitivity the the [ES=n^2] in numenator. if we use [ES] of node:
        % (a) "0" for both T&R ("correct") => Converges to Unphysical/wrong neff & field
        % (b) "+/-" for both T&R => Converges to correct neff but Field=Ex not Hy(!)->i.e discontinuues
        % (c) "+/-" for T and "0" for R (and vice-versa)=> NO convergence.
        %
        % My implementation of the FDM for the operator matches the (c) case. The FDM params 
        % suggested by Huang&Xu 1993 classic FD-BPM scheme match the (a) case. The only somewhat
        % "functional" case, (b), coincides with the operator used in a BPM for the Ex field of an 
        % infinite slab mode... In all cases the PML effect seems weak, ie it doesnt "hinder".
        FDM_Variation = 2;
        switch FDM_Variation
            case 1, %FDM application: 1/u -> 2/(u_next+u_prev) 
                T_next =     4 .* es ./ sx ./ ( sx_n + sx ) ./ ( es_n + es );  % "Trans" coeff, "next-node"
                T_prev =     4 .* es ./ sx ./ ( sx_p + sx ) ./ ( es_p + es );  % "Trans" coeff, "prev-node"
                R_next = 1 + 4 .* es ./ sx ./ ( sx_n + sx ) ./ ( es_n + es ) ;  % "Refl" coeff, "next-node"
                R_prev = 1 + 4 .* es ./ sx ./ ( sx_p + sx ) ./ ( es_p + es ) ;  % "Relf" coeff, "prev-node"
            case 2, %FDM application: 1/u -> 0.5*(1/u_next+1/u_prev) 
                T_next =     0.5 * es ./ sx .* ( 1./sx_n./es_n + 1./sx./es ); % "Trans" coeff, "next-node"
                T_prev =     0.5 * es ./ sx .* ( 1./sx_p./es_p + 1./sx./es ); % "Trans" coeff, "prev-node"
                R_next = 1 + 0.5 * es ./ sx .* ( 1./sx_n./es_n + 1./sx./es ); % "Refl" coeff, "next-node"
                R_prev = 1 + 0.5 * es ./ sx .* ( 1./sx_p./es_p + 1./sx./es ); % "Refl" coeff, "prev-node"
            case 3, %FDM applied on H-wave-equation suggested by Huang&Xu(1993)->eq(14)
                T_next = 1 - (es_n-es_p)/4./es ;
                T_prev = 1 + (es_n-es_p)/4./es ;
                R_next = 1 + 1 ; 
                R_prev = 1 + 1 ; 
        end
                
        %Sparse diagonals
        diag00 = (ns*k0).^2 + (1-R_next+1-R_prev)/dx^2; %(ns*k0).^2 + (-T_next-T_prev)/dx^2 
        diag1p = T_next/dx^2;
        diag1m = T_prev/dx^2;
        A = spdiags( [ diag1m(:) , diag00(:) , diag1p(:) ] , [-1 0 +1] , N , N );

    elseif strcmp( Formulation , 'hyb' ) %{Hy;Ex} dual-transverse component approach... No results
        %U/R diagonals:
        %diag1p = 0.5 * 1/dx^2 .* ( 1./sx_n./es_n + 1./sx./es  ) /wm/ep0 ; % "Trans" coeff, "next-node"
        %diag1m = 0.5 * 1/dx^2 .* ( 1./sx_p./es_p + 1./sx./es  ) /wm/ep0; % "Trans" coeff, "prev-node"
         diag1p = 4 ./ ( sx_n + sx ) ./ ( es_n + es ) /dx^2/wm/ep0 ; % "Trans" coeff, "next-node"
         diag1m = 4 ./ ( sx_p + sx ) ./ ( es_p + es ) /dx^2/wm/ep0 ; % "Trans" coeff, "prev-node"
         auxi = -diag1p-diag1m ;
        
%          diag1p = ( 1./es_n - (es_n-es_p)/4 ./es_n.^2 )/dx^2/wm/ep0;
%          diag1m = ( 1./es_p + (es_n-es_p)/4 ./es_p.^2 )/dx^2/wm/ep0;
%          auxi = -2./es/dx^2/wm/ep0;
        
        diag00 = sx.*wm*mu0 + auxi;

        %D/L diagonal:
        diag0m = ns.^2./sx*wm*ep0;

        %Prepare diagonals for sparse matrix
        ds(:,1) = [ zeros(N,1) ; diag1m(:) ];
        ds(:,2) = [ zeros(N,1) ; diag00(:) ];
        ds(:,3) = [ zeros(N,1) ; diag1p(:) ];
        ds(:,4) = [ diag0m(:) ; zeros(N,1) ];
        A = spdiags( ds , [ N+[-1 0 +1] -N ] , 2*N , 2*N );
        
    elseif strcmp( Formulation , 'alt' ) %[Ex,Ez'] instead of Hy only...
        
        %Form "semi-diagonals" of A-matrix (length=N)
        dAxx_0 = k0^2*es./sx;
        dAxz_n = +1j./sx/2/dx;
        dAxz_p = -1j./sx/2/dx;
        dAzz_n =  2./(sx_n+sx)/dx^2;
        dAzz_0 = -2./(sx_n+sx)/dx^2-2./(sx_p+sx)/dx^2 + k0^2*es.*sx;
        dAzz_p =  2./(sx_p+sx)/dx^2;
        
        %Insert all diagonals in 2*N length
        ds(:,1) = [ dAxx_0(:) ; dAzz_0(:) ];
        ds(:,2) = [ zeros(N,1) ; dAxz_n(:) ];
        ds(:,3) = [ zeros(N,1) ; dAxz_p(:) ];
        ds(:,4) = [ zeros(N,1) ; dAzz_n(:) ];
        ds(:,5) = [ zeros(N,1) ; dAzz_p(:) ];
        A = spdiags( ds , [ 0 , N+1 , N-1 , +1 , -1 ] , 2*N , 2*N ); clear ds
        
        %Form "semi-diagonals" of N-matrix (length=N)
        dBxx_0 = 1./sx;
        dBzx_n = -1j*1./sx_n/2/dx;
        dBzx_p = +1j*1./sx_p/2/dx;
        
        %Insert all diagonals in 2*N length
        ds(:,1) = [ dBxx_0(:) ; zeros(N,1) ];
        ds(:,2) = [ dBzx_n(:) ; zeros(N,1) ];
        ds(:,3) = [ dBzx_p(:) ; zeros(N,1) ];
        B = spdiags( ds , [ 0 , -N+1 , -N-1 ] , 2*N , 2*N );
        
    end

end

% For 'std' and 'hyb' formulations, B is an "eye" (unity) sparse-array
if ~strcmp( Formulation , 'alt' )
   B = speye(size(A)); 
end

%=========================================================================
% Numerically Solve EigenValue Problem
%=========================================================================

% Upper/lower bounds of neff
ub_neff = max( real(ns) ) - eps;
lb_subs = min( real( ns([1 end]) ) ); % mode should not "oscillate" inside L/R substrates
lb_neff = max( [min(real(ns)) , lb_subs] ) + eps;

% Solve w EIGS -- Inputs:{ How many, Close-to-what-eigval }
if strcmp( NumSolv , 'eigs' )

    if strcmp( Formulation , 'hyb' );
        beta0 = mean( [ub_neff lb_neff] )*k0;
        [V,D] = eigs_mod( A , B , NEV , beta0 );
        neff = diag(D,0)/k0;
    else
        beta0sqr = ( mean( [ub_neff lb_neff] ) *k0 )^2;
        %[V,D] = eigs_mod( A , B , NEV , beta0sqr );
        [V,D] = eigs( A , B , NEV , beta0sqr );
        neff = sqrt( diag(D,0) )/k0;
    end

end

% Solve w SPTARN, Inputs:{ upper- & lower-eigval bound }
if strcmp( NumSolv , 'sptarn' )

    if strcmp( Formulation , 'hyb' );
        lb = lb_neff * k0 ;
        ub = ub_neff * k0 ;
        [V,D] = sptarn_mod( A , B , lb , ub );%, 0 , 1e2*eps , 100 , 20 );
        neff = D / k0 ;
    else
        lb = (k0*lb_neff)^2 ;
        ub = (k0*ub_neff)^2 ;
        [V,D] = sptarn_mod( A , B , lb , ub );%, 0 , 1e2*eps , 100 , 1 );
        neff = sqrt(D)/k0;
    end

end


%-------------------------------------------------------------------------
% Post Process
%-------------------------------------------------------------------------

[dummy,sis] = sort(real(neff),'descend'); % High-neff modes first
sis = sis( dummy >= lb_neff ); % Keep only guided modes (neff >= nsub)
neff = neff(sis); % Re-order EigVals
V = V(:,sis); % Re-order EigVects
clear dummy sis;

if isempty(D), disp(' ## No EigenValues Converged...' ); return; 
elseif isempty(neff), disp(' ## No VALID (Re>1) EigenValues Found...' ); end

if nargin ~= 0
    return
end


%-------------------------------------------------------------------------
% Plots
%-------------------------------------------------------------------------

% Plot index-profile & PML profile
figure('NumberTitle','off','Name','FIDODIES | n(x) profile' );
subplot(2,1,1); plot( x , real(ns) , 'ko-' );
title( 'Refractice-Index Profile: Re\{ n(x) \} ' ); xlabel( 'x (um)' );
subplot(2,1,2); plot( x , -imag(sx) , 'ro-' );
title( 'PML absoprtion profile: -Im\{ sx(x) \}' ); xlabel( 'x (um)' );

% Plot mode (eigenvector) profiles
N_modes = length(neff); % How many plots?
figname = sprintf( 'Mode=%s | Formulation=%s | Solver=%s' , ...
    ModePol, Formulation , NumSolv );
figure('NumberTitle','off','Name',figname );
for k = 1 : N_modes

    % Prepare subplot axes
    nc = ceil( sqrt(N_modes) ); % subplot's number-of-collumns
    nr = ceil( N_modes/nc ); % subplot's number-of-rows
    subplot(nr,nc,k); hold on; box on;

    % Plot PML regions with patches
    patch( x( 1 )+wPML*[0 0 1 1] , [-1 1 1 -1] , 0.8*[1 1 1]); % Left-PML
    patch( x(end)-wPML*[1 1 0 0] , [-1 1 1 -1] , 0.8*[1 1 1]); % Right-PML
    
    % Plot x-Interfaces
    for xIF = xInterfaces
        plot( xIF*[1 1], [-1 1], 'k--' );
    end

    % Plot mode-profile
    VTP = real(V(1:N,k));
    if mean(VTP) < 0, R=-1; else R=1; end % "Reflect" negative-turned EigVector
    plot( x , R*VTP/max(abs(VTP)), 'ro' , 'MarkerSize' , 2 , 'MarkerFaceColor' , 'w' );

    % Cosmetics
    title( sprintf( '{\\bf{[ %02d ]}} n_{eff} = %6.4f' , k , neff(k) ) );
    plot( x , 0*x , 'k:' );
    set(gca,'YLim',1*[-1 1])
    set(gca,'XLim',[min(x) max(x)])
    xlabel( 'x-axis (um)');
    ylabel( 'Fy' );
    
end

% FixMultiFigsPos -- Function to arrange matlab figures (not included)

% % Evalute EigProblem solution, by directly replacing it in the PDE-to-solve
% if strcmp( ModePol , 'TM' )
%     figure
%     H=VTP.';
%     T1 = ( H(1:N-2)-2*H(2:N-1)+H(3:N) )/dx^2;
%     T2 = k0^2*(es(2:N-1)-neff(1)).*H(2:N-1);
%     T3 = -1./es(2:N-1)/4/dx^2.*( es(3:N)-es(1:N-2) ).*( H(3:N)-H(1:N-2) );
%     U  = T1+T2+T3;
%     plot( x(2:N-1) , real(U) , 'bo-' )
%     set(gca,'YLim',[-1 1]); title( 'PDE using EigVal & EigVect (should be all-zero!)' )
% elseif strcmp( ModePol , 'TE' )
%     figure
%     E=VTP.';
%     T1 = ( E(1:N-2)-2*E(2:N-1)+E(3:N) )/dx^2;
%     T2 = k0^2*(es(2:N-1)-neff(1)).*E(2:N-1);
%     U  = T1+T2;
%     plot( x(2:N-1) , real(U) , 'bo-' )
%     set(gca,'YLim',[-1 1]); title( 'PDE using EigVal & EigVect (should be all-zero!)' )
% end


