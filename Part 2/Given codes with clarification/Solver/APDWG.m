function [ neff_TE , neff_TM , Ey_TE , Ex_TM ] = APDWG( wl , h , n1 , n2 , n3 , x )

% Alexandros Pitilakis
% Thessaloniki, September/October 2010

%Test Input Vars:
if nargin == 0
        
    close all; clc;
    
%     n1 = 3.60; % Index : Guiding Layer (Silicon)
%     n2 = 3.24; % Index : Substratre (Oxide)
%     n3 = 3.24; % Index : Cladding (Air)
%     wl = 1.55; %[um] Wavelength
%     h  = 0.40; %[um] Slab Thickness
%     x  = linspace( -4*h , 5*h , 500 ); %[um] transverse-coordinate vector

    % Low index contrast
    n1 = 1.456;
    n2 = 1.445;
    n3 = n2;
    wl = 1.55;
    h = 6.0;
    x = linspace( -4*h , 5*h , 500 ); 
     
%     %DLSPP
%     n1 = 1.53 - 0*3e-4; n2 = 0.55-100j; n3 = 1.00; %
%     wl = 1.566; h  = 0.54; %[um] Slab Thickness
%     x  = linspace( -4*h , 5*h , 5000 ); %[um] transverse-coordinate vector
 
%     %Plamonics
%     nm = j*11.5; 
%     n2 = 1.00 ; n1 = nm ; n3 = 1.00; 
%     wl = 1.555; h  = 0.01; %[um] Slab Thickness
%     x  = linspace( -4*h , 5*h , 5000 ); %[um] transverse-coordinate vector
%   
%     %SOI-Wire
%     n1 = 3.2; n2 = 1.45; n3 = 1.00 ;
%     wl = 1.55; h  = 0.40; %[um]
%     x = linspace( -1 , 1 , 1024 );

%     %MDLMs
%     n1 = 1.5;
%     n2 = 0.55-j*11.5; %-sqrt( -17 );
%     n3 = n2;
%     wl = 0.6328;
%     h = 0.026;
%     x = h*linspace( -3 , 4 , 1024 );
    

%     %DLSPP w/ EIM
%     n1 = 1.4497 - 0.0028*j; %neff for TM-Slab:(air/PMMA-600nm/Palik-Gold) @ 1.55um
%     n2 = 1.0038 - 0.0004*j; %nspp for Palik-Gold/Air Interface
%     n3 = n2;    wl = 1.55;    h = 0.8; x  = linspace( -4*h , 4*h , 500 );
%     clc;

end

% ========================================================================
% Solve for n_effective
% ========================================================================

%Some Aux-Variables
k0 = 2*pi / wl; % [um^-1] Wavenumber, free-space
ne = linspace( max(real(n2),1)+0.01 , n1-0.01 , 10000 ); % A discrete n_eff range

%ne = linspace( 1 , 1.5 , 1e5 );

%Characteristic Equation Parameters
k  = k0 * sqrt( n1^2 - ne.^2 ); %common for TE/TM
g  = k0 * sqrt( ne.^2 - n2^2 ); %TE modes
gt = g * (n1/n2)^2;             %TM modes
d  = k0 * sqrt( ne.^2 - n3^2 ); %TE modes
dt = d * (n1/n3)^2;             %TM mode

%Form the Chacteristic Eqs
XE  = tan( k * h ) - k.*(g +d )./(k.^2-g.*d  ); %TE
XEt = tan( k * h ) - k.*(gt+dt)./(k.^2-gt.*dt); %TM

isreal(XE)
isreal(XEt)

%Find the Roots
if all( isreal( XE ) ) && all( isreal( XEt) ) % Classic "No-Complex" Solving

    neff_TE = interpinv( ne , XE  , 0 , 10);
    neff_TM = interpinv( ne , XEt , 0 , 10);
    
else % Complicated (SPP, Losses etc)
   
    %opts = optimset( optimset('fsolve') , 'MaxFunEvals' , 100 );
    
    %TE-Mode..............................................................
    fprintf( ' ** FSOLVE : Complex EIM , TE-Modes...  ' );
    candidate_x0 = [ mean(real([n1 n2])) , ...
        linspace( max([1 ,min(real([n1 n2 n3]))]) , max(real([n1 n2 n3])) , 100 ); ];
    % candidate_x0 = sqrt( linspace( min(abs(real([n1 n2 n3].^2))) , max(abs(real([n1 n2 n3].^2))) , 100 ) )
    %candidate_x0 = abs(candidate_x0);
    for ksg = 1 : length(candidate_x0)

%          [xroot1,fval1,flag1] = fsolve( @(x) APDWG_CharEq( x , n1 , n2 , n3 , wl , h , 1 ) , ...
%              candidate_x0(ksg) );

        %Newton-Raphson
        x0 = candidate_x0(ksg); Tol_neff = 1e-8; Tol_XE = 1e-8; MaxIter = 100 ;
        for knr = 1 : MaxIter
            [XEa dXEa] = APDWG_CharEq( x0 , n1 , n2 , n3 , wl , h , 1 );
            x0 = x0 - XEa./dXEa;
            if abs(XEa./dXEa) < Tol_neff && abs(XEa)< Tol_XE && ...
                abs(abs(x0)-abs(n2)) > 1e-5 && abs(abs(x0)-abs(n1)) > 1e-5 
                xroot1 = x0; flag1 = 1;
            elseif knr == MaxIter, xroot1 =[]; flag1 = 0;
            end
        end
        %plot( log10(yu) ); error
        %disp( [ abs(XEa./dXEa) abs(XEa)] )
        
        if flag1 == 1, disp(' --> Converged OK'); break; end
    end
    if flag1 == 0, disp( ' ### Did NOT converge' ); end
    
    %TM-Mode..............................................................
    fprintf( ' ** FSOLVE : Complex EIM , TM-Modes...  ' );
    candidate_x0 = [ mean(real([n1 n2])) , ...
        linspace( max([1 ,min(real([n1 n2 n3]))]) , max(real([n1 n2 n3])) , 100 ); ];
    %candidate_x0 = sqrt( linspace( min(abs(real([n1 n2 n3].^2))) , max(abs(real([n1 n2 n3].^2))) , 100 ) );
    %candidate_x0 = abs(candidate_x0);
    for ksg = 1 : length(candidate_x0)
        
        % [xroot2,fval2,flag2] = fsolve( @(x) APDWG_CharEq( x , n1 , n2 , n3 , wl , h , 2 ) , ...
        %     candidate_x0(ksg) );
         
        %Newton-Raphson
        x0 = candidate_x0(ksg); Tol_neff = 1e-2; Tol_XE = 1e-2; MaxIter = 100 ;
        for knr = 1 : MaxIter
            [XEa dXEa] = APDWG_CharEq( x0 , n1 , n2 , n3 , wl , h , 2 );
            x0 = x0 - XEa./dXEa;
            if abs(XEa./dXEa) < Tol_neff && abs(XEa)< Tol_XE && abs(abs(x0)-abs(n2)) > 1e-5 && abs(abs(x0)-abs(n1)) > 1e-5 
                xroot2 = x0; flag2 = 1;
            elseif knr == MaxIter, xroot2 =[]; flag2 = 0;
            end
        end
        
        if flag2 == 1, disp(' --> Converged OK'); break; end
    end
    if flag2 == 0, disp( ' ### Did NOT converge' ); end
        
    %Store n_effective
    if flag1 == 1, neff_TE = xroot1; else neff_TE = []; end
    if flag2 == 1, neff_TM = xroot2; else neff_TM = []; end
    
end

% if ~isempty(neff_TE) && ~isempty(neff_TM) && any( abs(neff_TE - n2) < 1e-5 || abs(neff_TM - n2) < 1e-5 )
% %     neff_TE
% %     neff_TM
% %     n2
%     error( 'Wrong Convergence: neff ~ n_substrate' )
% end


%some Error-Checking
Nmodes_TE = length( neff_TE );
Nmodes_TM = length( neff_TM );
if Nmodes_TE == 0 && Nmodes_TM == 0
    disp( ' ## No guided modes found/supported! Check n_eff discetization....' );
    return
end

% ========================================================================
% Calc Field Distributions
% ========================================================================

% Translate x-vector, to have its center at the middle of the guide
x = x + h/2;

% Region-Separation boolean x-vectors for Ey(x) and Hy(x).
is1 = x>=0 & x<=h ; %Guiding Layer
is2 = x<0         ; %Substrate
is3 = x>h         ; %Cladding

% TE-Modes
Ey_TE = NaN * zeros( Nmodes_TE , length(x) );
for kkm = 1 : Nmodes_TE
    neff = neff_TE(kkm);
    k  = k0 * sqrt( n1^2 - neff^2 );
    g  = k0 * sqrt( neff^2 - n2^2 );
    d  = k0 * sqrt( neff^2 - n3^2 );
    phi = atan( g / k );
    A1 = 1 ;
    A2 = A1 * cos( phi );
    A3 = A1 * cos( k*h -phi );
    Ey_TE( kkm , is1 ) = A1 * cos( k * x(is1) - phi );
    Ey_TE( kkm , is2 ) = A2 * exp( g * x(is2) );
    Ey_TE( kkm , is3 ) = A3 * exp( -d *(x(is3)-h) );
end

%The Refractive Indices across the x-vector cross-section
ns = NaN*zeros(1,length(x));
ns( is1 ) = n1;
ns( is2 ) = n2;
ns( is3 ) = n3;

% TM-Modes
Ex_TM = NaN * zeros( Nmodes_TM , length(x) );
Hy = NaN * zeros( 1 , length(x) );
for kkm = 1 : Nmodes_TM
    neff = neff_TM(kkm);
    k  = k0 * sqrt( n1^2 - neff^2 );
    g  = k0 * sqrt( neff^2 - n2^2 );
    d  = k0 * sqrt( neff^2 - n3^2 );
    phi = atan( g / k  * (n1/n2)^2 );
    A1 = 1 ;
    A2 = A1 * cos( phi );
    A3 = A1 * cos( k*h -phi );
    Hy( is1 ) = A1 * cos( k * x(is1) - phi );
    Hy( is2 ) = A2 * exp( g * x(is2) );
    Hy( is3 ) = A3 * exp( -d *(x(is3)-h) );
    Ex = Hy ./ ns.^2 ; %actually, its Proportional to, as derived from TM-Mode Electromagnetic theory
    Ex_TM( kkm , : ) = Ex ./ max(abs(Ex)) ;    
end

%Plot Fields:
if nargin == 0
    %close all; clc;
       
    neff_TE
    neff_TM

    subplot(2,1,1); hold on;
    title( sprintf( 'TE-mode : n_{eff} = %1.6f' , neff_TE ) ); ylabel( 'Re( Ey )' ); xlabel( 'x(um)' )
    patch( h/2*[-1 -1 1 1] , [-1 1 1 -1] , 0.7*[1 1 1] )
    if ~isempty(neff_TE)
        plot( x-h/2 , real(Ey_TE) , 'r' , 'LineWidth' , 2 );
    end
    box on;
    set( gca , 'Layer' , 'Top' , 'XLim' , -h/2+[x(1) x(end)] ) %  , 'YLim' , [-1 1]
    
    as = sprintf( 'wl=%1.2f | h=%1.3f | ns=[%1.4f %1.4f %1.4f]' , wl , h , [ n2 n1 n3 ] ); 
    text( x(1) , -0.5 , 1 ,  as , 'BackGround' , 'w' , 'EdgeColor' , 'k' );
    
    subplot(2,1,2); hold on;
    title( sprintf( 'TM-mode : n_{eff} = %1.6f' , neff_TM ) ); ylabel( 'Re( Hy )' ); xlabel( 'x(um)' )
    patch( h/2*[-1 -1 1 1] , [-1 1 1 -1] , 0.7*[1 1 1] )
    if ~isempty(neff_TM)
     plot( x-h/2 , real(Hy) , 'b' , 'LineWidth' , 2 );
    end
    box on;
    set( gca , 'Layer' , 'Top', 'XLim' , -h/2+[x(1) x(end)] ) %  , 'YLim' , [-1 1] 

    return
end

% clc
% close all
% ns = NaN*x;
% ns(is2) = n2;
% ns(is1) = n1;
% ns(is3) = n3;
% dx=mean(diff(x));
% AAA = diff(diff(Ey_TE))/dx^2 + (2*pi/wl)^2*( ns(2:end-1).^2 - neff_TE^2).*Ey_TE(2:end-1);
% 
% BBB = ns(2:end-1).^2.*diff(diff(1./ns.^2.*Hy))/dx^2 + ...
%     (2*pi/wl)^2*( ns(2:end-1).^2 - neff_TM^2).*Hy(2:end-1);
% 
% plot( x(2:end-1) , abs(AAA) , 'bo-' ); hold on;
% plot( x(2:end-1) , abs(BBB) , 'rs-' ); hold on;
% 
% set(gca,'YLim',[0 1e0])
% error


% ========================================================================
% Flip Result Order!
% ========================================================================

%interpinv returns the roots (neffs of supported modes) sort in ascending
%order. We, usually, prefer to call #1 the fundamental mode, who in the
%above ordering comes last, hence the flipping/reversion
neff_TE = fliplr( neff_TE );
neff_TM = fliplr( neff_TM );
Ey_TE = flipud( Ey_TE );
Ex_TM = flipud( Ex_TM );


%-------------------------------------------------------------------------
%Test-Plots
%-------------------------------------------------------------------------
if nargin == 0

    close all; figure('Name','APDWG','NumberTitle','off')
    
    % Characteristic Equation Solving -> n_eff
    subplot(2,1,1); hold on
    plot( real(ne) , real(XE) , 'bs' ); hold on
    plot( real(ne) , real(XEt) , 'ro' );
    plot( real(ne) , zeros(size(ne)) , 'k' , 'LineWidth' , 2 );
    axis( [min(ne) max(ne) 10*[-1 1]] )
    axis on; grid on;
    
    %Annotations
    str1 = 'function( n_{eff} ) = tan(\kappah) - \kappa(\gamma+\delta)/( \kappa^2 - \gamma\delta)';
    ylabel( str1 ); xlabel( 'n_{effective} = \beta / k_0' );

    title( sprintf( 'Asymetric Planar Dielectric Waveguide\nCharacteristic Equation - Numerical Solution' ) , ...
        'FontWeight' , 'Bold' , 'FontSize' , 14 )
    legend( 'TE Modes' , 'TM Modes' , 'Location' , 'NorthWest' );

    text( neff_TE(1) , -0.25 , sprintf( '\\uparrow neff(TE) ~ %6.4f' , neff_TE(1) ) , ...
        'FontWeight' , 'Bold' , 'FontSize' , 14 , 'VerticalAlignment' , 'Top' , ...
        'Color' , 'b' , 'EdgeColor' , 'b' , 'BackGroundColor' , 'w' )
    text( neff_TM(1) , +0.25 , sprintf( '\\downarrow neff(TM) ~ %6.4f' , neff_TM(1) ) , ...
        'FontWeight' , 'Bold' , 'FontSize' , 14 , 'VerticalAlignment' , 'Bottom' , ...
        'Color' , 'r' , 'EdgeColor' , 'r' , 'BackGroundColor' , 'w' )

    % Field-Profiles
    subplot(2,1,2); hold on; grid on;
    plot( x , abs(Ey_TE(1,:)) , 'LineWidth' , 2 , 'Color' , 'b' , 'LineStyle' , '-' );
    %plot( x , Ex_TM(1,:) , 'LineWidth' , 2 , 'Color' , 'r' , 'LineStyle' , '-' );
    %--Plots Ex_TM piecewise to display the discontinuity
    ib = find( abs(diff(Ex_TM(1,:))) > 0.05 );
    plot( x((      1):ib(1)) , abs(Ex_TM(1,(      1):ib(1))) , 'LineWidth' , 2 , 'Color' , 'r' , 'LineStyle' , '-' );
    plot( x((ib(1)+1):ib(2)) , abs(Ex_TM(1,(ib(1)+1):ib(2))) , 'LineWidth' , 2 , 'Color' , 'r' , 'LineStyle' , '-' );
    plot( x((ib(2)+1):end  ) , abs(Ex_TM(1,(ib(2)+1):end  )) , 'LineWidth' , 2 , 'Color' , 'r' , 'LineStyle' , '-' );
    

    legend( 'TE-0 Mode : Ey(x)' , 'TM-0 Mode : Ex(x)' )
    plot( x , 0*x , 'k' );
    if n3 == 1, cCla = 1; else cCla = 0.80; end
    if n2 ~= n3, cSub = 0.75; else cSub = cCla; end
    fill( [min(x) min(x) 0 0] , 1*[-1 1 1 -1] , cSub*[1 1 1] , 'FaceAlpha' , 0.5 ); %Substrate
    fill( [0 0 h h]           , 1*[-1 1 1 -1] , 0.50*[1 1 1] , 'FaceAlpha' , 0.5 ); %Guiding
    fill( [h h max(x) max(x)] , 1*[-1 1 1 -1] , cCla*[1 1 1] , 'FaceAlpha' , 0.5 ); %Cladding
    set( gca,'XLim',[min(x) max(x)], 'YLim' , [0 A1] );
    xlabel( 'x-axis [\mum]' ); ylabel( 'Amplitude E_y(x) and E_x(x)' )

    
%     % ========================================================================
%     % Power Calculations
%     % ========================================================================
% 
%     %It can be proven that Pz(x)=(n_eff/c0/mu0) * |Ey(x)|^2. Using this, is is
%     %easy to analytically integrate the Ey(x) form to get the guided power in
%     %each of the three layers.
%     GuiPow_Anal = NaN*ones(3,2);

%     %Guided-Power Analytical
%     GuiPow_Anal(1,kkk) = (2*h*k + sin(2*h*k-2*phi) + sin(2*phi)) / (4*k);
%     GuiPow_Anal(2,kkk) = 1/2/g * abs(A2/A1)^2;
%     GuiPow_Anal(3,kkk) = 1/2/d * abs(A3/A1)^2;
%     if kkk == 2 % TM
%         GuiPow_Anal(:,kkk) = GuiPow_Anal(:,kkk) ./ [ n1;n2;n3 ].^2;
%     end
%                 
%     %Now, normalize Powers in [%]
%     GuiPow_Anal(:,kkk) = GuiPow_Anal(:,kkk) / sum( GuiPow_Anal(:,kkk) ) * 100;

    
%     %Power-Guied Percentage in each layer
%     str = sprintf( 'Power in Layers [%%]\n%s\n       TE   |   TM\n%s\nP1 =  %4.1f  |  %4.1f\nP2 =  %4.1f  |  %4.1f\nP3 =  %4.1f  |  %4.1f' ,...
%         flwcs('=',20) , flwcs('=',20) , GuiPow_Anal' );
%     text( min(get(gca,'XLim'))+0.1 , max(get(gca,'YLim'))-0.1 , str , 'VerticalAlignment','Top',...
%         'HorizontalAlignment' , 'Left' , 'EdgeColor' , 'k' , 'BackgroundColor' , 'w' , ...
%         'FontSize' , 12 , 'FontName' , 'FixedWidth' , 'Margin' , 5 , 'LineWidth' , 2

end








