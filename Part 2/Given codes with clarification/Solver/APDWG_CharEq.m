function [ XE , dXE ] = APDWG_CharEq( ne , n1 , n2 , n3 , wl , h , Mode )

k0 = 2*pi / wl; % [um^-1] Wavenumber, free-space

%Characteristic Equation Parameters
k  = k0 * sqrt( n1^2 - ne.^2 ); %common for TE/TM
g  = k0 * sqrt( ne.^2 - n2^2 ); %TE modes
gt = g * (n1/n2)^2;             %TM modes
d  = k0 * sqrt( ne.^2 - n3^2 ); %TE modes
dt = d * (n1/n3)^2;             %TM mode

%Form the Chacteristic Eqs
switch Mode
    case 1
        XE  = tan( k * h ) - k.*(g +d )./(k.^2-g.*d  ); %TE
    case 2
        XE = tan( k * h ) - k.*(gt+dt)./(k.^2-gt.*dt); %TM
end

% ------------------------------------------------------------------------
% Used for a Newton-Raphson Method
% ------------------------------------------------------------------------

%Form the Derivatives of k,g,d with respect to neff (ne)
dk = -k0^2*ne./k;
dg = +k0^2*ne./g; %TE
dd = +k0^2*ne./d; %TE
dgt = (n1/n2)^2 * dg; %TM
ddt = (n1/n3)^2 * dd; %TM

if Mode == 2 %TM
    dg = dgt;
    dd = ddt;
end

%Form the parts of the derivative of XE with respect to neff (ne)
Part1 = 1./cos(k*h).^2 * h .* dk;
Part2 = +( dk.*(g+d) + k.*(dg+dd) ) .* ( k.^2 - g.*d );
Part3 = -( k.*(g +d ).*( 2*k.*dk - dg.*d - g.*dd ) ) ;
DN23  = ( k.^2 - g.*d ).^2;

%This is the total derivative
dXE = Part1 - (Part2+Part3)./DN23;





