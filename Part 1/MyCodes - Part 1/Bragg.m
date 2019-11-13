function rB = Bragg(wl, wlo, na, nH, nL, M, r_start)

    b = (pi/2)*(wlo./wl); % phase k_i*d_i
    
    % elementary reflection coefficients
    r = zeros(1,M);
    r(1) = (na-nH)/(na+nH);
    r_mid = (nH-nL)/(nH+nL);
    for i = 2:2:length(r)-1 % = M-1: even number
        r(i) = r_mid;
        r(i+1) = -r_mid;
    end
    
    rB = r_start; % == r_(M+1) (vector)
    for i = length(r):-1:1
        rB = (r(i)+rB.*exp(-2i*b))./(1+r(i)*rB.*exp(-2i*b));
    end
    
end

