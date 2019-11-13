function [rLC,n] = LC(wl, na, d, V, N, r_start)

    n = zeros(1,N);
    for i = 1:N
        n(i) = nLC(V,d,i*d/N);
    end
    
    r = zeros(1,N);
    r(1) = (na-n(1))/(na+n(1));
    for i = 2:N
        r(i) = (n(i-1)-n(i))/(n(i-1)+n(i));
    end
    
    rLC = r_start;
    for j = N:-1:1
        ph = (2*pi*d*n(j)/N)./wl;  % phase = ki*di
        rLC = (r(j)+rLC.*exp(-2i*ph))./(1+r(j)*rLC.*exp(-2i*ph));
    end


        
        