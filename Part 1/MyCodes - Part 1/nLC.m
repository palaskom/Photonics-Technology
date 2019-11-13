function neff = nLC(V,d,z)

    no = 1.5;  ne = 1.7; 
    theta_max = 40*V/3 + 10;
    theta = theta_max*sin(pi*z/d);
    
    neff = (no*ne)/sqrt((no*cosd(theta))^2+(ne*sind(theta))^2);

end

