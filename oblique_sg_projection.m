% Codigo para reproyectar lon,lat a coord oblique stereographic (Reerink et
% al. 2010)
function [x,y] = oblique_sg_projection(lambda,phi)
% coord lambda (lon), phi (lat) in degrees check if is 0-360

% there is the possibility to change beta, lamda_M, etc to input and this code
% work for everything..

% Needed variables for ANT domain

beta_deg = 71.0;
lambda_M_deg = 0.0;
phi_M_deg = -90;
earth_radius = 6.371221e6;

% convert longitude from -180 - 180 to 0 - 360 if needed
if lambda < 0
    lambda = 360 + lambda;
end
% Convert beta to alpha
    alpha_deg = 90 - beta_deg;

%  Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180) * phi ;
    lambda_P = (pi / 180) * lambda ;

%     Convert projection parameters to radians:
    lambda_M = (pi / 180) * lambda_M_deg;
    phi_M    = (pi / 180) * phi_M_deg;
    alpha    = (pi / 180) * alpha_deg ;

%    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1 + cos(alpha)) / (1 + cos(phi_P) * cos(phi_M) * cos(lambda_P - lambda_M) + sin(phi_P) * sin(phi_M));

%    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x =  earth_radius * (cos(phi_P) * sin(lambda_P - lambda_M)) * t_P_prime;
    y =  earth_radius * (sin(phi_P) * cos(phi_M) - (cos(phi_P) * sin(phi_M)) * cos(lambda_P - lambda_M)) * t_P_prime;
