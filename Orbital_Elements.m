%% Orbital Elements Calculator; Chandler Hutchens 

format compact ;
format long ;
close all ;
clear all ;
clc ;

%% Intial Position and Velocity Vectors

% Example Data taken from ISS Obit on 2022-03-01T08:53:46.784 When Passing over Ground Station 
r_ijk = [-6245.303275424280 2252.116482750380 1460.824266982890]*1000 ; % m
v_ijk = [-0.44413750167854 -5.00474016070674 5.78160877374074]*1000 ; % m/s

%% Calculations 

% Given fucntion from MATLAB
% Includes numerous Elements but excludes a few, which are calculated below
[a, ecc, incl, RAAN, argp, nu, truelon, arglat, lonper] = ijk2keplerian(r_ijk, v_ijk) ;
fprintf('Semi Major Axis (a) = %.5f kilometers \n', a/1000)
fprintf('Eccentricity (e) = %.5f \n', ecc)
fprintf('Inlcination (i) = %.5f degrees \n', incl)
fprintf('Right Ascenasion of Ascending Node (RAAN) = %.5f degrees \n', RAAN)
fprintf('Angle Between Object Ascending Node and Periapsis (arpg) = %.5f degrees \n', argp)
fprintf('Angle Between Periapsis and Current Position of Object (nu) = %.5f degrees \n', nu)
fprintf('Angle Between x-axis and Object Position Vector (truelon) = %.5f degrees \n', truelon)
fprintf('Angle Between Ascending Node and Object Position Vector (arglat) = %.5f degrees \n', arglat)
fprintf('Angle Between x-axis and Eccentricity Vector (lonper) = %.5f degrees \n', lonper)

r = (r_ijk/1000) ; % km
v = (v_ijk/1000) ; % km
mu = 3.986e5 ; % Mu_Earth in km^3/s^2

magv = norm(v) ; % Magnitude of velocity
magr = norm(r) ; % Magnitude of r

h = cross(r,v) ; % Anglular Momentum 
magh = norm(h) ; % Magnitude of h 
nh = h/magh ; % Normalize the angular momentum

ni = [0 0 1] ;
n = cross(ni,h); 
magn = norm(n) ; % Magnitude of n

e = (1/mu)*cross(v,h) - (r/magr) ; % Eccentricity 
mage = norm(e) ; % Magnitude of Eccentricty
fprintf('Confirmation: Eccentricity (e) = %.5f degrees \n', mage)

E = magv^2/2 - (mu/magr) ; % Total Energy

a2 = -1/(2*E) ; % Semi Major Axis
fprintf('Confirmation: Semi Major Axis (a) = %.5f kilometers \n', a/1000)

i = acosd(dot(nh,ni)) ; % Inclination
fprintf('Confirmation: Inlcination (i) = %.5f degrees \n', i)

omega = acosd(dot([1 0 0], n/magn)) ; % RAAN
check1 = dot([0 1 0],n) ; % Ambiugity check

    if check1 < 0
        omega = 360 - omega ;
    else omega = omega ;
    end 

fprintf('Confirmation: Right Ascenasion of Ascending Node (RAAN) = %.5f degrees \n', omega)

w = acosd((dot(n,e))/(magn*norm(e))) ; % Argument of Perigee
check2 = dot(r,v) ; % Ambiugity check

 if check2 < 0
        w = 360 - w ;
 else w = w ;
 end 

fprintf('Argument of Perigee (w) = %.5f degrees \n', w)

fo1 = acosd((dot(r,e)/(magr*mage))) ; % True Anamomly 

 if fo1 < 0 
        fo1 = 2*pi - fo1 ;
 else fo1 = fo1 ;
 end 

fprintf('True Anamomly (f) = %.5f degrees \n', fo1)

fo = fo1*(180/2*pi) ; % Convert to radians for E

Eo1 = 2*atan((sqrt((1-mage)/(1+mage)))*tan(fo/2)) ; % Eccentric anomaly 
Eo = Eo1*(2*pi/180) ;
fprintf('Eccentric Anamomly (E) = %.5f degrees \n', Eo)

Mo1 = Eo - mage*sin(Eo) ; % Mean anomaly 
Mo = Mo1*(2*pi/180) ;
fprintf('Mean Anamomly (M) = %.5f degrees \n', Mo)
