function Gspat = gpsSpat(Rnm, n, m, k, k0, E, dx, dy, sp, cp)
%GPSSPAT   Periodic Green's function for Helmholtz operator with 2D 
%          periodicity in 3D space. Term (n,m) in spatial Ewald series.
%          Singularities are subtracted.
%          Assumes time-dependency exp(-iwt).
%   Rnm: distance between observation and source points in array.
%   n,m: array source indices.
%   k: wavenumber.
%   k0: incident wave vector.
%   E: Ewald splitting parameter.
%   dx,dy: array periods.
%   sp,cp: sine and cosine of array axis angle.
%
%   Author: Jouni Mäkitalo (Tampere University of Technology).

rho = [n*dx*cp m*dy+n*dx*sp 0];

a = 0;
for l=[-1,1]
  a = a + exp(l*1i*k*Rnm).*erfcz(Rnm*E + l*1i*k/(2*E));
end

Gspat = a*exp(1i*dot(k0,rho))./(8*pi*Rnm);

% Remove singularities.
if(abs(n)<2 && abs(m)<2)
    Gspat = Gspat + (-1./(4*pi*Rnm) + k^2*Rnm/(8*pi))*exp(1i*dot(k0,rho));
end

% The first entry in Rnm may be zero, in which case a limiting value
% is used.
if(Rnm(1)==0)
    Gspat(1) = (2*1i*k*(erfcz(1i*k/(2*E)) - 1) ...
        -4*E/sqrt(pi)*exp(k^2/(4*E^2)))*exp(1i*dot(k0,rho))/(8*pi);
end

end