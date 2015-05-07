function [gradGspect gradGspecz] = gradGpsSpec(z, n, m, k, k0, E, dx, dy, sp, cp)
%GRADGPSSPEC   Gradient of periodic Green's function for Helmholtz operator
%              with 2D periodicity in 3D space. Term (n,m) in spectral
%              Ewald series.
%              Assumes time-dependency exp(-iwt).
%   z: distance between observation and source point z-coordinates.
%   n,m: array source indices.
%   k: wavenumber.
%   k0: incident wave vector.
%   E: Ewald splitting parameter.
%   dx,dy: array periods.
%   sp,cp: sine and cosine of array axis angle.

kt = [(k0(1) + 2*pi*(n/(dx*cp) - m*sp/(dy*cp))),...
    (k0(2) + 2*pi*m/dy), 0];

% Choose the proper branch of the square root.
kz = sqrt(k^2 - dot(kt,kt));

at = 0;
az = 0;
for l=[-1,1]
    at = at + exp(-l*1i*kz*z).*erfcz(-1i*kz/(2*E) + l*z*E);

    az = az + l*exp(-l*1i*kz*z).*(1i*kz*erfcz(-1i*kz/(2*E) + l*z*E)...
        +2*E/sqrt(pi)*exp(-(-1i*kz/(2*E) + l*z*E).^2));
end

A = dx*dy*cp;

gradGspect = at/(4*A*kz);
gradGspecz = az/(-4*A*1i*kz);

end