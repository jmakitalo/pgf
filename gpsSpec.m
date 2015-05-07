function Gspec = gpsSpec(z, n, m, k, k0, E, dx, dy, sp, cp)
%GPSSPEC   Periodic Green's function for Helmholtz operator with 2D
%          periodicity in 3D space. Term (n,m) in spectral Ewald series.
%          Assumes time-dependency exp(-iwt).
%   z: distance between observation and source point z-coordinates.
%   n,m: array source indices.
%   k: wavenumber.
%   k0: incident wave vector.
%   E: Ewald splitting parameter.
%   dx,dy: array periods.
%   sp,cp: sine and cosine of array axis angle.
%
%   Author: Jouni Mäkitalo (Tampere University of Technology).

kt = [k0(1) + 2*pi*(n/(dx*cp) - m*sp/(dy*cp)),...
    k0(2) + 2*pi*m/dy, 0];

% Choose the proper branch of the square root.
kz = sqrt(k^2 - dot(kt,kt));

a = 0;
for l=[-1,1]
    a = a + exp(-l*1i*kz*z).*erfcz(-1i*kz/(2*E) + l*z*E);
end

A = dx*dy*cp;
Gspec = a/(-4*A*1i*kz);

end