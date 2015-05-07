function gradGspat = gradGpsSpat(Rnm, n, m, k, k0, E, dx, dy, sp, cp)
%GRADGPSSPAT   Gradient of periodic Green's function for Helmholtz operator
%              with 2D periodicity in 3D space. Term (n,m) in spatial Ewald
%              series. Singularities are subtracted.
%              Assumes time-dependency exp(-iwt).
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
    a = a + exp(l*1i*k*Rnm).*((1 - l*1i*k*Rnm)./(Rnm.^3).*erfcz(Rnm*E + l*1i*k/(2*E))...
        + 2*E./(sqrt(pi)*Rnm.^2).*exp(-(Rnm*E + l*1i*k/(2*E)).^2));
end

gradGspat = a.*exp(1i*dot(k0,rho))/(8*pi);

% Remove singularities.
if(abs(n)<2 && abs(m)<2)
    gradGspat = gradGspat - (1./(4*pi*Rnm.^3) + k^2./(8*pi*Rnm))*exp(1i*dot(k0,rho));
end

% Impose evaluation of (r-r') as ((r-r')/R)*R.
gradGspat = gradGspat.*Rnm;

% The first entry in Rnm may be zero, in which case a limiting value
% is used.
if(Rnm(1)==0)
    gradGspat(1) = 0;
end

end