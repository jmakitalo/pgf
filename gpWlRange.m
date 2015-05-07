function gpWlRange(filename, dx, dy, dz, phi, pwtheta, pwphi, wlst, ri_name, ri0)
%GPWLRANGE   Compute and save smoothed periodic Green's function terms.
%   filename: name of output data file.
%   dx,dy: periods of the array along two array axes.
%   dz: maximum distance in z-direction between source and observation
%       points.
%   phi: angle (radians) between the first array axis and the x-axis.
%   pwtheta: angle (radians) between incident wave vector and z-axis.
%   pwphi: angle (radians) between incident wave vector projected to 
%          xy-plane and the x-axis.
%   wlst: list of wavelengths.
%   ri_name: name of file for medium refractive index data. Can also be
%            directly refractive index value if not wavelength dependent.
%   ri0: refractive index of the incident wave medium.
%
%   This function is used to precompute the periodic Green's function (GF)
%   of the Helmholtz equation (2D periodicity in 3D space).
%   The summants in the Ewald series are precomputed for linear
%   interpolation.
%   The GF is the smoothed by subtraction of first two odd Taylor series 
%   terms.
%
%   The output data is used by the project for90-mom2 hosted at
%   https://github.com/jmakitalo/for90-mom2
%   See the module greenprd.f90.
%
%   Author: Jouni Mäkitalo (Tampere University of Technology).

% Number of wavelengths to consider in the precomputation.
nwl = length(wlst);

sp = sin(phi);
cp = cos(phi);

% Refractive index data.
if(~isnumeric(ri_name))
  ridata = load(ri_name);
  for n=1:nwl
    ri_interp(n) = interp1(ridata(:,1), ridata(:,2), wlst(n)) +...
	1i*interp1(ridata(:,1), ridata(:,3), wlst(n));
  end
end

% Open file for writing.
fid = fopen(filename, 'w');

% Write header.
fprintf(fid, 'type 2dinterp1d\n');
fprintf(fid, 'periods %.15g %.15g %.15g %.15g\n', dx, dy, dz, phi);
fprintf(fid, 'planewave %g %g\n', pwtheta, pwphi);
fprintf(fid, 'nwl %d\n', nwl);

% Setup parameters for each wavelength.
for nw=1:nwl
    disp(['Wavelength ' num2str(nw) ' of ' num2str(nwl)]);
    
    % Splitting parameter of the Ewald series.
    E = sqrt(pi/(dx*dy*cp));
    
    % Wavelength (vacuum).
    wl = wlst(nw);

    % Set refractive index.
    if(isnumeric(ri_name))
      ri = ri_name;
    else
      ri = ri_interp(nw);
    end
    
    % Wavenumber in the medium.
    k = ri*2*pi/wl;
    
    % Wavevector (tangential) of the incident plane-wave.
    % This is the same in all periodic media and is determined
    % by the incident waves.
    k0 = [sin(pwtheta)*cos(pwphi),...
        sin(pwtheta)*sin(pwphi), cos(pwtheta)]*ri0*2*pi/wl;
    
    % Maximum number of terms in the Ewald series.
    newald = 2;
    
    %if(wlst(nw)<ri*min([dx,dy]))
    %    newald = 3;
    %end
    
    N = newald;
    M = newald;
    
    % Number of evaluation points.
    t = (wlst(nw)-4e-7)/12e-7;
    npoints = round(t*1000 + (1-t)*6000);
    np = (max([N M]) + 1)*npoints;
    npz = npoints;
    
    % Evaluation range.
    range = sqrt(((N+1)*dx)^2 + ((M+1)*dy)^2 + dz^2) + 10e-9;
    
    Rnm = linspace(0, range, np);
    z = linspace(-dz, dz, npz);
    
    % Write header for this wavelength.
    fprintf(fid, 'wl %d\n', nw);
    fprintf(fid, 'wavelength %.15g\n', wl);
    fprintf(fid, 'ri (%.15g, %.15g)\n', real(ri), imag(ri));
    fprintf(fid, 'ri0 (%.15g, %.15g)\n', real(ri0), imag(ri0));
    fprintf(fid, 'E %.15g\n', E);
    fprintf(fid, 'npoints %d\n', np);
    fprintf(fid, 'npointsz %d\n', npz);
    fprintf(fid, 'nterms %d %d\n', N, M);
    fprintf(fid, 'range %.15g\n', range);
    
    % Reserve data arrays.
    dspat = zeros(2,np,2*M+1,2*N+1);
    dspec = zeros(3,npz,2*M+1,2*N+1);
    
    % Compute data.
    for n=-N:N
        for m=-M:M
            ni = N+n+1;
            mi = M+m+1;
            
            dspat(1,:,mi,ni) = gpsSpat(Rnm, n, m, k, k0, E, dx, dy, sp, cp);
            dspat(2,:,mi,ni) = gradGpsSpat(Rnm, n, m, k, k0, E, dx, dy, sp, cp);
            
            dspec(1,:,mi,ni) = gpsSpec(z, n, m, k, k0, E, dx, dy, sp, cp);
            [dspec(2,:,mi,ni) dspec(3,:,mi,ni)] = gradGpsSpec(z, n, m, k, k0, E, dx, dy, sp, cp);
        end
    end
    
    writeComplexArray(fid, dspat);
    writeComplexArray(fid, dspec);
end

fclose(fid);

end
