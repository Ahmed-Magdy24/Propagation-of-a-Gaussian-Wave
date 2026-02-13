%%% Defining Parameters %%%
Wo     = 0.01;                   % Beam Waist at z = 0
lambda = 0.0046;                 % Wavelength
zo     = pi * Wo^2 / lambda;     % Rayleigh distance
k      = 2*pi/lambda;            % Wave number
N      = 512;                    % Number of sampling points
L      = 0.3;                    % Spatial domain size
z      = 5*zo;                      % the distance the beam propagates

%%% 2D Spatial Coordinate System %%%
x = linspace(-L/2, L/2, N);
y = x;
[X, Y] = meshgrid(x, y);
rho_square = X.^2 + Y.^2;

%%% Spatial Frequency Grid for Kx and Ky %%%
dk = 2*pi/L;
kx = (-N/2:N/2-1) * dk;
ky = (-N/2:N/2-1) * dk;
[KX, KY] = meshgrid(kx, ky);

Uo = exp(-rho_square/(Wo^2));  % Complex Amplitude at z = 0

%%% Performing the spatial Fourier transform for Uo %%%
U_kx_ky = fftshift(fft2(ifftshift(Uo)));

kz = k - (KX.^2 + KY.^2)/(2*k);            % kz using paraxial approximation
H = exp(1j * kz * z);                      % The system transfer function
Uz_r = fftshift(ifft2(ifftshift(U_kx_ky .* H)));    % Field distribution at z




%%% Plotting Intensity versus x and spot intensity %%%

fixed_Intensity_x = [-30 30];  % Fixed x-axis for comparing peak intensities
fixed_Intensity_y = [0 1.2];  % Fixed y-axis for comparing peak intensities

figure(1);

%Intensity versus x 
subplot(1,2,1);
plot(x*1000, abs(Uz_r(N/2, :)).^2, 'LineWidth', 1.7);    %Convert to mm for display  , %center row slice at  y= 0: Uz_r(N.2, :)
xlabel('x (mm)'); ylabel('I');
title('Intensity versus x at z');
grid on;
ylim(fixed_Intensity);
xlim(fixed_Intensity_x) ;

%Intensity Spot 
subplot(1,2,2);
imagesc(x*1000, y*1000, abs(Uz_r).^2);
axis image; 
axis xy;
xlabel('x (mm)'); ylabel('y (mm)');
title('Intensity spot at z');
colorbar;






%%% Reflection from Parabolic Mirror %%%
f = -4*zo;                               % Focus of the mirror
M = exp(1j * k * rho_square / (2*f));    % Transfer Function of the mirror
Uout_r = M .* Uz_r;
Uout_r = fftshift(fft2(ifftshift(Uout_r)));  % Output field after reflector

% Propagation after reflection
z1 = zo;
z2 = 4*zo;
z3 = 6*zo;

H1 = exp(1j * kz * z1);
Uz_zo = fftshift(ifft2(ifftshift(Uout_r.*H1)));    % Field at z = zo

H2 = exp(1j * kz * z2);
Uz_4zo = fftshift(ifft2(ifftshift(Uout_r.* H2)));  % Field at z = 4zo

H3 = exp(1j * kz * z3);
Uz_6zo = fftshift(ifft2(ifftshift(Uout_r.*H3)));   % Field at z = 6zo







%%% Plotting beam after reflector at z = zo, 4zo, 6zo %%%
fixed_ylim = [0 0.1];       % Fixed y-axis for comparing peak intensities
fixed_xlim = [-75 75];      % Fixed x-axis for comparing peak intensities

figure(2);

% z = zo
subplot(2,3,1);
plot(x*1000, abs(Uz_zo(N/2,:)).^2,'LineWidth', 1.7);             %Convert to mm for display  , %center row slice at  y= 0: Uz_zo(N.2, :)
xlabel('x (mm)'); ylabel('I (w/m^2)');
title('Intensity at z = zo after reflector');
grid on;
ylim(fixed_ylim);
xlim(fixed_xlim);

% z = 4zo
subplot(2,3,2);
plot(x*1000, abs(Uz_4zo(N/2,:)).^2,'LineWidth', 1.7);            %Convert to mm for display  , %center row slice at  y= 0: Uz_4zo(N.2, :)
xlabel('x (mm)'); ylabel('I (w/m^2)');
title('Intensity at z = 4zo after reflector');
grid on;
ylim(fixed_ylim);
xlim(fixed_xlim);

% z = 6zo
subplot(2,3,3);
plot(x*1000, abs(Uz_6zo(N/2,:)).^2, 'LineWidth', 1.7);           %Convert to mm for display  , %center row slice at  y= 0: Uz_6zo(N.2, :)
xlabel('x (mm)'); ylabel('I (w/m^2)');
title('Intensity at z = 6zo after reflector');
grid on;
ylim(fixed_ylim);
xlim(fixed_xlim);







