close all; clc; clearvars;

% Constants
c = 2.998e8; % meter/second
h = 6.626e-34; % meter^2 * kg / s
k_b = 1.38e-23; % J*K^-1
nm = 1e-9;

% Input reading
prompt = {'Enter the temperature in kelvin:','Enter the number of points to simulate:',['' ...
          'Enter noise the noise level in percent:'],['Enter the amount of line broadening (nm):' ...
          ''],'Enter the range of wavelengths (e.g.[x,y] from x to y):', 'Enter the error in the wavelength in percent:'};
dlgtitle = 'Planck''s Law Inputs';
fieldsize = [1 45; 1 45; 1 45; 1 45; 1 45; 1 45];
definput = {'3000', '100', '10', '1', '[1, 2500]', '0.5'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

% Converting inputs to numbers
T = str2num(answer{1});
npoints = str2double(answer{2});
noise_scale = str2double(answer{3});
sigma = str2double(answer{4});
wavelengths = str2num(answer{5});
lambda_err_percent = str2num(answer{6});


% Radiated power density plancks law via wavelength
I = @(lambda) (2*pi*c^2*h./(lambda.^5)).*(exp(h*c./(lambda*k_b.*T)) - 1).^-1;

% Gaussian function to convolve
gaussian = @(x) exp(-x.^2/(2*sigma^2));

% Generating an array of evenly spaced wavelengths
lambda = linspace(min(wavelengths), max(wavelengths), npoints); %  in nanometers
lambda_err = lambda + (lambda_err_percent/100) .* lambda .* randn(size(lambda));
lambda_m = lambda*nm;

% Calcualting power densities
power_density = I(lambda_m);
drift_power_density = I(lambda_err*nm);

% Adding normally distributed noise
noisy_signal = rand(size(power_density));
drift_noisy_signal = rand(size(drift_power_density));

% Creating an amplitude for the noise that is adjusted by noise_scale
amplitude = (noise_scale/100)*power_density;
drift_amplitude = (noise_scale/100)*drift_power_density;

% Adding the noise only signal to the clean signal
noisy_power_density = power_density + amplitude.*noisy_signal;
drift_noisy_power_density = drift_power_density + drift_amplitude.*drift_noisy_signal;

% Line broadening parameters
x = linspace(-5*sigma, 5*sigma, npoints);
gaussian_kernel = gaussian(x);
gaussian_kernel = gaussian_kernel/sum(gaussian_kernel); % Normalize

% Apply line broadening to ideal and noisy signal
broad_power_density = conv(power_density, gaussian_kernel, 'same');
broad_noisy_power_density = broad_power_density + amplitude.*noisy_signal;

broad_drift_power_density = conv(drift_noisy_power_density, gaussian_kernel,'same');
all_power_density = broad_drift_power_density + drift_amplitude.*drift_noisy_signal;

% Calculating Wien's constant for all configurations
ideal_b = wienconst(lambda_m, power_density, T)
noisy_b = wienconst(lambda_m, noisy_power_density, T)
broad_b = wienconst(lambda_m, broad_power_density, T)
broad_noisy_b = wienconst(lambda_m, broad_noisy_power_density, T)
drift_b = wienconst(lambda_err*nm, all_power_density, T)

% Calculating the Stefan-Boltzmann constant for all configurations
ideal_sb_const = sbconst(lambda_m, power_density, T)
noisy_sb_const = sbconst(lambda_m, noisy_power_density, T)
broad_sb_const = sbconst(lambda_m, broad_power_density, T)
broad_noisy_sb_const = sbconst(lambda_m, broad_noisy_power_density, T)
drift_sb_const = sbconst(lambda_err*nm, all_power_density, T)

% Ideal theoretical ouput
figure (1)
hold on
grid on
plot(lambda, power_density/10^13)
title('Theoretical Output')
xlabel('\lambda (nm)')
ylabel('Power Density (10^{13} watts/m^{3})')

% Noisy output
figure (2)
hold on
grid on
plot(lambda, noisy_power_density/10^13)
title('Noisy Output')
subtitle(['Noise level: ',num2str(noise_scale), '%'])
xlabel('\lambda (nm)')
ylabel('Power Density (10^{13} watts/m^{3})')

% Broadened theoretical output
figure (3)
hold on
grid on
plot(lambda, broad_power_density/10^13)
title('Broadened Output')
subtitle(['LB = ',num2str(sigma)])
xlabel('\lambda (nm)')
ylabel('Power Density (10^{13} watts/m^{3})')

figure (4)
hold on
grid on
plot(lambda, broad_noisy_power_density/10^13)
title('Broadened and Noisy Output')
subtitle(['Noise level: ', num2str(noise_scale), '%', ', LB = ', num2str(sigma)])
xlabel('\lambda (nm)')
ylabel('Power Density (10^{13} watts/m^{3})')

figure (5)
hold on
grid on
plot(lambda_err, all_power_density/10^13)
title(['\Delta\lambda: ',num2str(lambda_err_percent),'%'])
xlabel('\lambda (nm)')
ylabel('Power Density (10^{13} watts/m^{3})')


function [b_const] = wienconst(lambda, power_density, T)

    lambda_max = lambda(power_density == max(power_density))
    b_const = lambda_max.*T;

end

function [sb_const] = sbconst(lambda, power_density, T)

    flux = trapz(lambda, power_density);
    sb_const = flux./T.^4;

end