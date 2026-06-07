% This code demonstrates the generation of highly multimode soliton.

close all; clearvars;

%% Load folders
addpath('../../GMMNLSE algorithm/','../../user_helpers/');

% Add multimode fiber informaiton
fiber.MM_folder = '../../../BuildFiber/Fibers/step_wavelength1550nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_50modes.mat';

%% Setup fiber parameters
sim.lambda0 = 1550e-9; % the central wavelength
sim.pulse_centering = false;
sim.midx = 1:25; % use two spatial modes for fast demonstration; users are free to try more
sim.gpu_yes = true; % make it "false" to test CPU
num_modes = length(sim.midx);

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
sim.include_Raman = true;
sim.gain_model = 0;
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode'); % load default parameters

num_save = 100;
fiber.L0 = 10; % m
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^12; % the number of time points
time_window = 50; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
total_energy = 200; % nJ
tfwhm = 0.5; % ps
% input_field = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), Nt, {'ifft',0}, ones(1,length(sim.midx)), -time_window/2*0.4);
coeffs = zeros(1, num_modes);
coeffs(16:25) = 1;    % equal amplitude → equal energy (20 nJ each) in modes 16–25

input_field = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, Nt, ...
                                 {'ifft',0}, coeffs, -time_window/2*0.4);
input_field.fields = input_field.fields.*exp(1i*2*pi*num_modes);

%% Propagate
tic;
prop_output = GMMNLSE_propagate(fiber,input_field,sim);
sim_time = toc;
fprintf('Simulation time: %.2f s (%.2f min)\n', sim_time, sim_time/60);

%% Save prop_output
save_dir = fileparts(mfilename('fullpath'));
save(fullfile(save_dir, sprintf('prop_output_%dmodes.mat', num_modes)), 'prop_output', '-v7.3');

%% Plot
% Time
fig1 = figure;
plot(t,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',14);
savefig(fig1, fullfile(save_dir, sprintf('time_%dmodes.fig', num_modes)));
saveas(fig1, fullfile(save_dir, sprintf('time_%dmodes.png', num_modes)));

% Spectrum
fig2 = figure;
plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2*c./lambda.^2,'linewidth',2);
xlabel('\lambda (nm)');
ylabel('PSD (a.u.)');
% xlim([1500,1800]);
set(gca,'fontsize',14);
savefig(fig2, fullfile(save_dir, sprintf('spectrum_%dmodes.fig', num_modes)));
saveas(fig2, fullfile(save_dir, sprintf('spectrum_%dmodes.png', num_modes)));

% Comparison of time
fig3 = figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('Time (ps)');
ylabel('Propagation distance (m)');
set(gca,'fontsize',14);
savefig(fig3, fullfile(save_dir, sprintf('time_evolution_%dmodes.fig', num_modes)));
saveas(fig3, fullfile(save_dir, sprintf('time_evolution_%dmodes.png', num_modes)));

% Comparison of spectra
fig4 = figure;
[x,y] = meshgrid((f-sim.f0),prop_output.z(2:end));
tmp = 10*log10(permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2])); tmp = tmp - max(tmp(:));
pcolor(x,y,tmp);
shading interp; colormap(jet); caxis([-20,0]);
xlabel('\Deltaf (THz)');
ylabel('Propagation distance (m)');
set(gca,'fontsize',14);
savefig(fig4, fullfile(save_dir, sprintf('spectrum_evolution_%dmodes.fig', num_modes)));
saveas(fig4, fullfile(save_dir, sprintf('spectrum_evolution_%dmodes.png', num_modes)));