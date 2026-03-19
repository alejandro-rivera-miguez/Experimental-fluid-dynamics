%% =======================================================================
% SCRIPT: Lab2_NACA23012_FullAnalysis
% DESCRIPTION: Aerodynamic characterization of a NACA 23012 airfoil using 
%              surface pressure measurements and hot-wire wake data.
%              Calculates lift (Cl) and pressure drag (Cdp) via surface 
%              pressure integration, and estimates total drag (Cd) by 
%              adding the viscous contribution from the wake velocity deficit.
%
% AUTHORS: Alejandro Rivera Míguez
%          Alberto Rivero García
%          Mikel Segovia Díaz
% COURSE: Experimental Fluid Dynamics - Politecnico di Milano
%% =======================================================================

clear all; close all; clc;

addpath("Measures\")
%% 1. NACA 23012 AIRFOIL GEOMETRY GENERATION
% NACA 23012 Airfoil Parameters
x_m = 0.2025;          % Position of maximum camber (15% of the chord)
k1  = 15.957;          % Constant for design lift coefficient (Cl=0.3)
t   = 0.12;            % Maximum thickness (12% of the chord)
c   = 0.1;             % Chord length [m]

% Chord discretization
x = linspace(0, c, 200); 

% Array initialization
yc = zeros(size(x));
dyc_dx = zeros(size(x));

% Camber line and slope calculation
for i = 1:length(x)
    if x(i)/c < x_m
        yc(i) = c * ((k1 / 6) * ((x(i)/c)^3 - 3 * x_m * (x(i)/c)^2 + x_m^2 * (3 - x_m) * x(i)/c));
        dyc_dx(i) = ((k1 / 6) * (3 * (x(i)/c)^2 - 6 * x_m * x(i)/c + x_m^2 * (3 - x_m)));
    else
        yc(i) = c * (k1 * x_m^3 / 6) * (1 - x(i)/c);
        dyc_dx(i) = -(k1 * x_m^3 / 6);
    end
end

% Thickness distribution (y_t)
yt = c * 5 * t * (0.2969 * sqrt(x/c) - 0.1260 * x/c - 0.3516 * (x/c).^2 + 0.2843 * (x/c).^3 - 0.1015 * (x/c).^4);

%% 2. EXACT SURFACE DERIVATIVES (SYMBOLIC MATH)
syms x_espesor t_espesor
yt_espesor = c * 5 * t_espesor * (0.2969 * sqrt(x_espesor/c) - 0.1260 * x_espesor/c - 0.3516 * (x_espesor/c)^2 + 0.2843 * (x_espesor/c)^3 - 0.1015 * (x_espesor/c)^4);
dyt_dx = diff(yt_espesor, x_espesor);

disp('--- Symbolic Derivative of thickness (yt) with respect to x: ---');
disp(dyt_dx);

% Camber slope angle (theta) and surface coordinates
theta = atan(dyc_dx);
xu = x - yt .* sin(theta);  
yu = yc + yt .* cos(theta); 
xl = x + yt .* sin(theta);  
yl = yc - yt .* cos(theta); 

%% 3. PRESSURE TAPS DEFINITION
pex = c * [0.022, 0.031, 0.062, 0.155, 0.201, 0.33, 0.53, 0.681];
pez = c * [0.034, 0.04, 0.053, 0.071, 0.075, 0.077, 0.063, 0.047];
pix = c * [0.01, 0.028, 0.058, 0.147, 0.196, 0.326, 0.427, 0.837];
piz = c * [-0.01, -0.017, -0.023, -0.033, -0.038, -0.044, -0.043, -0.013];

%% 4. NORMAL VECTORS & SURFACE INTEGRATION PREPARATION
derivada_yc_x_ext = zeros(8,1); derivada_yc_x_int = zeros(8,1);
derivada_yt_x_ext = zeros(8,1); derivada_yt_x_int = zeros(8,1);

for i = 1:8
    if pex(i)/c < x_m
        derivada_yc_x_ext(i) = ((k1 / 6) * (3 * (pex(i)/c)^2 - 6 * x_m * pex(i)/c + x_m^2 * (3 - x_m)));
    else
        derivada_yc_x_ext(i) = - (k1 * x_m^3 / 6);
    end
    if pix(i)/c < x_m
        derivada_yc_x_int(i) = ((k1 / 6) * (3 * (pix(i)/c)^2 - 6 * x_m * pix(i)/c + x_m^2 * (3 - x_m)));
    else
        derivada_yc_x_int(i) = - (k1 * x_m^3 / 6);
    end
    derivada_yt_x_ext(i) = subs(dyt_dx, [x_espesor, t_espesor], [pex(i), 0.12]);
    derivada_yt_x_int(i) = subs(dyt_dx, [x_espesor, t_espesor], [pix(i), 0.12]);
end

parcialz_ext = derivada_yc_x_ext + derivada_yt_x_ext;
parcialz_int = derivada_yc_x_int - derivada_yt_x_int;

nx_ext = -parcialz_ext ./ sqrt(1 + parcialz_ext.^2);
nz_ext =  1 ./ sqrt(1 + parcialz_ext.^2);
nx_int =  parcialz_int ./ sqrt(1 + parcialz_int.^2);
nz_int = -1 ./ sqrt(1 + parcialz_int.^2);

%% 5. AERODYNAMIC FORCES CALCULATION (PRESSURE INTEGRATION)
% Load experimental pressure data (Add 'ReadVariableNames', false to avoid warnings if needed)
presiones = readtable('Measures/Pressure/gruppo5_2.dat', 'ReadVariableNames', false);
presiones_matriz = table2array(presiones);

alpha = -12:2:20;
indx_alpha_8 = find(alpha == 8);

Fz_ext = zeros(size(alpha)); Fz_int = zeros(size(alpha)); CFz = zeros(size(alpha));
Fx_ext = zeros(size(alpha)); Fx_int = zeros(size(alpha)); CFx = zeros(size(alpha));
Cl     = zeros(size(alpha)); Cdp    = zeros(size(alpha)); q   = zeros(size(alpha));

load("tunnelcoefficient.mat", "tunnelcoeff");

for j = 1:length(alpha)
    q(j) = tunnelcoeff * presiones_matriz(j, 18);
    
    Fz_ext(j) = trapz(pex, -presiones_matriz(j, 1:8) .* nz_ext');
    Fz_int(j) = trapz(pix, -presiones_matriz(j, 9:16) .* nz_int');
    CFz(j)    = (Fz_ext(j) + Fz_int(j)) / q(j) / c;
    
    Fx_ext(j) = trapz(pex, -presiones_matriz(j, 1:8) .* nx_ext');
    Fx_int(j) = trapz(pix, -presiones_matriz(j, 9:16) .* nx_int');
    CFx(j)    = (Fx_ext(j) + Fx_int(j)) / q(j) / c;
    
    Cl(j)  = CFz(j) * cosd(alpha(j)) - CFx(j) * sind(alpha(j));
    Cdp(j) = CFx(j) * cosd(alpha(j)) + CFz(j) * sind(alpha(j));
end

%% 6. WAKE VELOCITY PROFILE & VISCOUS DRAG CALCULATION (HOT-WIRE)
% Dynamic File Reading (Eliminates warnings and repetitive code)
Z_mm = [100, 119, 124, 126, 128, 130, 133, 137, 142, 146, 151, 152, 155, 170];
Z_filenames = Z_mm / 10; % Generates 10.0, 11.9, etc.
U_wake = zeros(size(Z_mm));

for i = 1:length(Z_mm)
    % Construct filename (Change path if files are inside a specific folder)
    % e.g., fname = sprintf('GRUPPO5/velocita/gruppo5_1_Z%.1f.dat', Z_filenames(i));
    fname = sprintf('Measures/Velocity/gruppo5_1_Z%.1f.dat', Z_filenames(i));
    
    % Read table suppressing header warnings, convert to array, and compute mean
    temp_table = readtable(fname, 'ReadVariableNames', false);
    temp_array = table2array(temp_table);
    U_wake(i) = mean(temp_array(:)); 
end

% Drag coefficient (Cd) calculation using Momentum Deficit
c_mm = 100; % Reference length in mm
Uinf = (U_wake(1) + U_wake(end)) / 2; % Free-stream velocity approximation
Cd_alpha8 = 2 * trapz(Z_mm / c_mm, 1 - (U_wake / Uinf).^2);

fprintf('\n--- HOT-WIRE WAKE ANALYSIS ---\n');
fprintf('The total drag coefficient (Cd) at alpha = 8 deg is: %.4f\n', Cd_alpha8);

% Viscous Drag Correction
Cdv = Cd_alpha8 - Cdp(indx_alpha_8); % Viscous contribution
Cd = Cdp + Cdv; % Total corrected drag

%% 7. AERODYNAMIC EFFICIENCY (L/D) CALCULATION
LD_ratio = abs(Cl ./ Cd);
[LD_max, idx_max_LD] = max(LD_ratio);
alpha_opt = alpha(idx_max_LD); 

fprintf('\n--- EFFICIENCY ANALYSIS ---\n');
fprintf('Maximum L/D = %.2f occurs at alpha = %.1f degrees.\n\n', LD_max, alpha_opt);

%% 8. PLOTTING RESULTS
% --- Plot 1: Aerodynamic Profile ---
figure('Name', 'Airfoil Geometry');
plot(xu, yu, 'b-', 'LineWidth', 2, 'DisplayName', 'Upper surface'); hold on;
plot(xl, yl, 'r-', 'LineWidth', 2, 'DisplayName', 'Lower surface');
plot(x, yc, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Camber line');
scatter(pex, pez, 70, 'b', 'filled', 'DisplayName', 'Pressure taps (dorsal)');
scatter(pix, piz, 70, 'r', 'filled', 'DisplayName', 'Pressure taps (ventral)');
axis equal; grid on;
xlabel('$x$ (Chord)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$y$ (Height)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('NACA 23012 Profile Geometry', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend('Location', 'bestoutside', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% --- Plot 2: Wake Velocity Profile ---
figure('Name', 'Wake Velocity Profile');
plot(U_wake, Z_mm, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.85 0.325 0.098])
grid on; grid minor;
xlabel('Mean velocity $U$ [m/s]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Height $Z$ [mm]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Wake Velocity Profile (Hot-Wire Anemometry)', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% --- Plot 3: Lift and Total Drag Coefficients ---
figure('Name', 'Aerodynamic Coefficients');
hold on;
plot(alpha, Cl, 'b-', 'LineWidth', 2, 'DisplayName', '$C_L$');
plot(alpha, Cd, 'r-', 'LineWidth', 2, 'DisplayName', '$C_{D}$ (Total)');
plot(alpha, Cdp, 'r--', 'LineWidth', 1.5, 'DisplayName', '$C_{D_p}$ (Pressure)');
xlabel('$\alpha$ (Angle of attack)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Aerodynamic coefficients', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Lift and Drag Coefficients vs $\alpha$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% --- Plot 4: Polar Plot (Cl vs Total Cd) ---
figure('Name', 'Aerodynamic Polar');
plot(Cd, Cl, 'k-', 'LineWidth', 2); hold on;
scatter(Cd, Cl, 60, 'b', 'filled'); 
xlabel('$C_{D}$ (Total drag coefficient)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$C_L$ (Lift coefficient)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Aerodynamic polar ($C_L$ vs $C_{D}$)', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% --- Plot 5: Lift-to-Drag Ratio vs Angle of Attack ---
figure('Name', 'Efficiency L/D');
plot(alpha, LD_ratio, 'k-', 'LineWidth', 2, 'DisplayName', '$L/D$'); hold on;
scatter(alpha_opt, LD_max, 70, 'r', 'filled', 'DisplayName', sprintf('Max L/D at $\\alpha = %.1f^\\circ$', alpha_opt));
xlabel('$\alpha$ (Angle of attack)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$L/D$ (Efficiency)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Lift-to-Drag Ratio vs Angle of Attack', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');