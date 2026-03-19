%% =======================================================================
% SCRIPT: Lab1_Force_Balance_Analysis
% DESCRIPTION: Wind tunnel calibration, external force balance calibration, 
%              and aerodynamic characterization of an airfoil. 
%              1. Computes the wind tunnel coefficient via Linear Regression.
%              2. Computes the balance calibration matrix via Least Squares.
%              3. Extracts the aerodynamic coefficients (CL, CD, CM) 
%                 across varying angles of attack.
%
% AUTHORS: Alejandro Rivera Míguez   
%          Alberto Rivero García
%          Mikel Segovia Díaz
% COURSE: Experimental Fluid Dynamics - Politecnico di Milano
%% =======================================================================

clear all; close all; clc;

addpath("Measures\")
addpath("Auxiliary functions\")

%% 1. WIND TUNNEL CALIBRATION
fprintf('--- 1. Computing Wind Tunnel Calibration Coefficient ---\n');

% Import calibration data (Ensure 'taratura_galleria.dat' is in the path)
TunnelData = readtable("gallery_calibration.dat");

St = 200; % Specific tunnel parameter [Pa/V]

% Calculate pressures in Pascals
P1mP2 = TunnelData.Var6 * St; 
DynP  = (TunnelData.Var5 - TunnelData.Var4) * St; 

% Perform Linear Regression (forced through origin: y = m*x)
tunnelcoeff = P1mP2 \ DynP;

% R^2 calculation to assess the regression quality
DynPmean = mean(DynP);
DynPpredict = tunnelcoeff * P1mP2;
R2 = 1 - sum((DynP - DynPpredict).^2) / sum((DynP - DynPmean).^2);

% Print results to console
fprintf('Tunnel Calibration Coefficient: %.4f\n', tunnelcoeff);
fprintf('Regression R-squared (R^2): %.4f\n\n', R2);

% Save coefficient for future use
save("Measures/tunnelcoefficient", "tunnelcoeff");

% --- Plot 1: Linear Regression for Wind Tunnel Calibration ---
figure('Name', 'Wind Tunnel Calibration', 'Color', 'w');
plot(P1mP2, DynP, 'bo', 'LineWidth', 1.5, 'DisplayName', 'Experimental Data'); hold on;
plot(P1mP2, DynPpredict, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Linear Fit ($R^2 = %.4f$)', R2));
grid on; grid minor;
xlabel('$\Delta P_{1-2}$ [Pa]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Dynamic Pressure $q$ [Pa]', 'Interpreter', 'latex', 'FontSize', 12);
title('Wind Tunnel Calibration Regression', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 11);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');


%% 2. BALANCE CALIBRATION MATRIX COMPUTATION
fprintf('--- 2. Computing Balance Calibration Matrix (Least Squares) ---\n');

% Import balance calibration raw data
ExpData = importfile("Leveling.dat", [2 12]);
ExpMat = zeros([11, 3]);
ExpMat(:,1) = ExpData.NI9237ch0VVch1VVch2VVNI9205Ch0Vch1Vch2V;
ExpMat(:,2) = ExpData.VarName2;
ExpMat(:,3) = ExpData.VarName3;

% Fx real definition
FxReal = zeros([11, 1]);
FxReal(1:3)   = 9.81 * [1, 0.5, 1.5]';
FxReal(10:11) = sqrt(2)/2 * 9.81 * [1, 0.5]';

% Fy real definition
FyReal = zeros([11, 1]);
FyReal(4:9)   = 9.81 * [1, 0.5, 1.5, 1, 0.5, 1.5]';
FyReal(10:11) = FxReal(10:11);

% M real definition
Mreal = zeros([11, 1]);
Mreal(7:9) = -0.2 * 9.81 * [1, 0.5, 1.5]';

% Apply Least Squares Method to find calibration coefficients
Row1 = lsqr(ExpMat, FxReal);
Row2 = lsqr(ExpMat, FyReal);
Row3 = lsqr(ExpMat, Mreal);

% Assemble Calibration Matrix
K = [Row1'; Row2'; Row3'];

% Test calibration matrix residues
norm_residue = ones([11, 1]);
for j = 1:11
    V_residue = [FxReal(j); FyReal(j); Mreal(j)] - K * [ExpMat(j,1); ExpMat(j,2); ExpMat(j,3)];
    norm_residue(j) = norm(V_residue);
end

% Save Calibration Matrix for external references if needed
save("Measures/CalMatrix", "K");
fprintf('Calibration Matrix K computed and saved successfully.\n\n');


%% 3. AERODYNAMIC FORCES CALCULATION
fprintf('--- 3. Processing Wind Tunnel Aerodynamic Data ---\n');

% Import experimental force data from the airfoil test
ForceExp = importfileforce("Measures/gruppo_5.dat", [2 16]);

% Parameter definitions
c = 0.1;       % Chord length [m]
b = 0.294;     % Wingspan [m]
S = c * b;     % Wing area [m^2]
rho = 1.195;   % Air density [kg/m^3]
mu = 1.8e-5;   % Dynamic viscosity of air [Pa·s]

% Angle of attack vector
alphav = -12:2:10;
alphav = [alphav, 11, 12, 13];

% Variable initialization
L = zeros([15, 1]);       % Lift
D = zeros([15, 1]);       % Drag
M = zeros([15, 1]);       % Moment
DynP = zeros([15, 1]);    % Dynamic Pressure

% Force and dynamic pressure calculation loop
for j = 1:15
    % Convert raw voltages to physical forces using the Calibration Matrix
    Forces = K * [ForceExp.Fx(j); ForceExp.Fy(j); ForceExp.M(j)];
    
    % Projection into Wind Axes
    D(j) = -(Forces(1) * cosd(alphav(j)) + Forces(2) * sind(alphav(j)));
    L(j) = -(Forces(2) * cosd(alphav(j)) - Forces(1) * sind(alphav(j)));
    M(j) = Forces(3);
    
    % Dynamic pressure computation (using the computed tunnelcoeff)
    DynP(j) = tunnelcoeff * St * ForceExp.P1mP2(j);
end

% Dimensionless coefficients
CD = D ./ DynP / S;
CL = L ./ DynP / S;
CM = M ./ DynP / S / c;

% Experimental velocity and Reynolds number calculation
Vexp = sqrt(2 * DynP / rho);
Re = rho * Vexp * c / mu;

fprintf('Aerodynamic coefficients computed successfully. Generating plots...\n');


%% 4. PLOTTING RESULTS
% --- Plot 2: Forces as a function of the angle of attack ---
figure('Name', 'Absolute Aerodynamic Forces', 'Color', 'w');
plot(alphav, L, '-ob', 'LineWidth', 1.5, 'DisplayName', 'Lift ($L$)')
hold on
plot(alphav, D, '-sr', 'LineWidth', 1.5, 'DisplayName', 'Drag ($D$)')
plot(alphav, M, '-^g', 'LineWidth', 1.5, 'DisplayName', 'Moment ($M$)')
hold off
grid on; grid minor;
xlabel('Angle of Attack ($\alpha$) [$^\circ$]', 'FontSize', 12, 'Interpreter', 'latex')
ylabel('Forces [N]', 'FontSize', 12, 'Interpreter', 'latex')
title('Forces (Lift, Drag, Moment) vs. Angle of Attack', 'FontSize', 14, 'Interpreter', 'latex')
legend('Location', 'best', 'FontSize', 11, 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

% --- Plot 3: Dimensionless coefficients as a function of angle of attack ---
figure('Name', 'Aerodynamic Coefficients', 'Color', 'w');
plot(alphav, CL, '-ob', 'LineWidth', 1.5, 'DisplayName', 'Lift Coefficient ($C_L$)')
hold on
plot(alphav, CD, '-sr', 'LineWidth', 1.5, 'DisplayName', 'Drag Coefficient ($C_D$)')
plot(alphav, CM, '-^g', 'LineWidth', 1.5, 'DisplayName', 'Moment Coefficient ($C_M$)')
hold off
grid on; grid minor;
xlabel('Angle of Attack ($\alpha$) [$^\circ$]', 'FontSize', 12, 'Interpreter', 'latex')
ylabel('Dimensionless Coefficients', 'FontSize', 12, 'Interpreter', 'latex')
title('Coefficients ($C_L$, $C_D$, $C_M$) vs. Angle of Attack', 'FontSize', 14, 'Interpreter', 'latex')
legend('Location', 'best', 'FontSize', 11, 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')

% --- Plot 4: Coefficient polar plot (C_D vs. C_L) ---
figure('Name', 'Aerodynamic Polar', 'Color', 'w');
plot(CD, CL, '-db', 'LineWidth', 1.5)
grid on; grid minor;
xlabel('Drag Coefficient ($C_D$)', 'FontSize', 12, 'Interpreter', 'latex')
ylabel('Lift Coefficient ($C_L$)', 'FontSize', 12, 'Interpreter', 'latex')
title('Aerodynamic Polar ($C_L$ vs. $C_D$)', 'FontSize', 14, 'Interpreter', 'latex')
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')