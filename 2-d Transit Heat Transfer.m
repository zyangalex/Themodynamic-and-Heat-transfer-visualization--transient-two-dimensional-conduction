% Program:
% Transient_Conduction_With_Finite_Differences.m
% Transient 2D Conduction Solver using Finite Difference Method.
%
% Description:
% Numerically solves the transient two dimensional conduction problem
% using the finite difference method and plots color contour plot. Assumes 
% transient 2D conduction with constant properties.
% 
% Variable List:
% T = Temperature (deg. Celsius)
% T1 = Boundary condition temperature 1 (deg. Celsius)
% T2 = Boundary condition temperature 2 (deg. Celsius)
% theta = Non-dimensionalized temperature difference = (T-T1)/(T2-T1)
% Lx = Plate length in x-direction (m)
% Ly = Plate length in y-direction (m)
% x = Create x-distance node locations
% y = Create y-distance node locations
% Nx = Number of increments in x-direction
% Ny = Number of increments in y-direction
% dx = Increment size in x-direction (m)
% dy = Increment size in y-direction (m)
% dT = Temperature step between contours
% Lmax = Maximum number of time steps before stopping
% Told = Stores temperature array for previous time step
% diff = Percent difference at x = y = 0.4 m compared to theoretical
% k = Thermal conductivity (W/mK)
% rho = Density (kg/m^3)
% Cp = Specific heat (J/kgK)
% alpha = Thermal diffusivity (m^2/s)
% deltaT = Time step (s)
% deltaT_stable_max = Maximum stable time step (s)
% Fo = Fourier number 
% tol = Stop when T at x = y = 0.4 m reaches within this percent of tol_ans
% tol_ans = Steady-state temperature at x = y = 0.4m (deg. Celsius)
% Tplot = Stores T (deg. Celsius) at x = y = 0.4 m at each time step
% L = Loop counter
% p = Current iteration
% i = Current column
% j = Current row
% v = Sets temperature levels for contours
% Nc = Number of contours for plot

clear, clc              % Clear command window and workspace

Lx = 1;                 % Plate length in x-direction (m)
Ly = 1;                 % Plate length in y-direction (m)
Nx = 20;                % Number of increments in x-direction

Ny = Nx;                % Number of increments in y-direction
dx = Lx/Nx;             % Increment size in x-direction (m)
dy = Ly/Ny;             % Increment size in y-direction (m)

T1 = 0;                 % BC temperature 1 (deg. Celsius)
T2 = 100;               % BC temperature 2 (deg. Celsius)

k = 222;               % Thermal conductivity (W/mK)
rho = 2800;             % Density (kg/m^3)
Cp = 896;               % Specific heat (J/kgK)
alpha = k/(rho*Cp);     % Thermal diffusivity (m^2/s)

deltaT = 5;                          % Time step (seconds)
deltaT_stable_max = .25*dx^2/alpha;     % Maximum stable time step
Fo = alpha*deltaT/dx^2;                 % Fourier number 

tol = 0.1;              % Stop when T at x = y = 0.4m reaches percentage of tol_ans
tol_ans = 16.8568;      % Steady-state temperature at x = y = 0.4m
Lmax = 10^4;            % Maximum number of time steps before stopping

T = T1*ones(Nx+1,Ny+1); % Initialize T array to T1 everywhere
T(2:Nx,Ny+1) = T2;      % Initialize top row to T2 boundary condition
T(1,Ny+1)    = T1;      % Initialize top left
T(Nx+1,Ny+1) = T1;      % Initialize top right
Tplot = ones(Lmax,1);   % Initialize Tplot to allocate memory

x = 0:dx:Lx;            % Create x-distance node locations
y = 0:dy:Ly;            % Create y-distance node locations

for L = 1:Lmax          % Loop through time steps
    Told = T;           % Store previous T array as Told for next time step
    
    for j = 2:Ny        % Loop through rows
        for i = 2:Nx    % Loop through columns
            
            % Calculates temperatures for new time step
            T(i,j) = Fo*(Told(i-1,j) + Told(i+1,j) + Told(i,j-1)...
                + Told(i,j+1)) + (1-4*Fo)*Told(i,j);
        end 
    end    
    
    Tplot(L) = T(Nx/2-1,Ny/2-1); % Store T at x = y = 0.4m
    
    % Percent difference at x = y = 0.4 m compared to theoretical
    diff = abs((T(Nx/2-1,Ny/2-1) - tol_ans)/tol_ans*100); ...
        
    fprintf('Time step = %8.0f - Diff. = %10.6f percent\n', L, diff);
    if (diff < tol)     % Exit iteration loop because reached steady-state
        break 
    end
        
    Nc = 50;                    % Number of contours for plot
    dT = (T2 - T1)/Nc;          % Temperature step between contours
    v = T1:dT:T2;               % Sets temperature levels for contours
    colormap(jet)               % Sets colors used for contour plot
    contourf(x, y, T',v, 'LineStyle', 'none') 
    colorbar                    % Adds a scale to the plot
    axis equal tight            % Makes the axes have equal length
    title(['Contour Plot of Temperature in deg. C at time = ',...
        num2str(deltaT*L/3600),' h']) 
    xlabel('x (m)') 
    ylabel('y (m)')        
    set(gca,'XTick',0:.1:Lx)    % Sets the x-axis tick mark locations
    set(gca,'YTick',0:.1:Ly)    % Sets the y-axis tick mark locations
    pause(0.001)                % Pause between time steps to display graph
        
    %if L == 55 || L == 65 || L == 80 % Chosen time steps to save plot
    %    saveas(gcf, ['Transient_Plot_Unstable_',num2str(L)], 'jpg'); % save plot
    %end    
end

fprintf('Number of time steps = \t %8.0f \n\n', L) % Print time steps

if (L == Lmax)      % Warn if number of iterations exceeds maximum
    disp('Warning: Maximum time steps exceeded')   
    fprintf('\n') 
end

figure

plot(Tplot(1:L))
xlabel('Timestep')
ylabel('Temperature (C)')
    