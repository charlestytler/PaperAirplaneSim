%  Paper Airplane Simulation
%  Based on simulation obtained from:
%  http://www.princeton.edu/~stengel/PaperPlane.html
%	Example 1.3-1 Paper Airplane Flight Path
%	Copyright 2005 by Robert Stengel
%	August 23, 2005


%% User-set Parameters

    % Initial Airspeed (m/s)
    V0 = 15.0;
    % Initial Throw Angle (rad) (0 = Level)
    Gam0 = deg2rad(10);
    % Initial Height (meters)
    H0 = 2;
    
    % Option to plot all numerical solution methods (ODE, RK4, EULER)
    plotAllMethods = 1;  %Set to 1 for TRUE, 0 for FALSE
    
    % Time step for RK4 and EULER fixed-step solution methods
    dt = 0.01; %(sec)


%%   Paper Airplane configuration paramters
    global CL CD S m g rho
	S		=	0.017;			% Reference Area, m^2
	AR		=	0.86;			% Wing Aspect Ratio
	e		=	0.9;			% Oswald Efficiency Factor;
	m		=	0.003;			% Mass, kg
	g		=	9.81;			% Gravitational acceleration, m/s^2
	rho		=	1.225;			% Air density at Sea Level, kg/m^3	
	
%%   Aerodynamic Coefficients
    CLa		=	pi * AR/(1 + sqrt(1 + (AR / 2)^2));
							% Lift-Coefficient Slope, per rad
	CDo		=	0.02;			% Zero-Lift Drag Coefficient
	epsilon	=	1 / (pi * e * AR);      % Induced Drag Factor	
	CL		=	sqrt(CDo / epsilon);	% CL for Maximum Lift/Drag Ratio
	CD		=	CDo + epsilon * CL^2;	% Corresponding CD
	
    
%%	a) Equilibrium conditions for Glide at Maximum Lift/Drag Ratio
    LDmax	=	CL / CD;			% Maximum Lift/Drag Ratio
	Gam_eq		=	-atan(1 / LDmax);	% Corresponding Flight Path Angle, deg
	V_eq	=	sqrt(2 * m * g /(rho * S * (CL * cos(Gam_eq) - CD * sin(Gam_eq))));
							% Corresponding Velocity, m/s
	Alpha	=	CL / CLa;			% Corresponding Angle of Attack, rad
	H_eq	=	2;			% Initial Height, m
	R_eq	=	0;			% Initial Range, m
	
    
    % Equilibrium solutions:
    %Gam_eq = -0.1794; V_eq = 3.5512; H_eq = 2; R_eq = 0
    
    
    
%%   Initial Conditions

    % Initial Range (meters)
    R0 = 0;
    
    t0		=	0;			% Initial Time, sec
	tf		=	6;			% Final Time, sec
	tspan	=	[t0 tf];
	x0 = [V0; Gam0; H0; R0];
    
    
%%   Time History solution of Equations of Motion
    
	% Method 1: ODE Solver (variable time step, low order)
    [tODE,xODE]	=	ode23(@EqMotion,tspan,x0);
    
    % Method 2: Fixed step 4th Order Runge Kutta
    % Step size set in user parameters
    tRK = [t0:dt:tf];
    xRK = zeros(length(tRK),4);
    xRK(1,:) = x0';
    for k = 2:length(tRK)
        xRK(k,:) = RK4(xRK(k-1,:)',dt)';
    end
    
    % Method 3: Fixed step Euler's Method (1st Order)
    % Step size set in user parameters
    tEM = [t0:dt:tf];
    xEM = zeros(length(tEM),4);
    xEM(1,:) = x0';
    for k = 2:length(tRK)
        xdot = EqMotion(-1,xEM(k-1,:)')';
        xEM(k,:) = xEM(k-1,:) + xdot*dt;
    end
	
    
%%   Plot time histories
    
    plotTimeHist(tODE,xODE,tRK,xRK,tEM,xEM,plotAllMethods);

    

%% Subfunctions

	function xdot = EqMotion(t,x)
    %	Fourth-Order Equations of Aircraft Motion
    
    % Get necessary global config parameters
    global CL CD S m g rho
	
    % Extract state values
	V       =	x(1);  % Airspeed
	Gamma   =	x(2);  % Flight Path Angle
    h       =   x(3);  % Height
    r       =   x(4);  % Range
    
    % Calculate dynamic pressure for lift/drag expressions
	q	=	0.5 * rho * V^2;	% Dynamic Pressure, N/m^2
    
    % Calculate forces
    L = CL * q * S;
    D = CD * q * S;
	
    % Equations of Motion update
	xdot	=	[(-D - m*g*sin(Gamma)) / m
				 (L - m*g*cos(Gamma)) / (m * V)
				 V * sin(Gamma)
				 V * cos(Gamma)];
             
    end

    
    function x_next = RK4(x,dt)
    % 4th Order Runge Kutta Calculation
    % Returns estimate of x(t + dt) given x(t) and u(t)
    
    t = -1; %Meaningless number, not used by EqMotion
    
    K1 = EqMotion(t,x);
    K2 = EqMotion(t,x + K1*dt/2);
    K3 = EqMotion(t,x + K2*dt/2);
    K4 = EqMotion(t,x + K3*dt);
    x_next = x + 1/6*(K1 + 2*K2 + 2*K3 + K4)*dt;
    
    end

    
    function plotTimeHist(tODE,xODE,tRK,xRK,tEM,xEM,plotAllMethods)
    
    % Create large figure;
	fig = figure;
    figPos = fig.Position;
    fig.Position(3) = 1.25*figPos(3);
    fig.Position(4) = 1.5*figPos(4);
    fig.Position(2) = figPos(2) - 0.5*figPos(4);
    
    
    % Plot Vertical Profile
    subplot(4,1,[1,2]);
	plot(xODE(:,4),xODE(:,3),'lineWidth',1.5);
    ax = gca; xax = ax.XLim; yax = ax.YLim;
    hold on;
    if plotAllMethods
        plot(xRK(:,4),xRK(:,3),'g:',...
            xEM(:,4),xEM(:,3),'r:','lineWidth',1.5)  %r vs. h
    end
	xlabel('Range, m'), ylabel('Height, m'), grid on; axis equal;
    ax.XLim = xax; ax.YLim = yax;
    legend('ODE23','RK4','EULER');
    title('Paper Airplane Simulation Time History');
    
    % Plot time history of Flight Path Angle
    subplot(4,1,3);
    plot(tODE,wrapTo180(rad2deg(xODE(:,2))),'lineWidth',1.5);
    ax = gca; xax = ax.XLim; yax = ax.YLim;
    hold on;
    if plotAllMethods
        plot(tRK,wrapTo180(rad2deg(xRK(:,2))),'g:',...
            tEM,wrapTo180(rad2deg(xEM(:,2))),'r:','lineWidth',1.5);     %t vs. gamma
    end
    ax.XLim = xax; ax.YLim = yax;
    xlabel('Time (sec)'); ylabel('Gamma (Deg)'); grid on;
    
    % Plot time history of Airspeed
    subplot(4,1,4);
    plot(tODE,xODE(:,1),'lineWidth',1.5);
    hold on;
    ax = gca; xax = ax.XLim; yax = ax.YLim;
    if plotAllMethods
        plot(tRK,xRK(:,1),'g:',...
            tEM,xEM(:,1),'r:','lineWidth',1.5);     %t vs. V
    end
    ax.XLim = xax; ax.YLim = yax;
    xlabel('Time (sec)'); ylabel('Velocity (m/s)'); grid on;

    end
    
