% Mike Hennessy and Jack Michaelis
% Hydraulic System Design Project
clear
clc

% Knowns

% Fluid constants
g = 9.8; % m/s^2
density = 1000; % [kg/m^3]
gamma = 980; % rho * g
viscosity = .00084; % [Pa*s]
epsilon = .00000015; % [m]

% NPS converted to m
d0 = 4 * 0.0254;
d1 = 2 * 0.0254;
d2 = d1;
d3 = d1;
d4 = d1;

% Pipe Areas [m^2]
a0 = (pi*(d0/2)^2);
a1 = (pi*(d1/2)^2);
a2 = a1;
a3 = a1;
a4 = a1;

% Pipe Lengths [m]
l0 = 34;
l1 = 10;
l2 = 23;
l3 = 18;
l4 = 0;

% Elevations [m]
z0 = 0;
z1 = 4;
z2 = 8;
z3 = 12;

% Pressure Values [Pa]
p0 = 15000;
p1 = 0;
p2 = 0;
p3 = 0;

% K Values
k_elbow = 1.5; % threaded (overestimate)
k_ball_valve = .05;
k_tee = 2; % branched, threaded (overestimate)
k_reducer = .43; % where A_2 / A_1 = 0.25

% Quantities of Valves, Unions, Tees, and Reducer
ball_valves_0 = 3; % quantity of ball valves for pipe 0
ball_valves_1 = 1;
ball_valves_2 = ball_valves_1;
ball_valves_3 = ball_valves_1;
ball_valves_4 = 0;
elbow_valves_0 = 2;
elbow_valves_1 = 0;
elbow_valves_2 = 1;
elbow_valves_3 = elbow_valves_2;
elbow_valves_4 = 0;
tee_valves_0 = 0;
tee_valves_1 = 1;
tee_valves_2 = tee_valves_1;
tee_valves_3 = tee_valves_1;
tee_valves_4 = tee_valves_1;
reducer = 1;

% Flow Rate and Velocity Calculations
supply = .02 + ((3 + 4) / 2) * 10^(-4); % supplied flow rate with March and April birthday months
v0 = supply / a0; % initial velocity

% Initial friction factor and Reynolds values
Re_0 = reynolds(density, v0, d0, viscosity);
f_0 = frictionFactor(Re_0, epsilon, d0);

% Unknowns
syms v1 v2 v3 v4 ph

% Matrix of all the diameters
d_all = [d0, d1, d2, d3, d4];

% Preallocating Reynolds Number and friction numbers
Re = [Re_0, Re_0, Re_0, Re_0, Re_0];
f = [f_0, f_0, f_0, f_0, f_0];

% Equation 1
Eq1 = p0/gamma + (v0^2)/(2*g) + z0 + ph == p1/gamma + (v1^2)/(2*g) + z1 + ((v0^2)/(2*g))*((f_0*l0)/d0 + ball_valves_0 * k_ball_valve + elbow_valves_0 * k_elbow + reducer * k_reducer) + ((v4^2)/(2*g))*((f(5)*l4)/d4 + tee_valves_4*k_tee) + ((v1^2)/(2*g))*((f(2)*l1)/d1 + tee_valves_1*k_tee + ball_valves_1*k_ball_valve);

% Equation 2
Eq2 = p0/gamma + (v0^2)/(2*g) + z0 + ph == p2/gamma + (v2^2)/(2*g) + z2 + ((v0^2)/(2*g))*((f_0*l0)/d0 + ball_valves_0 * k_ball_valve + elbow_valves_0 * k_elbow + reducer * k_reducer) + ((v4^2)/(2*g))*((f(5)*l4)/d4 + tee_valves_4*k_tee) + ((v2^2)/(2*g))*((f(3)*l2)/d2 + tee_valves_2*k_tee + ball_valves_2*k_ball_valve + elbow_valves_2*k_elbow);

% Equation 3
Eq3 = p0/gamma + (v0^2)/(2*g) + z0 + ph == p3/gamma + (v3^2)/(2*g) + z3 + ((v0^2)/(2*g))*((f_0*l0)/d0 + ball_valves_0 * k_ball_valve + elbow_valves_0 * k_elbow + reducer * k_reducer) + ((v3^2)/(2*g))*((f(4)*l3)/d3 + tee_valves_3*k_tee + ball_valves_3*k_ball_valve + elbow_valves_3*k_elbow);

% Equation 4
Eq4 = v0*a0 == v3*a3 + v4*a4;

% Equation 5
Eq5 = v4*a4 == v1*a1 + v2*a2;

% Solve initially for v1, v2, v3, v4, and ph
[v1, v2, v3, v4, ph] = solve([Eq1, Eq2, Eq3, Eq4, Eq5], [v1, v2, v3, v4, ph]);

response = filter(v1, v2, v3, v4, ph);

v1 = response(1);
v2 = response(2);
v3 = response(3);
v4 = response(4);
ph = response(5);

% Assign found variables
vguess = [v0, v1, v2, v3, v4];

% while loop variables
percent_change = 1;
old_values = zeros();
new_values = zeros();

% Loops through preselected iterations that runs the equations any number of times
while percent_change > .05
    % Store old values to calculate percent change at end of loop
    old_values = [v1, v2, v3, v4, ph];

    % Calculate Reynolds based on velocities
    for n = 2:5
        Re(n,:) = reynolds(density, vguess(:,n), d_all(:,n), viscosity);
        f(n, :) = frictionFactor(Re(n,:), epsilon, d_all(:,n));
    end

    syms v1 v2 v3 v4 ph

    % Equation 1
    Eq1 = p0/gamma + (v0^2)/(2*g) + z0 + ph == p1/gamma + (v1^2)/(2*g) + z1 + ((v0^2)/(2*g))*((f_0*l0)/d0 + ball_valves_0 * k_ball_valve + elbow_valves_0 * k_elbow + reducer * k_reducer) + ((v4^2)/(2*g))*((f(5)*l4)/d4 + tee_valves_4*k_tee) + ((v1^2)/(2*g))*((f(2)*l1)/d1 + tee_valves_1*k_tee + ball_valves_1*k_ball_valve);
    
    % Equation 2
    Eq2 = p0/gamma + (v0^2)/(2*g) + z0 + ph == p2/gamma + (v2^2)/(2*g) + z2 + ((v0^2)/(2*g))*((f_0*l0)/d0 + ball_valves_0 * k_ball_valve + elbow_valves_0 * k_elbow + reducer * k_reducer) + ((v4^2)/(2*g))*((f(5)*l4)/d4 + tee_valves_4*k_tee) + ((v2^2)/(2*g))*((f(3)*l2)/d2 + tee_valves_2*k_tee + ball_valves_2*k_ball_valve + elbow_valves_2*k_elbow);
    
    % Equation 3
    Eq3 = p0/gamma + (v0^2)/(2*g) + z0 + ph == p3/gamma + (v3^2)/(2*g) + z3 + ((v0^2)/(2*g))*((f_0*l0)/d0 + ball_valves_0 * k_ball_valve + elbow_valves_0 * k_elbow + reducer * k_reducer) + ((v3^2)/(2*g))*((f(4)*l3)/d3 + tee_valves_3*k_tee + ball_valves_3*k_ball_valve + elbow_valves_3*k_elbow);
    
    % Equation 4
    Eq4 = v0*a0 == v3*a3 + v4*a4;
    
    % Equation 5
    Eq5 = v4*a4 == v1*a1 + v2*a2;
    
    [v1, v2, v3, v4, ph] = solve([Eq1, Eq2, Eq3, Eq4, Eq5], [v1, v2, v3, v4, ph]);

    response = filter(v1, v2, v3, v4, ph);
    
    v1 = response(1);
    v2 = response(2);
    v3 = response(3);
    v4 = response(4);
    ph = response(5);

    % Reassigns vguess to contain new velocities
    vguess= [v0, v1, v2, v3, v4];
    new_values = [v1, v2, v3, v4, ph];

    percent_change = .05;
    for j = 1:5
        percent_change = max(percent_change, ((new_values(j) - old_values(j)) / old_values(j)));
    end

end

disp(new_values)

function [Re] = reynolds(density, v, d, viscosity)
    Re = (density*v*d)/viscosity;
end

function [f] = frictionFactor(Re, epsilon, d)
    if(Re > 2100)
        syms f;
        f = vpa(solve(-2*log10((epsilon/(d*3.7)) + (2.51./(Re*sqrt(f)))) == 1./sqrt(f), f));
    else
        f = 64 / Re;
    end
end

function [values] = filter(v1, v2, v3, v4, ph)
    for n = 1:4
        if(v1(n) >= 0 && v2(n) >= 0 && v3(n) >= 0 && v4(n) >= 0 && ph(n) >= 0)
            values = [v1(n), v2(n), v3(n), v4(n), ph(n)];
        end
    end
end