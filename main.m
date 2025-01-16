% Reset MatLab
clf; clearvars; clc;

%% Define the parameterization of the surface
%% E.g: Parameterization of a Torus
% 0 < r < R
r = 1;
R = 2;

% Define the individual components of the parameterization
sx = @(u, v) r * sin(v);
sy = @(u, v) (R + r * cos(v) ) .* sin(u);
sz = @(u, v) (R + r * cos(v) ) .* cos(u);

% Define range of evaluation
au = -pi; bu = pi;
av = -pi; bv = pi;

%% E.g: Parameterization of a sphere w/ radius R as a surface of revolution
% $$$ R = 2;
% $$$ 
% $$$ sx = @(u, v) R * sin(v);
% $$$ sy = @(u, v) R * cos(v) .* sin(u);
% $$$ sz = @(u, v) R * cos(v) .* cos(u);
% $$$ 
% $$$ % Define region of evaluation 
% $$$ au = -pi; bu = pi;
% $$$ av = -pi / 2; bv = pi /2;

%% Stereographic Parameterization of the sphere (of radius R)
% $$$ R = 2;
% $$$ sx = @(u, v) 2 * u ./ (u.^2 + v.^2 + 1);
% $$$ sy = @(u, v) 2 * v ./ (u.^2 + v.^2 + 1);
% $$$ sz = @(u, v) 1 - ( 2 ./ (u.^2 + v.^2 + 1) );

%% E.g.: Parameterization of a cyllinder
% $$$ R = 2;
% $$$ sx = @(u, v) R * cos(v);
% $$$ sy = @(u, v) R * sin(v);
% $$$ sz = @(u, v) u;
% $$$ 
% $$$ % Define range of parameterization
% $$$ h = 2;
% $$$ au = -h; bu = h;
% $$$ av = -pi; bv = pi;

%% Parameterization of an ellipsoid
% $$$ Rx = 2; Ry = 2.5; Rz = 1.5;
% $$$ 
% $$$ sx = @(u, v) Rx * sin(v);
% $$$ sy = @(u, v) Ry * cos(v) .* sin(u);
% $$$ sz = @(u, v) Rz * cos(v) .* cos(u);
% $$$ 
% $$$ % Define region of evaluation 
% $$$ au = -pi; bu = pi;
% $$$ av = -pi / 2; bv = pi /2;


%% Define the parameteric surface & its PDs
% Use the individually defined components
s = @(u, v) [sx(u, v); sy(u, v); sz(u, v)];

% Define derivatives numerically w/ D.m
su = D(s, 1, 'u');
sv = D(s, 1, 'v');

suu = D(s, 2, 'u');
svv = D(s, 2, 'v');
% The uv PD has to be done component wise
sxuv = D(sx, 1, 'uv');
syuv = D(sy, 1, 'uv');
szuv = D(sz, 1, 'uv');
suv = @(u, v) [sxuv(u, v); syuv(u, v); szuv(u, v)];

%% Numerically solve for the path of the geodesic
% Define the differential equation to solve for as per GeodesicEquations.m
F = @(X) GeodesicEquation( X(1), X(2), X(3), X(4), s, su, sv, suu, suv, svv);

% Define initial conditions of the problem
uin = -pi/3; vin = -pi/3;
% Choose a starting "vector"
theta = -pi / 4;
utin = cos(theta); vtin = sin(theta);

% Construct initial vector
Fin = [uin; vin; utin; vtin];

% Specify time range of evaluation
tIn = 0;
tFinal = 50;

% Decide the # of steps the program should take
nstep = ceil( 400 * tFinal );

% Use Runge-Kutta 4 to update scheme as defined in rk4aut.m
% Note the Geodesic Equations are independent of time
[tvals, Fvals] = rk4aut(F, tIn, tFinal, Fin, nstep);

% Store the evaluation of the path in the uv-plane on s
gamma = s( Fvals(1,:), Fvals(2,:) );

%% Plot the results

% First we need to sample points on s
% Define number of grid cells used
Nu = 100; Nv = 100;

% Generte grids in R2 for which s should be evaluated
[uvals, vvals] = meshgrid( linspace(au,bu,Nu + 1), linspace(av, bv, Nv + 1) );

% Store the value of the evaluations
X = sx(uvals, vvals);
Y = sy(uvals, vvals);
Z = sz(uvals, vvals);

hold on
% Plot the geodesic based on numerical partial derivatives
plot3(gamma(1,:), gamma(2,:), gamma(3,:), 'r-', 'LineWidth', 8);
% Plot the surface
surf(X, Y, Z);
% Specified 
axis([-3 3 -3 3 -3 3]);
% Legend
legend('Exact \sigma_{uv}');
hold off
