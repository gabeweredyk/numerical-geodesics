% Reset MatLab
clf; clearvars; clc;

%% Define the parameterization of the surface
% E.g: Parameterization of a Torus
% $$$ % 0 < r < R
% $$$ r = 1;
% $$$ R = 2;
% $$$ 
% $$$ % Define the individual components of the parameterization
% $$$ sx = @(u, v) r * sin(v);
% $$$ sy = @(u, v) (R + r * cos(v) ) .* sin(u);
% $$$ sz = @(u, v) (R + r * cos(v) ) .* cos(u);
% $$$ 
% $$$ % Define the ample partial derivatives analytically
% $$$ sue = @(u, v) [ 0; (R + r * cos(v)) .* cos(u); -(R + r * cos(v)) .* sin(u) ];
% $$$ sve = @(u, v) [ r * cos(v) ; - r * sin(v) .* sin(u) ; - r * sin(v) .* cos(u) ];
% $$$ 
% $$$ suue = @(u, v) [ 0; - (R + r * cos(v)) .* sin(u); - (R + r * cos(v)) .* cos(u) ];
% $$$ suve = @(u, v) [ 0; - r * sin(v) .* cos(u) ; r * sin(v) .* sin(u) ];
% $$$ svve = @(u, v) [ - r * sin(v) ;  - r* cos(v) .* sin(u);   - r* cos(v) .* sin(u)] ;

% E.g: Parameterization of a sphere w/ radius R
R = 2;

sx = @(u, v) R * sin(v);
sy = @(u, v) R * cos(v) .* sin(u);
sz = @(u, v) R * cos(v) .* cos(u);

% Define partial derivatives exactly
sue = @(u, v) [ 0; R * cos(v) .* cos(u); -R * cos(v) .* sin(u) ];
sve = @(u, v) [R * cos(v); -R * sin(v) .* sin(u) ; -R * sin(v) .* cos(v)];

suue = @(u, v) [0; -R * cos(v) .* sin(u); -R * cos(v) .* cos(u) ];
suve = @(u, v) [0; -R * sin(v) .* cos(u); -R * sin(v)  .* sin(u) ];
svve = @(u, v) [-R * sin(v); R * cos(v) .* sin(u); R * cos(v) .* cos(u)];

% $$$ % E.g.: Parameterization of a cyllinder
% $$$ R = 2;
% $$$ sx = @(u, v) R * cos(v);
% $$$ sy = @(u, v) R * sin(v);
% $$$ sz = @(u, v) u;

% Use the individually defined components
s = @(u, v) [sx(u, v); sy(u, v); sz(u, v)];

% Define derivatives numerically if not previously given
suc = D(s, 1, 'u');
svc = D(s, 1, 'v');

suuc = D(s, 2, 'u');
suvc = D(s, 1, 'uv');
svvc = D(s, 2, 'v');

% Define region of evaluation of surface
au = -pi; bu = pi;
av = -pi; bv = pi;

% Define number of grid cells used
Nu = 50; Nv = 50;

% Generte grids in R2 for which s should be evaluated
[uvals, vvals] = meshgrid( linspace(au,bu,Nu + 1), linspace(av, bv, Nv + 1) );

% Store the value of the evaluations
X = sx(uvals, vvals);
Y = sy(uvals, vvals);
Z = sz(uvals, vvals);

% Define the differential equation to solve for as per GeodesicEquations.m
Fe = @(X) GeodesicEquation( X(1), X(2), X(3), X(4), s, sue, sve, suue, suve, svve);
Fc = @(X) GeodesicEquation( X(1), X(2), X(3), X(4), s, suc, svc, suuc, suvc, svvc);

% Define initial conditions of the problem
uin = 1; vin = 2;
% Choose a starting "vector"
theta = pi / 4;
utin = cos(theta); vtin = sin(theta);

Fin = [uin; vin; utin; vtin];

tIn = 0;
tFinal = 3;

nstep = ceil( 1000 * tFinal );

[tvals, Fevals] = rk4aut(Fe, tIn, tFinal, Fin, nstep);
%[tvals, Fcvals] = rk4aut(Fc, tIn, tFinal, Fin, nstep);

gammae = s( Fevals(1,:), Fevals(2,:) );
%gammac = s( Fcvals(1,:), Fcvals(2,:) );

% Plot the results
hold on
% Plot the geodesic based on exact partial derivatives
plot3(gammae(1,:), gammae(2,:), gammae(3,:), 'g-', 'LineWidth', 8);
% Plot the geodesic based on numerical partial derivatives
%plot3(gammac(1,:), gammac(2,:), gammac(3,:), 'r-', 'LineWidth', 8);
% Plot the surface
surf(X, Y, Z);
% Specified 
axis([-3 3 -3 3 -3 3]);
hold off
