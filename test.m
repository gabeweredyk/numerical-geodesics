
%% Define the parameterization of the surface
%% E.g: Parameterization of a Torus
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
% $$$ 
% $$$ % Define range of evaluation
% $$$ au = -pi; bu = pi;
% $$$ av = -pi; bv = pi;

%% E.g: Parameterization of a sphere w/ radius R as a surface of revolution
% $$$ R = 2;
% $$$ 
% $$$ sx = @(u, v) R * sin(v);
% $$$ sy = @(u, v) R * cos(v) .* sin(u);
% $$$ sz = @(u, v) R * cos(v) .* cos(u);
% $$$ 
% $$$ % Define partial derivatives exactly
% $$$ sue = @(u, v) [ 0; R * cos(v) .* cos(u); -R * cos(v) .* sin(u) ];
% $$$ sve = @(u, v) [R * cos(v); -R * sin(v) .* sin(u) ; -R * sin(v) .* cos(u)];
% $$$ 
% $$$ suue = @(u, v) [0; -R * cos(v) .* sin(u); -R * cos(v) .* cos(u) ];
% $$$ suve = @(u, v) [0; -R * sin(v) .* cos(u); R * sin(v) .* sin(u) ];
% $$$ svve = @(u, v) -[R * sin(v); R * cos(v) .* sin(u); R * cos(v) .* cos(u)];
% $$$ 
% $$$ % Define region of evaluation of surface
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
% $$$ % With analytic derivatives
% $$$ sue = @(u, v) [0; 0; 1];
% $$$ sve = @(u, v) R * [ -sin(v); cos(v); 0 ];
% $$$ 
% $$$ suue = @(u, v) [0; 0; 0];
% $$$ suve = @(u, v) [0; 0; 0];
% $$$ svve = @(u, v) -R * [cos(v); sin(v); 0];
% $$$ 
% $$$ % Define range of parameterization
% $$$ h = 2;
% $$$ au = -h; bu = h;
% $$$ av = -pi; bv = pi;



% $$$ f = @(t) t.^5 + t.^3 - t + 3;
% $$$ 
% $$$ ftexact = @(t) 5*t.^4 + 3*t.^2 - 1;
% $$$ 
% $$$ ftcompute = D(f, 1, 't');

f = @(t) atan(t);

ftexact = @(t) 1 ./ ( 1 + t.^2);

ftcompute = D(f, 1, 't');

tvals = linspace(-2 * pi, 2 * pi, 101);

ue = ftexact(tvals);
un = ftcompute(tvals);
% $$$ 
% $$$ plot(tvals, un, 'b-o', tvals, ue, 'r-x');
plot(tvals, ue - un, 'k-x');
legend('Error');

%legend('Computed', 'Exact');

