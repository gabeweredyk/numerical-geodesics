
%Define a parameterized torus

% 0 < r < R
% $$$ r = 2;
% $$$ R = 5;
% $$$ 
% $$$ s = @(u, v) [ r * sin (v) ; (R + r * cos(v)) * sin(u); (R + r * cos(v)) * cos(u) ];
% $$$ 
% $$$ % define the ample partial derivatives
% $$$ 
% $$$ su = @(u, v) [ 0; (R + r * cos(v)) * cos(u); -(R + r * cos(v)) * sin(u) ];
% $$$ sv = @(u, v) [ r * cos(v) ; - r * sin(v) * sin(u) ; - r * sin(v) * cos(u) ];
% $$$ 
% $$$ suu = @(u, v) [ 0; - (R + r * cos(v)) * sin(u); - (R + r * cos(v)) * cos(u) ];
% $$$ suv = @(u, v) [ 0; - r * sin(v) * cos(u) ; r * sin(v) * sin(u) ];
% $$$ svv = @(u, v) [ - r * sin(v) ;  - r* cos(v) * sin(u);   - r* cos(v) * sin(u)] ;


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

