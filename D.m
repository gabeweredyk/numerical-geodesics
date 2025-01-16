% Estimate Derivate with Divided Difference operator
% Inputs:
% f - The function to be partly differentiated
%   - NOTE: f must output a column vector
% n - The degree of the derivative (how many times)
% x - Variable to which f is being differentiated
% h - Small gap by default 10^-5
% m - The degree of accuracy will be 2*m + 1 (in a symmetric scheme) If not specified, m will be 3
% b - The center of the symmetric scheme, if not specified will be 0
% Outputs:
% ft - The derivative

function fxn = D(f, n, variable, h, m, b)
    
    % Assign default values
    if (nargin < 6) b = 0; end
    if (nargin < 5) m = 3; end
    if (nargin < 4) h = 1e-3; end

    % Define the # of unknowns 
    M = 2*m + 1;
    
    % Form vector of shifts
    shifts = linspace(-m + b,m + b,M);

    % allocate space for matrix
    A = ones(M, M);

    % Form matrix based on Taylor-Series coefficients
    for k = 2:M
        A(k,:) = shifts.^(k - 1);
    end
    
    % Set up the equation to solve for
    b = zeros(M,1);
    b(n + 1) = factorial(n);

    % Solve the equation w/ backslash
    c =  A \ b;
    c = (1 / h^n) * c;

    % Define a matrix to back multiply a shifted f with
    C = diag(c);

    % Define a vector to "normalize" the matlab operation
    N = ones(size(shifts));
    
    % Define the new function
    switch variable
      case 't'
        fxn = @(t) c' * f(t + h*shifts);
      case 'u'
        fxn = @(u, v) sum( ( f(u + h*shifts, v .* N) * C), 2);
      case 'v'
        fxn = @(u, v) sum((f(u.*N, v + h*shifts) * C), 2);
      case 'uv'
        [Su, Sv] = meshgrid(shifts);

        subindex = @(A, r, c) A(r, c);

        C = c * c';
        
        f1 = @(u, v, n) sum( subindex(f(u + h*Su, v + h*Sv) * C, (M*(n - 1) + 1):(n*M), 1:M ) , 'all') ;
        fxn = @(u, v) [ f1(u,v, 1) ; f1(u, v, 2); f1(u, v, 3) ]; 
    end 
end
