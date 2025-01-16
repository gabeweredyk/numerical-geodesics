% Runge-Kutta scheme for an autonoumus system of differential equations
function [t, y] = rk4aut(f, tIn, tFinal, yIn, nstep)

    % Determine the dimension of the input space of f
    m = max(size(yIn));
    
    % Allocate space for storing time & solution
    t = zeros(1,nstep+1);
    y = zeros(m,nstep+1);

    % Define the time step
    dt = (tFinal - tIn)/nstep;

    % Initiate the loop
    t(1) = tIn;
    y(:,1) = yIn;

    % Implement Runge-Kutta Update Scheme
    for j = 1:nstep
        t(j + 1) = t(j) + dt;
        k1 = f( y(:,j) );
        k2 = f( y(:,j) + 0.5*dt*k1);
        k3 = f(y(:,j) + 0.5*dt*k2);
        k4 = f(y(:,j) + dt*k3);
        y(:,j+1) = y(:,j) + dt*(k1/6 + k2/3 + k3/3 + k4/6);
    end

end
