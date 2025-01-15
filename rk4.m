function [t, y] = rk4(f, a, b, ya, nstep)

m = max(size(ya));

t = zeros(1,nstep+1);
y = zeros(m,nstep+1);
dt = (b - a)/nstep;
t(1) = a;
y(:,1) = ya;

for j = 1:nstep
    t(j + 1) = t(j) + dt;
    k1 = f(t(j), y(:,j));
    k2 = f(t(j) + dt/2, y(:,j) + 0.5*dt*k1);
    k3 = f(t(j) + dt/2, y(:,j) + 0.5*dt*k2);
    k4 = f(t(j) + dt, y(:,j) + dt*k3);
    y(:,j+1) = y(:,j) + dt*(k1/6 + k2/3 + k3/3 + k4/6);
end

end