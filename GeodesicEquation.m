% Define the differential equation for a Geodesic

function X = GeodesicEquation(u, v, ut, vt, s, su, sv, suu, suv, svv)

    % Define the First Fundamental Form and their derivatives
    E = su(u, v)' * su(u, v);
    Eu = 2 * suu(u, v)' * su(u, v);
    Ev = 2 * suv(u, v)' * su(u, v);

    F = su(u, v)' * sv(u, v);
    Fu = suu(u, v)' * sv(u, v) + su(u, v)' * suv(u, v);
    Fv = suv(u, v)' * sv(u, v) + su(u, v)' * svv(u, v);

    G = sv(u, v)' * sv(u, v);
    Gu = 2 * suv(u, v)' * sv(u, v);
    Gv = 2 * svv(u, v)' * sv(u, v);

    % Define the "Magnitude" of the surface
    N = E * G - F.^2;

    % Define the inverse of First Fundamental Form matrix
    FFFinv = [G, -F; -F, E] / N;

    % Define commonly appearing matricies in the update scheme
    Q1 = -[Eu, Ev; Ev, 2 * Fv - Gu];
    Q2 = -[2 * Fu - Ev, Gu; Gu, Gv];

    Q = zeros(4,4);
    Q(1:2,1:2) = Q1;
    Q(3:4,3:4) = Q2;
    
    % Update the differential equation
    u = ut;
    v = vt;

    % Solve quadratic systems for ut and vt
    y = 0.5 * FFFinv * [ ut, vt, 0, 0; 0, 0, ut, vt] * Q * [ut; vt; ut; vt];

    % Pack solution into vector X
    X = [u; v; y];
end
