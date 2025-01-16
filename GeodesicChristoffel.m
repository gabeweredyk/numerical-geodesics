% Define the differential equation for a Geodesic

function X = GeodesicChristoffel(u, v, ut, vt, s, su, sv, suu, suv, svv)

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

    % Define the "norm" of the surface
    N = 2 * (E .* G - F.^2);
    
    % Define the Christofel Symbols
    Gamma111 = ( G * Eu - 2 * F * Fu + F * Ev ) / N;
    Gamma112 = ( G * Ev - F * Gu ) / N;
    Gamma122 = ( 2 * G * Fv - G * Gu - F * Gv) / N;
    Gamma211 = ( 2 * E * Fu - E * Ev - F * Eu) / N;
    Gamma212 = (E * Gu - F * Ev) / N;
    Gamma222 = (E * Gv - 2 * F * Fv + F * Gu) / N;

    X = -[-ut; -vt; Gamma111 * ut^2 + 2 * Gamma112 * ut * vt + Gamma122 * vt^2; Gamma211 * ut^2 + 2 * Gamma212 * ut * vt + Gamma222 * vt^2];
    
end
