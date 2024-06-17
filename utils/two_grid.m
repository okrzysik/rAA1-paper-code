% Update x that solves A*x = b using a two-grid method.
function x = two_grid(x, b, A, A_coarse, P, R, options)

    x        = wJacobi(A, x, b, options);
    r        = b - A*x;
    r_coarse = R*r;
    e_coarse = A_coarse \ r_coarse;
    x        = x + P*e_coarse;
    x        = wJacobi(A, x, b, options);
end