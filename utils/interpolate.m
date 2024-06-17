function I2hh = interpolate(N)

% Builds the interpolation matrix from level 2h to level h.
% There are (N2h-1)x(N2h-1) interior points on the coarse
% level (level 2h), and (N-1)x(N-1) interior points on the fine
% level (level h). This means that I2hh is a
% ((N-1)x(N-1)) x ((N2h-1)x(N2h-1)) matrix
% (many more rows than columns). Use the formulas
% on slides 56 and 57. I2hh is again a sparse matrix.

% Note that this script assumes N2h = N/2, i.e. that the coarse spacing is
% twice the fine grid spacing.


e = 0.5*ones(N-1,1); % Dummy matrix elements.
I2hh = spdiags([e 2*e e], -1:1, N-1, N-1); % A square matrix to extract 
% every second column from.
I2hh = I2hh(:, 2:2:N-1); % The 1-dimensional interpolation operator.
I2hh = kron(I2hh, I2hh); % The 2-dimensional operator is the kronecker 
% tensor product with itself (see Trottenberg et al. Remark 2.3.2)

end