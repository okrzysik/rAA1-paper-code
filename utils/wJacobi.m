function [vnew, Romega] = wJacobi(A, vold, f, options)
% Performs one weighted Jacobi iteration on Au=f with weight omega.
% vold is the approximation before the Jacobi iteration, and vnew after. 
% Make sure to implement the method in matrix form (this is much
% faster than implementing it per component), and to use 
% sparse matrices for all steps.

omega = options.omega;

n = size(A, 1); % Size of square matrix A.
diagA = diag(A); % Diagonal components of A.
D = sparse(1:n, 1:n, diagA, n, n); % Diagonal component of A in a matrix.
invD = sparse(1:n, 1:n, 1./diagA, n, n); % Inverse of diagonal component matrix of A.

% Note that: A = D - L - U => L + U = D - A
RJ = invD*(D - A); % Jacobi iteration matrix.
Romega = (1 - omega)*speye(n) + omega*RJ; % Weighted Jacobi iteration matrix.

vnew = Romega*vold + omega*invD*f; % One sweep of weighted Jacobi on v.

end