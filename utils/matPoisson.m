function A = matPoisson(N)

% Builds the problem matrix A for (N-1)x(N-1) unknowns 
% (interior points) on a grid with grid spacing h=1/N.

% This seems to be slower than the code below for moderate N, e.g. N ~ 300.
% I = speye(N-1); % Identity matrix.
% E = sparse(2:N-1, 1:N-2, 1 , N-1, N-1); % Ones along sub diagonal.
% D = 2*I - E - E.'; % 2's on diag, -1's on sub and super diagonals. 
% A = N^2*(kron(D,I) + kron(I,D));

A1dDiag = sparse(1:N-1, 1:N-1, 4*ones(1,N-1), N-1, N-1);
A1dSubDiag = -sparse(2:N-1, 1:N-2, ones(1,N-2) ,N-1 ,N-1);
A1d = A1dSubDiag + A1dDiag + A1dSubDiag'; % A for the 1-dimensional problem.
blockSubDiag = sparse( N:(N-1)*(N-1), 1:(N-1)*(N-2), 1, (N-1)^2, (N-1)^2); % The 1's down the middle of the block subdiagonal
A = N^2*(kron(speye(N-1), A1d) - blockSubDiag - blockSubDiag'); % Remember to multiply by N^2 = 1/h^2.

end