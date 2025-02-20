# Ahshhsha% Jacobi Method to Solve Ax = b

% Input: Coefficient matrix A and RHS vector b A = input('Enter the coefficient matrix A: '); b = input('Enter the right-hand side vector b: ');

% Check for diagonal dominance if ~all(2 * abs(diag(A)) >= sum(abs(A), 2)) warning('Matrix A is not strictly diagonally dominant. The method may not converge.'); end

% Input initial guess, tolerance, and maximum iterations x0 = input('Enter the initial guess vector x0: '); tol = input('Enter the tolerance for convergence (e.g., 1e-6): '); max_iter = input('Enter the maximum number of iterations: ');

% Initialize variables n = length(b); x = x0; x_old = x0;

% Jacobi iteration loop for k = 1:max_iter for i = 1:n % Compute the summation for the i-th equation sum1 = A(i, :) * x_old - A(i, i) * x_old(i); x(i) = (b(i) - sum1) / A(i, i); end

% Check for convergence
if norm(x - x_old, Inf) < tol
    fprintf('Converged in %d iterations.\n', k);
    disp('Solution:');
    disp(x);
    return;
end

% Update the old solution
x_old = x;

end

% If the loop completes without convergence fprintf('Did not converge within %d iterations.\n', max_iter); disp('Last computed solution:'); disp(x);

