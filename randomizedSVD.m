function [U, S, V] = randomizedSVD(A, k, p)
    % Fix the seed
    seedValue = 1234; 
    rng(seedValue, 'twister');
    % A: input matrix
    % k: target rank
    % p: oversampling parameter

    [m, n] = size(A);
    Omega = randn(n, k + p); % Step 1: Generate a random matrix
    Y = A * Omega;           % Step 2: Multiply A by random matrix
    [Q, ~] = qr(Y, 0);       % Compute the QR decomposition, we're interested in Q
    B = Q' * A;              % Form a smaller matrix
    [U_tilda, S, V] = svd(B, 'econ'); % Step 3: Compute the SVD of the smaller matrix
    U = Q * U_tilda;         % Step 4: Approximate U of the original matrix A

end
