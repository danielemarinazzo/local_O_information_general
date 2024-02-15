function prob = ProbabilityLocalGaussian(x)

[T,N] = size(x);

% find covariance and average
Sigma = cov(x);
mu = mean(x,1);

% vectorize
mu_vec = repmat(mu,[T,1]);

% normalization constant
c = sqrt((2*pi)^N * det(Sigma));

% Gaussian prob for all samples
prob = exp(-0.5 * sum(  ( (x - mu_vec) / Sigma ) .* (x -mu_vec) ,2) ) / c;

end

