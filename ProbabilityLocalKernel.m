function prob = ProbabilityLocalKernel(x)

[T,N] = size(x);

% Rule of thumb for kernel width
s = std(x);
b = (4/(T*(2+N)))^(1/(N+4)) * s;

% Gaussian kernels (slow)
prob = mvksdensity(x,x,'Bandwidth',b);

end

