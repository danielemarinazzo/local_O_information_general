function prob = ProbabilityLocalDiscrete(data)

[T,N] = size(data);

% Find max value, needed for the choice of the base
M = max(max(data)) + 1;

% Each configuration is a number in base M
dec = VecToDecimal(data,M);
edges = 0:M^N;
counts = histcounts(dec,edges);

% Probs of individual configurations
% prob = reshape(counts / T, M*ones(1,N));
% prob = permute(prob, N:-1:1);

% Prob for every event of the process 
prob = zeros(T,1);

for n = 0:(M^N-1)
    prob(dec == n) = counts(n+1)/T;
end

end

