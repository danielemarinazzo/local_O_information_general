function O_inf = OI_Local(x,method)

% Calculates the dynamic local O-information for multivariate stochastic processes.
% Works in the case of discrete systems (with positive integer values), of
% non-Gaussian continuous systems (with kernel methods for estimating the
% probability) or for Gaussian systems.
% 
% INPUT 
% x         -   Time series [T,N] T samples, N variables
% method    -   'discrete','continous','gaussian'

[T,N] = size(x);

% 
switch lower(method)
    case 'discrete'
        probFunction = @(x) ProbabilityLocalDiscrete(x);
    case 'continous'
        probFunction = @(x) ProbabilityLocalKernel(x);        
    case 'gaussian'
        probFunction = @(x) ProbabilityLocalGaussian(x);        
    otherwise
        O_inf = NaN;
        disp('Tipo di dati non corretto');
        return
end

% Initialize output
O_inf = zeros(T,1);

% Loop on series, compute entropies
for i = 1:N
    % Local prob for X_i
    p_i = probFunction(x(:,i));
    
    % Prob for all variables but X_i
    p_not_i = probFunction(x(:,setdiff(1:N,i)));
    
    % Compute the two entropies
    h_i = -log2(p_i);
    h_not_i = -log2(p_not_i);
    
    % Compoute difference and sum
    O_inf = O_inf + (h_i - h_not_i);
end

% Last piece with total entropy
p = probFunction(x);
h = -log2(p);

O_inf = O_inf + (N-2)*h;
end

