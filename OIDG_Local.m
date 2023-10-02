function O_inf_dyn = OIDG_Local(x,p,method)

% Calculate the GLOBAL dynamic local O-information for
% Multivariate stochastic processes. By global I mean there is no target, but
% calculation a sort of dynamic integrated information.
% Works in the case of discrete systems (with positive integer values), of
% non-Gaussian continuous systems (with kernel methods for estimating the
% probability) or for Gaussian systems.
% 
% INPUT 
% x         -   Time series size [T,N] with T samples and N variables 
% p         -   Length of temporal embedding
% method    -   data type: 'discrete','continous','gaussian'

% Samples and variables
[T,N] = size(x);

% create vectors with the past 
x_past = [];

for i = 1:p
    x_past = cat(3,x_past,circshift(x,i,1));
end

% Cut to avoid border effects
x = x((p+1):end,:);
x_past = x_past((p+1):end,:,:);
T = size(x,1);

% Compute probabilities
switch lower(method)
    case 'discrete'
        probFunction = @(x) ProbabilityLocalDiscrete(x);
    case 'continous'
        probFunction = @(x) ProbabilityLocalKernel(x);        
    case 'gaussian'
        probFunction = @(x) ProbabilityLocalGaussian(x);        
    otherwise
        disp('Tipo di dati non corretto');
        return
end

% Initialize output
O_inf_dyn = zeros(T,1);

% Compute MI probabilities
p_x = probFunction(x);

% Loop across series and MI
for i = 1:N
    % Subtract
    Xi = squeeze(x_past(:,setdiff(1:N,i),:)); 

    % All lag in the same dimension
    Xi = reshape(Xi,[T,p*(N-1)]);

    % MI between the past of all variables but i, and the future of all
    p_x_Xi = probFunction([Xi,x]);
    p_Xi = probFunction(Xi);

    mi_not_i = log2( p_x_Xi ./ (p_Xi .* p_x) );

    % Adding up
    O_inf_dyn = O_inf_dyn + mi_not_i;
    
    
    
    
end

% Last step with total MI
X = reshape(x_past,[T,p*N]);

p_x_X = probFunction([X,x]);
p_X = probFunction(X);

mi = log2( p_x_X ./ (p_X .* p_x) );

% Last contribution to dynamic O info
% O_inf_dyn = O_inf_dyn + (1-N) * mi; 
O_inf_dyn = O_inf_dyn + (1-N) * mi; 

end

