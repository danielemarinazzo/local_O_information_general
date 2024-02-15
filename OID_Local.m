function O_inf_dyn = OID_Local(x,p,target,method)

% Calculates the dynamic local O-information for multivariate stochastic processes.
% Works in the case of discrete systems (with positive integer values), of
% non-Gaussian continuous systems (with kernel methods for estimating the
% probability) or for Gaussian systems.
% 
% INPUT 
% x         -   Time series size [T,N] with T samples and N variables 
% p         -   Length of temporal embedding
% target    -   One or more target variables
% method    -   data type: 'discrete','continuous','gaussian'


% Samples and variables
[T,N] = size(x);

% Past vectors
x_past = [];

for i = 1:p
    x_past = cat(3,x_past,circshift(x,i,1));
end

% Cut to avoid border effect
x = x((p+1):end,:);
x_past = x_past((p+1):end,:,:);
T = size(x,1);

% Choose method to compute probs
switch lower(method)
    case 'discrete'
        probFunction = @(x) ProbabilityLocalDiscrete(x);
    case 'continuous'
        probFunction = @(x) ProbabilityLocalKernel(x);        
    case 'gaussian'
        probFunction = @(x) ProbabilityLocalGaussian(x);        
    otherwise
        O_inf_dyn = NaN;
        disp('unknown data type');
        return
end

% Initialize output
O_inf_dyn = zeros(T,1);

% Target time series
y = x(:,target); % Target present
Y = squeeze(x_past(:,target,:)); % Target past

% Reshape for more than one target
Y = reshape(Y,[T,p*size(target,2)]);

% Compute probs for MI
p_y_Y = probFunction([y,Y]);
p_Y = probFunction(Y);

% Consider only the time series different from target
X = x_past(:,setdiff(1:N,target),:);
n = size(X,2);

% Loop over all time series, compute MI
for i = 1:n
    % X - X_i term
    Xi = squeeze(X(:,setdiff(1:n,i),:)); 

    % All lags in one dimension
    Xi = reshape(Xi,[T,p*(n-1)]);

    % Compute the two other probs needed for
    % i(y ; X_i | Y) = log2(  (p(y,X_i,Y) * p(Y))  /  (p(y,Y) * p(X_i,Y))  )
    p_y_Xi_Y = probFunction([y,Xi,Y]);
    p_Xi_Y = probFunction([Xi,Y]);

    % Local MI
    mi_not_i = log2( (p_y_Xi_Y .* p_Y) ./ (p_Xi_Y .* p_y_Y) );

    % Add up
    O_inf_dyn = O_inf_dyn + mi_not_i;
end

% Last part with total MI
X = reshape(X,[T,p*n]);

p_y_X_Y = probFunction([y,X,Y]);
p_X_Y = probFunction([X,Y]);

mi = log2( (p_y_X_Y .* p_Y) ./ (p_X_Y .* p_y_Y) );

% Last contribution to dynamic O info
O_inf_dyn = O_inf_dyn + (1-n) * mi; 

end

