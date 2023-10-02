function [d_tc, d_dtc, d_o] = deltaOI_Local(x,method)

% Calculates the local O-information for multivariate stochastic processes.
% Works in the case of discrete systems (with positive integer values), of
% non-Gaussian continuous systems (with kernel methods for estimating the
% probability) or for Gaussian systems.
% 
% INPUT 
% x         -   Time series size [T,N] with T samples and N variables
% OUTPUT    
%           -   d_tc, size [T,N], d_tc(t,i) local variation of total 
%                     correlation at time t of variable i
%           -   d_dtc, size [T,N], d_tc(t,i) è local variation of dual 
%                     total correlation at time t of variable i        
% method    -   'discrete','continous','gaussian'

[T,N] = size(x);

switch lower(method)
    case 'discrete'
        probFunction = @(x) ProbabilityLocalDiscrete(x);
    case 'continous'
        probFunction = @(x) ProbabilityLocalKernel(x);        
    case 'gaussian'
        probFunction = @(x) ProbabilityLocalGaussian(x);        
    otherwise
        oi = NaN;
        disp('Tipo di dati non corretto');
        return
end



% Last bit with total entropy
p = probFunction(x);
h = -log2(p);

% Init output
d_o   = zeros(T,N);
d_tc  =  zeros(T,N);
d_dtc =  zeros(T,N);




% Loop on series, compute entropies
for i = 1:N
    var_not_i =  setdiff(1:N,i);
    % Local prob for X_i
    p_i = probFunction(x(:,i));
    
    % Prob for all variables but X_i
    p_not_i = probFunction(x(:,setdiff(1:N,i)));
    
    % Two entropies
    h_i = -log2(p_i);
    h_not_i = -log2(p_not_i);
    
    % Compute differences and sum
    d_o(:,i)  =   d_o(:,i) + (2-N)*( h_i + h_not_i - h);
    d_tc(:,i)  = h_i + h_not_i - h;
    d_dtc(:,i) = (N-1)*d_tc(:,i);
    for jj = 1:N-1
        j  = var_not_i(jj);
              
        p_not_ij   =  probFunction(x(:,setdiff(1:N,[i,j])));
        p_not_j    =  probFunction(x(:,setdiff(1:N,j))); 

        h_not_ij   =  -log2(p_not_ij);
        h_not_j    =  -log2(p_not_j);
        
        d_dtc(:,i) =   d_dtc(:,i) - (h_i + h_not_ij - h_not_j);
        d_o(:,i)   =   d_o(:,i) + (h_i + h_not_ij - h_not_j);
    end
end




end