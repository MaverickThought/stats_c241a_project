function [L, theta, mu_hat, k, theta_disp] = IPF(X,C,maxIter,threshold)
% This runs IPF on a matrix of observations X given clique set C as a cell
% array. I assume that the data X is binary and the cliques are pairwise

% Inputs
% X - data vector of n iid samples of d dim r.vs
% C - list of cliques input as a cell vector (assumed to be edges)
% maxIter - the max number of iterations
% threshold - the threshold at which to stop the Log Likelihood updates

% Outputs
% L - Vector of log likelihoods indexed by iteration k
% theta - model parameters on the cliques (edges in this case)
% mu_hat - empirical marginals on cliques
% k - total number of iterations 
% theta_disp - a matrix display of theta - cliques by values of(s,t)

% define the problem size
[d,n] = size(X);
[~,c] = size(C);

% intialize mu_hat, the empirical marginals
mu_hat = zeros(c,2,2);

% Calculate empirical marginals
for i = 1:c
    % set the start and end nodes of the edge in clique i
    s = C{i}(1);
    t = C{i}(2);
    
    for x_s = 1:2
        for x_t = 1:2
            % Create indicators if node s = value x_s-1
            I_s__x_s = (X(s,:) == (x_s-1));
            I_t__x_t = (X(t,:) == (x_t-1));
            % Assign to mu_hat
            mu_hat(i,x_s,x_t) = (1/n) * sum(I_s__x_s.*I_t__x_t);
        end
    end      
end

% Check the empirical marginals
for i = 1:c
    if sum(sum(mu_hat(i,:,:))) < .999
        error('empirical marginals do not add up to 1 in each clique')
    else
        % do nothing
    end
end

% Intialize model parameters
% From Jordan's reader, I know that the Z is the same over all iterations
% if we intialize the thetas to 0 (equivalent to phi as 1) then it will be
% easy to find the correct value of Z. It will just be the total number of
% possible permutations (since the phi's imply they are all equally
% possible)
theta = zeros(c,2,2); 
% Theta(i,x_s,x_t) is theta for clique i = (s,t) when s=x_s-1 and t=x_t-1
Z = 2^d;
mu = zeros(c,2,2); % intialize the model marginals

% Initialize the log Likelihood L
L = zeros(1,maxIter);
L(1) = n * sum(sum(sum(mu_hat .* theta))) - n*log(Z);

k = 2;

% While not at maxIter and the difference above the threshold, loop
while (k < maxIter) && (k == 2 || abs(L(k-1) - L(k-2)) > threshold)
    % over cliques, setting s and t
    for i = 1:c
        s = C{i}(1);
        t = C{i}(2);
        
        % over values of nodes s and t
        for x_s = 1:2
            for x_t = 1:2
                % Compute the marginals based on the model (mu)
                % I will do it directly as suggested on Piazza
                % For cleanliness, I will do it in a separate function
                mu(i,x_s,x_t) = marginalize(theta,d,s,t,x_s,x_t,Z,C,c);
                
                % Update the thetas
                theta(i,x_s,x_t) = theta(i,x_s,x_t) + log(mu_hat(i,x_s,x_t)) - log(mu(i,x_s,x_t)) ; 
            end
        end
    end
    
    % Now that all of the thetas are updated, Update the Log Likelihood
    L(k) = n * sum(sum(sum(mu_hat .* theta))) - n*log(Z);

    k = k+1;
    % iterate again until the while loop conditions are broken
       % either we hit the iteration limit
       % or we have reached the threshold
end

% we already iterated the counter, so we need to decrease by 1
k = k - 1;

% to aid in display, create a 2 dim representation of theta
theta_disp = zeros(c,4);
for i=1:c
    % i will order the values (0,0), (0,1), (1,0), (1,1)
    theta_disp(i,:) = [theta(i,1,1) theta(i,1,2) theta(i,2,1) theta(i,2,2)];
end

end

function mu = marginalize(theta,d,s,t,x_s,x_t,Z,C,c)
% This subfunction finds the marginal value of a clique directly

% Input arguments 
% theta - model parameters
% d - dim of variables
% s - start node of edge
% t - end node of edge
% x_s - value of node s
% x_t - value of node t
% Z = normilization constant
% C = clique list in cell array form (assumed to be pairwise)
% c = number of cliques

% Outputs
% mu - the value of the marginal

% Find all possible permutations of the nodes that are not s,t
% Can do this iteratively by looping over dimension d-2
% Intialize as empty
Permu = [1;
         0];
for i=2:(d-2)
    % In each step, add a colume of ones and zeros to the existing Permu
    Permu = [Permu ones(2^(i-1),1) ; 
             Permu zeros(2^(i-1),1)];
end
      
% Now generate vectors cooresponding to fixed values for nodes s and t
% For d-dim problem, there are d-2 free variables in each set
% So that means there are 2^(d-2) permutations
m = 2^(d-2); % the number of permuations

if x_s == 1
    S = zeros(m,1);
else
    S = ones(m,1);
end

if x_t == 1
    T = zeros(m,1);
else
    T = ones(m,1);
end

% Now we need to insert the S and T vectors into permutations
% First insert S after node s-1
Permu = [Permu(:,1:s-1), S, Permu(:,s:d-2)]; 

% now insert T after node t-1
Permu = [Permu(:,1:t-1), T, Permu(:,t:d-1)]; 

% In order to translate this to addresses, need to add 1 to Permu
Permu = Permu + 1;

mu = 0; % intialize marginal value mu as 0
% Now calculate the marginals by iterating over each permuations
for j = 1:m
    % I will do this in log probability space so everything is addition
    log_p = -log(Z); %first include the normalizing constant
    % Over each clique
    for i = 1:c
        % identify the start and end nodes of clique C{i}
        s = C{i}(1);
        t = C{i}(2);
        
        % Take existing log_p and add model value for this permuation
        log_p = log_p + theta(i,Permu(j,s),Permu(j,t));        
    end
    % Take from log probability space back into probability space
    mu = mu + exp(log_p);
end

end


