% Erik P Bertelli 
% Dec 1 2016
% UC Berkeley

% Runs the COVSEL precsion matrix approximation on supreme court data

X_raw = csvread(strcat(pwd,'\2010_2015_SCOTUS_term.csv'));

[n,d] = size(X_raw);
exclude_records = [];

% exclude the incomplete and unanimous records
for i = 1:n
    if sum(X_raw(i,:) == 0) > 0
        % add that row to missing_records
        exclude_records = [exclude_records, i];
    elseif sum(X_raw(i,:)) == 9
        % exclude 
        exclude_records = [exclude_records, i];
    elseif sum(X_raw(i,:)) == 18
        % exclude 
        exclude_records = [exclude_records, i];
    end
end

X = X_raw;
X(exclude_records,:) = [];
X = X - 1;
X(X==0) = -1;

% empirical correlation matrix
S = cov(X);

% Code from Banjeree, El Ghaoui, and d'Aspremont
% this code finds the l-1 penalized precision matrix
rho=.2;             % Controls sparsity
prec=1e-3;          % Numerical precision
maxiter=10;        % Maximum Number of sweeps in coordinate descent
algot='nest';       % Algorith for solving BoxQP: 'sedumi' uses SEDUMI, 'nest' uses smooth minimization (usually faster)
maxnest=1000;       % Maximum number of iterations in the smooth BoxQP solver
% return estimate of the precision matrix
[A,~] = spmlcdvec(S,rho,maxiter,prec,maxnest,algot);
% end code from Banjeree, El Ghaoui, and d'Aspremont


% Create the outputs
A_thresh = A;
A_thresh(abs(A_thresh) < .01) = 0;

% inversion does away with sparsity, so I have to threshold
theta = inv(A_thresh);
theta_thresh = theta;
theta_thresh(abs(theta_thresh) < .03) = 0;

justices = {'AScalia','CThomas','SAAlito','JGRoberts','AMKennedy','SGBreyer','EKagan','RBGinsburg','SSotomayor'};

figure
G = graph(theta_thresh,justices,'OmitSelfLoops');
H = plot(G);

% define negative edges to be red
[s_neg,t_neg] = find(theta_thresh<0);
for k = 1:length(s_neg)
    highlight(H,[s_neg(k), t_neg(k)],'EdgeColor', 'r')
end

% make edge weights proportional to theta value
abs_theta = abs(theta_thresh);
smallest_edge = min(min(abs_theta(abs_theta > 0)));
for i = 1:9
    for j = 1:9
           highlight(H,[i,j],'LineWidth', .5 * (abs(theta(i,j))/smallest_edge))
    end
end

title('Ising Model fit via Gaussian Relaxation')


