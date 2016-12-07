% Erik P Bertelli 
% Nov 22 2016
% UC Berkeley

% This code runs the IPF algoritm on historical SCOTUS votes to fit a 
% simple Pairwise Markov Random Field model
X_raw = csvread(strcat(pwd,'\2010_2015_SCOTUS_term.csv'));
maxIter = 100;
threshold = 0.0001;

[n,d] = size(X_raw);
missing_records = [];

% exclue incomplete votes 
for i = 1:n
    if sum(X_raw(i,:) == 0) > 0
        % add that row to missing_records
        missing_records = [missing_records, i];
    end
end
X = X_raw;
X(missing_records,:) = [];

% translate to (0,1) instead of (1,2) binary variables 
X = X - 1;
 
% Create potential edge sets
E = cell(1,36);
k = 1;
for i = 1:8
    for j = i+1:9
        E{k} = [i,j];
        k = k+1;
    end
end

X = X'; % transpose the matrix to match IPF algorithm

% Run IPF
[L, theta, mu_hat, k, theta_disp] = IPF(X,E,maxIter,threshold);

% define justices labels
justices = {'AScalia','CThomas','SAAlito','JGRoberts','AMKennedy','SGBreyer','EKagan','RBGinsburg','SSotomayor'};

% Construct 9-0 conservative case
theta_9_0 = theta(:,1,1);
A = zeros(9);
for i = 1:9
    for j = 1:9
        if i == j
            A(i,j) = 0;
        elseif i > j
            A(i,j) = 0;
        else
            A(i,j) = theta_9_0(1);
            theta_9_0(1) = [];
        end
    end
end

% convert triangular A to symmetric and threshold
A = A + A';
A_thresh = A;
A_thresh(A_thresh<.1) = 0;
G = graph(A_thresh , justices);

% plot
H = plot(G)%,'EdgeLabel',G.Edges.Weight)

% make line width proportional to parameter value
[s,t] = find(A_thresh>0)
for k = 1:length(s)
    highlight(H,[s(k), t(k)],'LineWidth', .5 * (A_thresh(s(k),t(k))/0.1))
end
title('MRF via IPF - 9-0 Conservative Decision')

% Create a 5-4 Conservative case
theta_5_4 = [];
for i = 1:size(theta,1)
    edge = E{i};
    if edge(1) < 6
        % conservative justice
        if edge(2) < 6
            % conservative justice
            theta_5_4 = [theta_5_4, theta(i,1,1)];
        else
            % liberal justice
            theta_5_4 = [theta_5_4, theta(i,1,2)];
        end
    else
        % liberal justice
        if edge(2) < 6
            % conservative justice
            theta_5_4 = [theta_5_4, theta(i,2,1)];
        else
            % liberal justice
            theta_5_4 = [theta_5_4, theta(i,2,2)];
        end
    end
end

A_5_4 = zeros(9);
for i = 1:9
    for j = 1:9
        if i == j
            A_5_4(i,j) = 0;
        elseif i > j
            A_5_4(i,j) = 0;
        else
            A_5_4(i,j) = theta_5_4(1);
            theta_5_4(1) = [];
        end
    end
end

A_5_4 = A_5_4 + A_5_4';
A_5_4_thresh = A_5_4;
A_5_4_thresh(A_5_4_thresh<.1) = 0;
G = graph(A_5_4_thresh , justices);

figure 
H = plot(G)%,'EdgeLabel',G.Edges.Weight)
[s,t] = find(A_5_4_thresh>0)
for k = 1:length(s)
    highlight(H,[s(k), t(k)],'LineWidth', .5 * (A_5_4_thresh(s(k),t(k))/0.1))
end
title('MRF via IPF - 5-4 Conservative Decision')

% Create a 4-5 liberal case
theta_4_5 = [];
for i = 1:size(theta,1)
    edge = E{i};
    if edge(1) < 5
        % conservative justice
        if edge(2) < 5
            % conservative justice
            theta_4_5 = [theta_4_5, theta(i,1,1)];
        else
            % liberal justice
            theta_4_5 = [theta_4_5, theta(i,1,2)];
        end
    else
        % liberal justice
        if edge(2) < 5
            % conservative justice
            theta_4_5 = [theta_4_5, theta(i,2,1)];
        else
            % liberal justice
            theta_4_5 = [theta_4_5, theta(i,2,2)];
        end
    end
end

A_4_5 = zeros(9);
for i = 1:9
    for j = 1:9
        if i == j
            A_4_5(i,j) = 0;
        elseif i > j
            A_4_5(i,j) = 0;
        else
            A_4_5(i,j) = theta_4_5(1);
            theta_4_5(1) = [];
        end
    end
end

A_4_5 = A_4_5 + A_4_5';
A_4_5_thresh = A_4_5;
A_4_5_thresh(A_4_5_thresh<.1) = 0;
G = graph(A_4_5_thresh , justices);

figure 
H = plot(G)%,'EdgeLabel',G.Edges.Weight)
[s,t] = find(A_4_5_thresh>0)
for k = 1:length(s)
    highlight(H,[s(k), t(k)],'LineWidth', .5 * (A_4_5_thresh(s(k),t(k))/0.1))
end
title('MRF via IPF - 5-4 Liberal Decision')

% Construct 9-0 liberal case
theta_0_9 = theta(:,2,2);
A_0_9 = zeros(9);
for i = 1:9
    for j = 1:9
        if i == j
            A_0_9(i,j) = 0;
        elseif i > j
            A_0_9(i,j) = 0;
        else
            A_0_9(i,j) = theta_0_9(1);
            theta_0_9(1) = [];
        end
    end
end

A_0_9 = A_0_9 + A_0_9';
A_0_9_thresh = A_0_9;
A_0_9_thresh(A_0_9_thresh<.1) = 0;
G = graph(A_0_9_thresh , justices);

figure
H = plot(G)%,'EdgeLabel',G.Edges.Weight)
[s,t] = find(A_0_9_thresh>0)
for k = 1:length(s)
    highlight(H,[s(k), t(k)],'LineWidth', .5 * (A_0_9_thresh(s(k),t(k))/0.1))
end
title('MRF via IPF - 9-0 Liberal Decision')



