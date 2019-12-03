%% Age-Weight Driver
% This file explores the effects of age-weighting on Barabasi-Albert style network 
% generation. Age functions are considered as follows:

% PART 1: linear, power-law:  f(age) = 1/age^alpha
% PART 2: threshold:          f(age) = (1 if age<threshold, 0 if age>threshold)
% PART 3: normalized Poisson: f(age) = Poisson(age,lambda)./max(long Poisson)

%for all networks, the following will be found
% (i) power-law coeff of degree distribution
% (ii) mean shortest path
% (iii) assortativity coefficient
%% PART 1: none/linear/power-law age-weighting
% for a given set of parameters, each of the above will be explored. A distribution 
% will be generated for each parametrization 

%set parameters to be constant across all of part 1:
N=200;
M=N; %try round(N*1.25);
%% 1A: linear

rng(1)

%set parameters
nreps=1; %SET TO 10+ FOR FULL RUN

linear_nets={};
for net_idx=1:nreps
    %Initialize adjacency matrix at time t = 0.
    A = zeros(1,1);

    %Add a node at each timestep.
    for t = 1:N-1

        %Get a list of the degrees of each node (Sum the columns in A).
        deg = sum(A,1);
        %Get the age of each node.
        age = [t:-1:1];

        %Compute LINEAR age_weight. Weight ranges from 0 to 1.
        age_weight = arrayfun(@(x) 1-(x-1)/(N-1), age);
        %Compute degree_weight, normalized. Weight ranges from 0 to 1.
        deg_weight = arrayfun(@(x) x/max(deg),deg);
        deg_weight(isnan(deg_weight)) = 1;

        %Compute total weight of each node. Weight ranges from 0 to 1.
        weight = (age_weight.*deg_weight);

        %Choose the node that will gain a neighbor.
        r = rand*sum(weight);
        for n = 1:t
            if sum(weight(1:n)) >= r
                %Chooses node n.
                break;
            end
        end

        %Add the node to the adjacency matrix.
        A(n,size(A,2)+1) = 1;  
        A(size(A,1)+1,n) = 1;

    end

    linear_nets{net_idx}=A;
    G = graph(A);
    figure
    plot(G,'layout','force');
    title('linear')
end
%% 1B: constant + power-law

rng(1)

%set parameters
nreps=1; %SET TO 10 FOR FULL
alphas=[0:.5:2]; %SET TO 0:.1:2 FOR FULL RUN


powerlaw_nets={};

%loop through each alpha
for alpha_idx=1:length(alphas)
    alpha=alphas(alpha_idx);
    
    for net_idx=1:nreps
        %Initialize adjacency matrix at time t = 0.
        A = zeros(1,1);

        %Add a node at each timestep.
        for t = 1:N-1

            %Get a list of the degrees of each node (Sum the columns in A).
            deg = sum(A,1);
            %Get the age of each node.
            age = [t:-1:1];

            %Compute LINEAR age_weight. Weight ranges from 0 to 1.
            age_weight = arrayfun(@(x) 1/x^(-alpha), age);
            %Compute degree_weight, normalized. Weight ranges from 0 to 1.
            deg_weight = arrayfun(@(x) x/max(deg),deg);
            deg_weight(isnan(deg_weight)) = 1;

            %Compute total weight of each node. Weight ranges from 0 to 1.
            weight = (age_weight.*deg_weight);

            %Choose the node that will gain a neighbor.
            r = rand*sum(weight);
            for n = 1:t
                if sum(weight(1:n)) >= r
                    %Chooses node n.
                    break;
                end
            end

            %Add the node to the adjacency matrix.
            A(n,size(A,2)+1) = 1;  
            A(size(A,1)+1,n) = 1;

        end

        powerlaw_nets{alpha_idx,net_idx}=A;
% % %         G = graph(A);
% % %         figure
% % %         plot(G,'layout','force');
% % %         title(['Power-law: \alpha=-' num2str(alpha)])
    end
end
%% PART 2: AGE THRESHOLD
%% 2C: threshold

rng(1)

%set parameters
nreps=1; %SET TO 10 FOR FULL
thresholds=[5 25 100 150 200]; %SET TO 0:.1:2 FOR FULL RUN


threshold_nets={};

%loop through each alpha
for thresh_idx=1:length(thresholds)
    thresh=thresholds(thresh_idx);
    
    for net_idx=1:nreps
        %Initialize adjacency matrix at time t = 0.
        A = zeros(1,1);

        %Add a node at each timestep.
        for t = 1:N-1

            %Get a list of the degrees of each node (Sum the columns in A).
            deg = sum(A,1);
            %Get the age of each node.
            age = [t:-1:1];

            %Compute LINEAR age_weight. Weight ranges from 0 to 1.
            age_weight = arrayfun(@(x) x<=thresh, age);
            %Compute degree_weight, normalized. Weight ranges from 0 to 1.
            deg_weight = arrayfun(@(x) x/max(deg),deg);
            deg_weight(isnan(deg_weight)) = 1;

            %Compute total weight of each node. Weight ranges from 0 to 1.
            weight = (age_weight.*deg_weight);

            %Choose the node that will gain a neighbor.
            r = rand*sum(weight);
            for n = 1:t
                if sum(weight(1:n)) >= r
                    %Chooses node n.
                    break;
                end
            end

            %Add the node to the adjacency matrix.
            A(n,size(A,2)+1) = 1;  
            A(size(A,1)+1,n) = 1;

        end

        threshold_nets{thresh_idx,net_idx}=A;
% % %         G = graph(A);
% % %         figure
% % %         plot(G,'layout','force');
% % %         title(['Power-law: t_{threshold}=' num2str(thresh)])
    end
end
%% PART 3: NORMALIZED POISSON AGING

%TODO
%TODO
%TODO
%% PART 4: NETWORK STATISTICS

%Power-law statistics.
for k = 1:nreps
    for i=1:length(alphas)
        %Degree Distribution.
        G = graph(powerlaw_nets{i,1});
        D = degree(G);
        degDist = histcounts(D,max(D));
        X = [];
        Y = [];
        for j = 1:size(degDist,2)
            if not(degDist(j) == 0)
                X = [X,j];
                Y = [Y,degDist(j)];
            end
        end
        p = polyfit(log(X),log(Y),1);
        plExp = p(1);
        
        %Average Shortest Path.
        sPaths = distances(G);
        avgSPath = mean2(sPaths);
        
        %Assortativity Coefficient
        assort = assortativity(powerlaw_nets{i,k},0);
        
        stats = [plExp,avgSPath,assort];
        powerlaw_stats{i,k} = stats;
    end
end

%Threshold statistics.
for k = 1:nreps
    for i=1:length(thresholds)
        %Degree Distribution.
        G = graph(threshold_nets{i,1});
        D = degree(G);
        degDist = histcounts(D,max(D));
        X = [];
        Y = [];
        for j = 1:size(degDist,2)
            if not(degDist(j) == 0)
                X = [X,j];
                Y = [Y,degDist(j)];
            end
        end
        p = polyfit(log(X),log(Y),1);
        plExp = p(1);
        
        %Average Shortest Path.
        sPaths = distances(G);
        avgSPath = mean2(sPaths);
        
        %Assortativity Coefficient
        assort = assortativity(threshold_nets{i,k},0);
        
        stats = [plExp,avgSPath,assort];
        threshold_stats{i,k} = stats;
    end
end

% FOR ALL NETWORKS:
%    -find power-law exponent of degree distribution's CCDF
%    -find average shortest path between all nodes (using distances() func)
%    -find coefficient of assortativity

% plot these variables as a function of the parameters:
%    -for power-law f(age), x axis is alpha
%    -for threshold f(age), x axis is threshold
%    -for normalized poisson f(age), x axis is lambda
%% PLOT NETS SO FAR

%plot powerlaw
figure
len_a=length(alphas);
for i=1:len_a
    subplot(1,len_a,i)
    plot(graph(powerlaw_nets{i,1}),'layout','force');
    title(['Power-law: \alpha=-' num2str(alphas(i))])
end

%plot threshold
figure
len_th=length(thresholds);
for i=1:len_th
    subplot(1,len_th,i)
    plot(graph(threshold_nets{i,1}),'layout','force');
    title(['Threshold: t=' num2str(thresholds(i))])
end
%%
adj2gephilab('toGephi',A);
%%
function EdgeL=adj2gephilab(filename,ADJ,parameters)
    % Convert ana adjacency matrix of a graph to 2 spreadhseets csv files
    % one for the edge table and the other for node table.
    % The files _node.csv and _edge.csv have to be open 
    %in Gephi via Data Laboratory.
    % INPUTS:
    %           filename: string for the prefix name of the two files .csv
    %           ADJ: the adjacency matrix
    %           parameters: vector as properties of the node to use as
    %                              attributes of them.
    % OUTPUTS:
    %            two csv spreadsheet files: 
    %                       filename_node.csv
    %                       filename_edge.csv
    %             EdgeL = it returns the edge list corresponing to the
    %             adjacency matrix. (it can be saved to be open in Gephi too)
    %
    % The two files must be open in Gephi via the Data Laboratory
    nodecsv=[filename,'_node.csv'];
    edgecsv=[filename,'_edge.csv'];
    n=size(ADJ,1); % square adjacency matrix
    if nargin<3 
        parameters=ones(n,1);% all nodes have the same attributes
    end
    ps=parameters;
    % Node Table:
    % header for node csv
    fidN = fopen(nodecsv,'w','native','UTF-8'); % best format for gephi
    fprintf(fidN,'%s\n','Id;Label;Attribute');
    % 
    for i=2:n+1,
        fprintf(fidN,'%s\n',[ num2str(i-1) ';"Node ' num2str(i-1) '"'...
            ';' num2str(ps(i-1))]);
    end
    fclose(fidN);
    % Edge Table
    EdgeL=conv_EdgeList(ADJ);
    S=EdgeL(:,1); % sources
    T=EdgeL(:,2); % targets
    W = EdgeL(:,3); % weights
    fidE = fopen(edgecsv,'w','native','UTF-8');
    % header for edge csv
    fprintf(fidE,'%s\n','Source;Target;Label;Weight');
    for i=2:length(S)+1,
          fprintf(fidN,'%s\n',[ num2str(S(i-1)) ';' num2str(T(i-1)) ';'...
              '"Edge from ' num2str(S(i-1)) ' to ' num2str(T(i-1)) '"' ';'...
              num2str(W(i-1))]);
    end
     fclose(fidE);
end
%% Aux function

function EdgeL=conv_EdgeList(adj)
    % convert adj matrix to edge list
    n=size(adj,1); % number of nodes
    edges=find(adj>0); % indices of all edges
    n_e=length(edges);
    EdgeL=zeros(n_e,3);
    for e=1:n_e
      [i,j]=ind2sub([n,n],edges(e)); % node indices of edge e  
      EdgeL(e,:)=[i j adj(i,j)];
    end
end