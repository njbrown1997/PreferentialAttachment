%Specify parameters.
M=130;
N=100;

p=(M-N)/N;
if p<0
    p=0;
end

timesteps = N; % This is also the number of nodes.

%Initialize adjacency matrix at time t = 0.
A = zeros(1,1);

%Add a node at each timestep.
for t = 1:timesteps-1

    %Get a list of the degrees of each node (Sum the columns in A).
    deg = sum(A,1);
    %Get the age of each node.
    age = [t:-1:1];
    
    %Compute age_weight based on some function. Weight ranges from 0 to 1.
    age_weight = arrayfun(@(x) 1/x, age);
    %Compute degree_weight based on some function. Weight ranges from 0 to 1.
    deg_weight = arrayfun(@(x) x/max(deg),deg);
    deg_weight(isnan(deg_weight)) = 1;
    
    %Compute total weight of each node. Weight ranges from 0 to 1.
    %There might be a better way to do this.
    weight = (age_weight + deg_weight)/2;
    
    %Choose the node that will gain a neighbor.
    r = rand*sum(weight);
    for n = 1:t
        if sum(weight(1:n)) >= r
            %Choose node n.
            break;
        end
    end

    %Add the node to the adjacency matrix.
    A(n,size(A,2)+1) = 1;  
    A(size(A,1)+1,n) = 1;
    
    if rand()<p
        [x,y]=find(A==0);
        if ~isempty(x)
            idx = datasample(1:length(x),1);
            
            % add new edge
            A(x(idx),y(idx))=1;
            A(y(idx),x(idx))=1;
        end
    end
end

G = graph(A);

disp(G.numnodes)
disp(G.numedges)
%
%MATLAB visualizations aren't that great, this is just for testing.
deg = sum(A,1);
age = [t+1:-1:1];
figure
plot(G,'layout','force');

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
