%Specify parameters.
timesteps = 10; % This is also the number of nodes.

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
end

%MATLAB visualizations aren't that great, this is just for testing.
deg = sum(A,1)
age = [t+1:-1:1]
G = graph(A);
plot(G);