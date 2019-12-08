function A = build_network(N,age_fcn)
    A=zeros(1,1);
    
    for t = 1:N-1 %add each node
        
        %Get a list of the degrees of each node (Sum the columns in A).
        deg = sum(A,1);
        %Get the age of each node.
        age = t:-1:1;

        %Compute LINEAR age_weight. Weight ranges from 0 to 1.
        age_weight = arrayfun(age_fcn, age);
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

end