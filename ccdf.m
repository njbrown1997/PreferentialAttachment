function [f,x]=ccdf(distribution)
    %[f,x]=ecdf(distribution);
    %f=1-f;
    n=size(distribution,1);
    x=1:n;
    f=zeros(1,n);
    for i=x
        f(i)=sum(distribution>=i)/200;
    end
end
