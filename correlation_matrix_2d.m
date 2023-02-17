function [p] = correlation_matrix_2d(I,J,dVert,dHor)

%t = cputime;
nPoints=I*J;

p = sparse(nPoints,nPoints);

parfor m = 1:nPoints
    %porcentagem = 100*(m/nPoints)
       
    
    j = double(int32(((m-1)/I)- 0.5+1));
    i = double(int32(((m/I) - double(int32((m/I)-0.5)))*I));
    if i == 0 
        i=I; 
    end
    
    for n =1:nPoints
        
        k = double(int32(((n-1)/I)- 0.5+1));
        l = double(int32(((n/I) - double(int32((n/I)-0.5)))*I));
        if l == 0 
            l=I; 
        end

        %p(m,n) = exp(-((i-l)/(dVert))^2)*exp(-((j-k)/(dHor))^2);
        value = exp(-sqrt(((i-l)/(dVert))^2))*exp(-sqrt(((j-k)/(dHor))^2));
        if value >0.005
            p(m,n) = value;
        end
    end
    
end

%tempo = cputime - t
