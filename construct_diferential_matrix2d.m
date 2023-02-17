function [ D2 ] = construct_diferential_matrix2d( I,J )

D2 = sparse((I-1)*J,I*J);

aux=1;
aux1=0;
for j = 1:J
    for i=1:I-1
        D2(aux,aux+aux1) = -2;
        D2(aux,aux+aux1+1) = 1;
        if aux+aux1+I<=I*J
            D2(aux,aux+aux1+I) = 1; 
        end
        if aux+aux1>I*J-I
            D2(aux,aux+aux1) = -1;
        end
        aux=aux+1;
    end
    aux1 = aux1 + 1;
end
for j = 1:J-1
    D2(aux,j*I) = -1;
    D2(aux,j*I+I) = 1;
    aux=aux+1;
end

end

