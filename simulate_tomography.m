function [d,G] = simulate_tomography(body,n_source,angle_source,n_detectors,angle_beam,R)
% ver se esta copiando um mais de uma vez !!!!!
L = size(body,1);
%R = sqrt(2)*L/2;
%R=(75/16)*L/2;
%R = 160;
data = zeros(n_detectors,n_source);

d_angle_source = 1*pi/n_source;
operator_line = 1;
aux=1;
for source=1:n_source   
    r_source = [R*cos(angle_source) R*sin(angle_source)];
    d_angle = 2*angle_beam/n_detectors;
    angle_detector = angle_source + pi - angle_beam+d_angle/2;
        for detector = 1:n_detectors

        r_detector =  [R*cos(angle_detector) R*sin(angle_detector)];
        angle_detector = angle_detector + d_angle;

        dl = (r_detector - r_source)/norm(r_detector - r_source);

            r_step = r_source;
            intensity = 0;
            for step=1:int32(2*R)
                % Caminha da fonte em direcao ao detector em passo unitario dl
                r_step = r_step + dl;
                r_stepInt = int64(r_step);
                % Verifica se ta dentro da rede quadrada                
                if r_stepInt(1)<L/2 && r_stepInt(1)>=-L/2
                    if r_stepInt(2)<L/2 && r_stepInt(2)>=-L/2
                        %transforma cordenada i,j em indice unico do vetor empilhado
                        i = r_stepInt(1)+L/2+1;
                        j = r_stepInt(2)+L/2+1;
                        % IF e ELSE para evitar somar duas vezes os mesmo sítio
                        if aux==1
                            %G(operator_line, (j-1)*L+i) = 1;
                            indexI(aux) = operator_line;
                            indexJ(aux) = (j-1)*L+i;
                            value(aux) = -1;
                            aux=aux+1;
                            intensity = intensity - body(i,j);
                        else
                            if (j-1)*L+i ~= indexJ(aux-1)
                            indexI(aux) = operator_line;
                            indexJ(aux) = (j-1)*L+i;
                            value(aux) = -1;
                            aux=aux+1;
                            intensity = intensity - body(i,j);
                            end
                        end
                    end
                end
            end
        data(detector,source) = intensity;
        operator_line = operator_line + 1;        
        end
    angle_source = angle_source - d_angle_source;
end
d = reshape(data,n_detectors*n_source,1);
G = sparse(double(indexI),double(indexJ),double(value),n_source*n_detectors,L*L);