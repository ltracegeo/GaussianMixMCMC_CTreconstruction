function [INVERSION] =  linear_Gaussian_mixture_MCMC(G,d,prior,sgm_m,P,PRIOR,signal2noise)

% TODO:
% 1) ATUALIZACAO GCmG' - Feito, não há vantagem nenhuma
% 2) corrigir signatl to noise - tentei, não deu
% 3) acceptance rate - feito
% 4) Random path completo - feito
% 5) Atualização em janela - Será que vale?
% 4) Tempratura? Annealing? ver JIFI??? nao sei o que


I = size(prior,1);
J = size(prior,1);
nd = length(d);
n_facies = size(P,1);


probability_map = zeros(I,J,n_facies);


P = [zeros(n_facies,1) P];
P = [0 ones(1,n_facies)/n_facies ;P];

%sgm_d2 = (1*std(highPassFilter2(d,4,110,7)))^2;
sgm_d2 = var(d)/signal2noise;
%sgm_m2 = (5*sqrt(8e-06))^2;
sgm_m2 = sgm_m^2;

posterior_mean = zeros(I,J);

model = ones(I+2,J+2);
model(2:end-1,2:end-1) = round(rand(I,J)*n_facies+0.5)+1;
model_atenuationMean = construct_mean(model(2:end-1,2:end-1)-1,PRIOR);
log_likelyhood = sum( (d-G*model_atenuationMean (:)).^2 )./sgm_d2;

random_path = randperm(J*I);
[X,Y] = meshgrid(1:I,1:J);

X = X(:);
Y = Y(:);

aux=0;
n_acceptance = 0;
n_step = 130;
step_convergence = 30;
log_likelyhood_array = zeros(n_step,1);
log_likelyhood_array(1) = log_likelyhood;
for step = 1:n_step
    %for up = 1:J*I
    %i = unidrnd(I) + 1;
    %j = unidrnd(J) + 1;
    for position = 1:length(random_path)
        i = X(random_path(position))+1;
        j = Y(random_path(position))+1;                
                
        probabilities = ones(1,n_facies);       
        Probs =  P(model(i-1,j),:).*P(model(i+1,j),:).*P(model(i,j-1),:).*P(model(i,j+1),:);                                      
        
        Probs = Probs./sum(Probs);
        class_sort = find( rand < cumsum(Probs));
        
        model_proposal = model;
        model_proposal(i,j) = class_sort(1);
        model_atenuationMean_proposal = construct_mean(model_proposal(2:end-1,2:end-1)-1,PRIOR);
        
        log_likelyhood_proposal = sum( (d-G*model_atenuationMean_proposal(:)).^2 )./sgm_d2;        
        
        if log_likelyhood_proposal < log_likelyhood
            model = model_proposal;
            model_atenuationMean = model_atenuationMean_proposal;
            log_likelyhood = log_likelyhood_proposal ;
            n_acceptance = n_acceptance + 1;
        else
            if rand<exp(-(log_likelyhood_proposal-log_likelyhood) )
                model = model_proposal;
                model_atenuationMean = model_atenuationMean_proposal;
                log_likelyhood = log_likelyhood_proposal ;
                n_acceptance = n_acceptance + 1;
            end
        end
        
    end
    
    
    if step>step_convergence
        aux = aux+1
        posterior_mean = posterior_mean + model_atenuationMean./(n_step-step_convergence);
        for class = 1:n_facies
            indicator_aux = zeros(I,J);
            indicator_aux((model(2:end-1,2:end-1)-1)==class) = 1;
            probability_map(:,:,class) = probability_map(:,:,class) + indicator_aux./(n_step-step_convergence);
        end
    end    
    
    
    subplot(3,1,2)
    imagesc(probability_map(:,:,2))
    title(step/n_step)
    drawnow
    subplot(3,1,3)
    imagesc(model(2:end-1,2:end-1))
    drawnow
    
    log_likelyhood_array(step) = log_likelyhood;
end

probability_map = reshape(probability_map,I*J,n_facies);

[~,MAP_class_model] = max(probability_map,[],2);

probability_map = reshape(probability_map,I,J,n_facies);
MAP_class_model = reshape(MAP_class_model,I,J);

INVERSION.CLASS.map = MAP_class_model;
INVERSION.CLASS.prob =  probability_map;
INVERSION.ATENUATION.mean = posterior_mean;
INVERSION.log_likelyhood = log_likelyhood_array;

INVERSION.acceptance_rate = n_acceptance/(n_step*length(random_path));

end


function [model_atenuationMean] = construct_mean(model,PRIOR)

model_atenuationMean = zeros(size(model));

for class = 1:length(PRIOR)
    model_atenuationMean(model==class) = PRIOR(class).MU;
end

end

