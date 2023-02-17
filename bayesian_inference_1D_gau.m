function [prob_map, facies_map] = bayesian_inference_1D_gau(Vp, PRIOR, prior_proportion)
%%%% TODO, ARRUMAR COLOCANDO TREND (MÉDIA MOVEL) 

prob_map = zeros(size(Vp,1)*size(Vp,2),length(PRIOR));

if nargin == 2
    prior_proportion = ones(length(PRIOR),1)/length(PRIOR);
end

for facie=1:length(PRIOR)
    if length(fieldnames(PRIOR))==2        
        prob_map(:,facie) = prior_proportion(facie)*mvnpdf(Vp(:),PRIOR(facie).MU(1)*ones(1,size(Vp(:),1))',PRIOR(facie).C);
    else
        prob_map(:,facie) = prior_proportion(facie)*mvnpdf(Vp(:),PRIOR(facie).MU(1)*ones(1,size(Vp(:),1))',PRIOR(facie).C);
%        prob_map(:,facie) = prior_proportion(facie)*mvnpdf([Vp(:)],[PRIOR(facie).VP.trend  PRIOR(facie).VS.trend PRIOR(facie).RHOB.trend ],PRIOR(facie).C);
    end
end

soma = sum(prob_map,2);
for facie=1:length(PRIOR)
    prob_map(:,facie) = prob_map(:,facie)./soma;
end

prob_map = reshape(prob_map,size(Vp,1),size(Vp,2),length(PRIOR));

prob_map = reshape(prob_map,size(prob_map,1)*size(prob_map,2),length(PRIOR));
    
[~, facies_map] = max(prob_map, [], 2);
    
prob_map = reshape(prob_map,size(Vp,1),size(Vp,2),length(PRIOR));
facies_map = reshape(facies_map,size(Vp,1),size(Vp,2));

