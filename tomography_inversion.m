function [body_inversion] = tomography_inversion(G,d,prior,sgm_m,signal2noise)
% var(body_resized(:,125)) =  7.9893e-06;
L = sqrt(size(G,2));
nd = length(d);

sgm_d2 = var(d)/signal2noise;
%sgm_m2 = (5*sqrt(8e-06))^2;
sgm_m2 = sgm_m^2;

% for i = 1:L
%     for j = 1:L
%         sgm_m(i,j) = exp( -(  (i-L/2)^2 +(j-L/2)^2  )/(2*(L/4)^2) );
%     end
% end
%sgm_m2A  = sgm_m2*(reshape(sgm_m,L^2,1)).^2;
%C_m = sparse(1:L^2,1:L^2,sgm_m2A,L^2,L^2);

um = reshape(prior,L^2,1);
%uml = sgm_m2*G'*((sgm_m2*(G*G') + sgm_d2*eye(nd,nd))\d);
uml = um + G'*(((G*G') + (sgm_d2/sgm_m2)*speye(nd,nd))\(d-G*um));
%uml = C_m*G'*((G*C_m*G' + sgm_d2*eye(nd,nd))\d);

body_inversion =  reshape(uml,L,L);

% figure
% imagesc(body_inversion)
% caxis([0 0.01])
% colormap('bone')