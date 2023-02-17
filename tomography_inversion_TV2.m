function [body_inversion] = tomography_inversion_TV2(G,d,prior,C_m,signal2noise)
%Cm = D2'*D2;
% var(body_resized(:,125)) =  7.9893e-06;
L = sqrt(size(G,2));
nd = length(d);
%sgm_d2 = (1*std(highPassFilter2(d,0.004,110,7)))^2;
sgm_d2 = var(d)/signal2noise;
%sgm_m2 = (5*sqrt(8e-06))^2;
%sgm_m2 = sgm_m^2;

%C_m = sgm_m2*inv(D2'*D2);

um = reshape(prior,L^2,1);
%uml = C_m*G'*((G*C_m*G' + sgm_d2*eye(nd,nd))\d);
uml = um + C_m*G'*((G*C_m*G' + sgm_d2*speye(nd,nd))\(d-G*um));
%uml = um + G'*(((G*G') + (sgm_d2/sgm_m2)*sparse(eye(nd,nd)))\(d-G*um));

body_inversion =  reshape(uml,L,L);

% figure
% imagesc(body_inversion)
%caxis([0 0.01])
%colormap('bone')