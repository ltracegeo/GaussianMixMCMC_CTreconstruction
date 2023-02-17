
%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO... MAXIMO DE SIGNAL2NOISE 100


load('data_toy_model_berea.mat')

data_to_use = data_filtrado_toy;

signal2noise = 10;
%signal2noise = 1000;
image_size = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Define Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialize=1;
if initialize == 1
    fprintf('Inicializando ...  ')    
    body = imresize(double(data_to_use),image_size/size(data_to_use,1), 'nearest');
    %body = 0.01*body/max(max(body));
    
    %[ D2 ] = construct_diferential_matrix2d( image_size ,image_size );
    [p] = correlation_matrix_2d(image_size,image_size,2,2);
    
    %C_m_d2 = 10000*inv(D2'*D2 + 1e-5*speye(size(image_size*image_size)));   
    sgm_m = 100;
    C_m = p*sgm_m^2;
    
    initialize =0;
    fprintf('INICIADO ! \n')
end


[~, segmented_body] = bayesian_inference_1D_gau(body, PRIOR_KH);

figure
imagesc(segmented_body)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Full Angle Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data
D = image_size;
%dsensor2 = 1.5;
%inc = 1;
dsensor2 = 0.75;
inc = 2;
[F_full_noNoise, sensor_pos_full, fan_rot_angles_full] = fanbeam(body,D,'FanSensorSpacing',dsensor2,'FanRotationIncrement',inc);
noise01 = randn(size(F_full_noNoise(:)));
F_full = F_full_noNoise + reshape(noise01,size(F_full_noNoise)).*sqrt(var(F_full_noNoise)/signal2noise);


theta = sensor_pos_full(end) - sensor_pos_full(1);
number_sources = size(F_full,2);
number_dedector = size(F_full,1);

fprintf('Construindo G ...  ')
[d_full_noNoise,G_full] = simulate_tomography(body,number_sources,0,number_dedector,(theta*pi/180),D);
fprintf('FEITO ! \n')
d_full = d_full_noNoise + reshape(noise01,size(d_full_noNoise)).*sqrt(var(d_full_noNoise)/signal2noise);

%inversions
fprintf('Invertendo Back projection...  ')
inversion_ifan_full = ifanbeam(F_full,D,'FanSensorSpacing',dsensor2,'FanRotationIncrement',inc,'OutputSize',image_size);
fprintf('FEITO ! \n')
inversion_ifan_full_resize = imresize(inversion_ifan_full,[image_size image_size]);
[~, segmented_BackProj_full] = bayesian_inference_1D_gau(inversion_ifan_full_resize , PRIOR_KH);


fprintf('Invertendo BLI ...  ')
[inversion_BLI_full] = tomography_inversion(G_full,d_full,zeros(image_size,image_size),sgm_m,200);
fprintf('FEITO ! \n')
[~, segmented_BLI_full] = bayesian_inference_1D_gau(inversion_BLI_full, PRIOR_KH);

fprintf('Invertendo BLI corr...  ')
[inversion_BLIcorr_full] = tomography_inversion_TV2(G_full,d_full,zeros(image_size,image_size),C_m,200);
fprintf('FEITO ! \n')
[~, segmented_BLIcorr_full] = bayesian_inference_1D_gau(inversion_BLIcorr_full, PRIOR_KH);

% para imagesize 20 
%P = [ 0.5 0.25 0.25;
%      0.25 0.5 0.25;
%      0.25 0.25 0.5 ];
  
% para imagesize 50   
 P = [ 0.7 0.15 0.15;
       0.15 0.7 0.15;
       0.15 0.15 0.7];  
%   
% para imagesize 100 
% P = [ 0.9 0.05 0.05;
%       0.05 0.9 0.05;
%       0.05 0.05 0.9];

tic
fprintf('Invertendo Linear GaussMix MCMC ...  ')
subplot(3,1,1)
imagesc(segmented_body)
title('Reference')
[INVERSION] =  linear_Gaussian_mixture_MCMC(G_full,d_full,zeros(image_size,image_size),sgm_m,P,PRIOR_KH,signal2noise);
%[INVERSION] =  linear_Gaussian_mixture_MCMC_v2(G_full,d_full,zeros(image_size,image_size),C_m,P,PRIOR_KH,signal2noise);
fprintf('FEITO ! \n')
toc


INVERSION_ = INVERSION;

subplot(3,1,2)
imagesc(INVERSION_.CLASS.map)
title('Inverted')



figure
plot(INVERSION_.log_likelyhood)
grid



figure
subplot(2,3,1)
imagesc(body)
title('Reference')
caxis([0 200])
subplot(2,3,2)
imagesc(segmented_body)
title('Reference Segmented')
subplot(2,3,3)
imagesc(segmented_BackProj_full)
title('Back Projection')
subplot(2,3,4)
imagesc(segmented_BLI_full)
title('Inversion BLI white')
subplot(2,3,5)
imagesc(segmented_BLIcorr_full)
title('Inversion BLI Corr')
subplot(2,3,6)
imagesc(INVERSION_.CLASS.map)
title('Inversion GaussianMix MCMC ')



 figure
 subplot(2,3,1)
 imagesc(body)
 title('Reference')
 caxis([0 200])
 subplot(2,3,2)
 imagesc(segmented_body)
 title('Reference Segmented')
 colormap('bone')
 subplot(2,3,3)
 imagesc(inversion_ifan_full)
 title('Back Projection')
 caxis([0 200])
 colormap('bone')
 subplot(2,3,4)
 imagesc(inversion_BLI_full)
 title('Inversion BLI white')
 caxis([0 200])
 colormap('bone')
 subplot(2,3,5)
 imagesc(inversion_BLIcorr_full)
 title('Inversion BLI Corr')
 caxis([0 200])
 colormap('bone')
 subplot(2,3,6)
 imagesc(INVERSION_.ATENUATION.mean)    
 title('Inversion GaussianMix MCMC ')
 caxis([0 200])
 colormap('bone')



figure
histogram(body(segmented_body==1),'Normalization','pdf')
hold all
histogram(body(segmented_body==2),'Normalization','pdf')
histogram(body(segmented_body==3),'Normalization','pdf')

plot(axis,normpdf(axis, PRIOR_KH(1).MU, sqrt(PRIOR_KH(1).C(1,1)) ),'LineWidth',2)
plot(axis,normpdf(axis, PRIOR_KH(2).MU, sqrt(PRIOR_KH(2).C(1,1)) ),'LineWidth',2)
plot(axis,normpdf(axis, PRIOR_KH(3).MU, sqrt(PRIOR_KH(3).C(1,1)) ),'LineWidth',2)
grid