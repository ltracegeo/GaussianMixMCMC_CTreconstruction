

%% LOAD AND FILTER DATA
load('berea_raw_uint8.mat')
tomo = single(im8);

cut_data = tomo(200:end-300,200:end-200,565);
cut_data = tomo(1:end,1:end,565);

I = size(cut_data,1);
J = size(cut_data,2);

[X,Y] = meshgrid([1:J],[1:I]);
order = 2;
L=2;
filter = exp( -( (X(:)-I/2).^order + (Y(:)-J/2).^order ) / (L.^order) );
filter = fftshift(reshape(filter,size(cut_data,1),size(cut_data,2)));
filter = single(filter/sum(filter(:)));
filtrado = fftshift(ifftn( fftn(filter).*fftn(fftshift(cut_data)) ));


%% SEGMENTATION K means
n_facies=3;

segmented_km = kmeans(filtrado(:),n_facies);
segmented_km = reshape(segmented_km, size(cut_data));

figure
histogram(cut_data(:),'Normalization','pdf')
hold all
axis_x = linspace(min(cut_data(:)),max(cut_data(:)),100);
sum_pdf = 0;
for class=1:length(n_facies)
    PRIOR_KM(class).MU = mean(cut_data(segmented_km==class));
    PRIOR_KM(class).C(class) = var(cut_data(segmented_km==class));
    plot(axis_x,0.33*normpdf(axis_x,PRIOR_KM(class).MU,sqrt(PRIOR_KM(class).C(class))))
    sum_pdf = sum_pdf + 0.33*normpdf(axis_x,PRIOR_KM(class).MU,sqrt(PRIOR_KM(class).C(class)));
    %histogram(cut_data(segmented_km==class),'Normalization','pdf')
end
plot(axis_x,sum_pdf)

%% PRIOR NA MÃO
PRIOR_KH(1).MU = 108;
PRIOR_KH(1).C = 11^2;
PRIOR_KH(1).W = 0.2;
PRIOR_KH(2).MU = 138;
PRIOR_KH(2).C = 8^2;
PRIOR_KH(2).W = 0.7;
PRIOR_KH(3).MU = 190;
PRIOR_KH(3).C = 11^2;
PRIOR_KH(3).W = 0.1;

figure
histogram(cut_data(:),'Normalization','pdf')
hold all
axis_x = linspace(min(cut_data(:)),max(cut_data(:)),100);
sum_pdf = 0;
for class=1:length(PRIOR_KH)
    plot(axis_x,PRIOR_KH(class).W*normpdf(axis_x,PRIOR_KH(class).MU,sqrt(PRIOR_KH(class).C)))
    sum_pdf = sum_pdf + PRIOR_KH(class).W*normpdf(axis_x,PRIOR_KH(class).MU,sqrt(PRIOR_KH(class).C));   
end
plot(axis_x,sum_pdf)

[~, segmented_KH] = bayesian_inference_1D_gau(filtrado, PRIOR_KH);


%% Cut to toy model

%% DISPLAY
figure
subplot(2,2,1)
imagesc(filtrado) 
subplot(2,2,2)
imagesc(cut_data)
subplot(2,2,3)
imagesc(segmented_km)
subplot(2,2,4)
imagesc(segmented_KH)


I_cut = 653:872;
%I_cut = 262:445;
J_cut = 346:565;
%J_cut = 295:497;

data_toy = cut_data(I_cut,J_cut);
data_filtrado_toy = cut_data(I_cut,J_cut);
segmented_KH_toy = segmented_KH(I_cut,J_cut);


figure
subplot(2,2,1)
imagesc(data_toy )
subplot(2,2,2)
%imagesc(tomo(262:445,295:497,900))
imagesc(data_filtrado_toy )
subplot(2,2,3)
imagesc(segmented_km(653:872,346:565))
subplot(2,2,4)
imagesc(segmented_KH_toy)









