

load('berea_raw_uint8.mat')
tomo = single(im8);

[X,Y,Z] = meshgrid([1:1024],[1:1024],[1:1024]);
order = 2;
filter = exp( -( (X(:)-512).^order + (Y(:)-512).^order + (Z(:)-512).^order ) / (2.^order) );
filter = reshape(filter,1024,1024,1024);
filter = single(filter/sum(filter(:)));
filtrado = fftshift(ifftn( fftn(filter).*fftn(tomo) ));