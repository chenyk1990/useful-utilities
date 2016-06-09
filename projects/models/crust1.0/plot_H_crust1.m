% Apr. 17, 2015, Yunfeng Chen, plot the Moho variation of crust 1.0 model
lat = 48:0.5:58;
lon = -121:0.5:-108;
moho = zeros(length(lat),length(lon));
elevation = zeros(length(lat),length(lon));
tic
for i = 1:length(lat)
    for j = 1:length(lon)
        [~, elevation(i,j),moho(i,j)] = get_crust1_model(lat(i),lon(j),'rocks');
    end
end
toc
% mesh grid to create the corrdinates
[X,Y] = meshgrid(lon,lat);
surf(X,Y,-moho);hold on;
surf(X,Y,elevation*10);
% save
outdir = pwd;
m = [X(:),Y(:),moho(:),elevation(:)];
save([outdir,'/crust1.0.model'],'m','-ascii');
