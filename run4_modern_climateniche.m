addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% download the 2015 year data
% This section of data may take sometime to run. Be patient and wait for
% around 10 minutes. But beware. This section of data does not need to be
% run everytime.

% documentation and data can be found here https://chelsa-climate.org/downloads/
% this is downloaded by clicking the data and using wget url.
% documentation of data can be found here 
% https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf
% or searching in zotero "chelsa"

% % Base directory
% dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate/';
% 
% % List of variables
% variables = {'pr', 'tas', 'sfcWind', 'hurs'};
% 
% % Loop over each variable
% for v = 1:length(variables)
%     % Get the current variable
%     var = variables{v};
%     
%     % Loop over each month
%     for month = 1:12
%         % Construct the URL and filename
%         url = sprintf('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/%s/CHELSA_%s_%02d_2015_V.2.1.tif', var, var, month);
%         filename = sprintf('CHELSA_%s_%02d_2018_V.2.1.tif', var, month);
%         full_path = fullfile(dir_path, filename);
% 
%         % Download the file
%         websave(full_path, url);
%     end
% end
%% %%%%%
% load and rescale the 2015 climate data you just downloaded into 1800*3600
% this data has high resolution. (I think) it can be rescaled into higher
% resolution
% 
% this section of code can also take sometime to run.
% the outcome of this section is 4 new .mat files in the main directory
% named "aggregate_hurs.mat", "aggregate_pr.mat", etc.

% % Specify your directory
% dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/climate/';
% 
% % List of variables to process
% %variables = {'pr', 'tas', 'sfcWind', 'hurs'};
% variables = {'hurs'};
% 
% % Scaling factors for each variable
% scaling_factors = [100, 10, 1000, 100];
% offsets = [0, -273.15, 0, 0];
% 
% % Loop over all variables
% for i = 1:length(variables)
%     % Current variable and scaling factor
%     var = variables{i};
%     scaling_factor = scaling_factors(i);
%     offset = offsets(i);
%     
%     % Load the first file
%     filename = sprintf('CHELSA_%s_01_2018_V.2.1.tif', var);
%     full_path = fullfile(dir_path, filename);
%     
%     % Load the data and get the size
%     data = imread(full_path);
%     data = imresize(data, [1740, 3600]); % Resize data
%     [nrows, ncols] = size(data);
%     
%     % Initialize the 3D matrix with the correct size
%     data_all = zeros([nrows, ncols, 12]);
%     data_all(:,:,1) = data / scaling_factor + offset; % Add the first month to the matrix
%     
%     % Load the rest of the data
%    for month = 2:12 % tas, pr, sfcWind
%            % for month = [2,3,5:12] % hurs data
% 
%         % Construct the filename
%         filename = sprintf('CHELSA_%s_%02d_2018_V.2.1.tif', var, month);
%         full_path = fullfile(dir_path, filename);
%         
%         % Load the data into the matrix
%         data = imread(full_path);
%         data = imresize(data, [1740, 3600]); % Resize data
%         data_all(:, :, month) = data / scaling_factor + offset;
%     end
%       Northplaceholder=zeros(60,3600)
% 
%     % Compute the mean or sum across the 3rd dimension (i.e., across months)
%     if strcmp(var, 'pr') % If the variable is 'pr', sum up monthly precipitations
%         aggregate_data = [Northplaceholder;sum(data_all, 3)];
%      else % For all other variables, compute the mean
%         aggregate_data = [Northplaceholder;mean(data_all, 3)];
%     end
%     
%     % Save the variable to a .mat file
%     save(['/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_' var '.mat'], 'aggregate_data');
% end
%% Load the data
% if you are running the code not for the first time, you only need to
% start with this section
% this is the 2015 year climate data
% the humidity values seem to have issues
% this is the chelsa data for present climate
load('./grazingniche/matdata/aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('./grazingniche/matdata/aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('./grazingniche/matdata/aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
aggregate_sfcWind(aggregate_sfcWind==66)=0
load('./grazingniche/matdata/aggregate_tas.mat')
aggregate_tas=aggregate_data;
% %% let's see the raw data
% figure;
% t = tiledlayout(2, 2);
% 
% % Plot mean annual temperature
% nexttile;
% imagesc(aggregate_tas);
% colorbar;
% title('Mean Annual Temperature');
% 
% % Plot precipitation
% nexttile;
% imagesc(aggregate_pr);
% caxis([0 6000]);
% colorbar;
% title('Precipitation');
% 
% % Plot near surface wind speed
% nexttile;
% imagesc(aggregate_sfcWind);
% colorbar;
% caxis([0 20]);
% title('Near Surface Wind Speed');
% 
% % Plot relative humidity
% nexttile;
% imagesc(aggregate_hurs);
% colorbar;
% title('Relative Humidity');
% 
% % Adjust the spacing and title of the tiled layout
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% title(t, 'Climate Data');
% 
% % Save the figure as a PNG
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/subplot_2015_climate.png');

%% new land use data 
% in case I forget, this data is from the EarthEnv data, which seems to be
% not so widely used, therefore is subject to quite some problems. 
% Also, I think it is wise to unify current and future land use data. 
% different land use dataset could yield very different results. 
% Also it is worth noting that the main difference of the niche depends
% mostly on the land use data and not the climatology data...
%% the below processing of shrubs/herbaceous data can be neglected
%% just skip to the grassland_env data
% this section does not need to be run twice
% if you are running the code for a second time, please run the below section code
% and not this one
% please see the documentation of the data here: https://www.earthenv.org/landcover
% shrubs=readgeoraster('consensus_full_class_5.tif')
% (save('shrubs.mat', 'shrubs'))
% herbaceous=readgeoraster('consensus_full_class_6.tif')
% (save('herbaceous.mat', 'herbaceous'))
% cultivated=readgeoraster('consensus_full_class_7.tif')
% (save('cultivated.mat', 'cultivated'))
% builtup=readgeoraster('consensus_full_class_9.tif')
% (save('builtup.mat', 'builtup'))
% barren=readgeoraster('consensus_full_class_11.tif')
% (save('barren.mat', 'barren'))
%
% load('shrubs.mat')
% load('herbaceous.mat')
% load('cultivated.mat')
% load('builtup.mat')
% load('barren.mat')
% Here I am processing the new land use data. 
% This dataset only covers 56S-90N. I am adding the south grid cells
% herbaceous=double(herbaceous);
% shrubs=double(shrubs);
% cultivated=double(cultivated);
% builtup=double(builtup);
% barren=double(barren);
% 
% gap=zeros(4080,43200)

% herbaceous_full=[herbaceous;gap];
% save('./grazingniche/matdata/herbaceous_full.mat', 'herbaceous_full', '-v7.3');
% shrubs_full=[shrubs;gap];
% save('./grazingniche/matdata/shrubs_full.mat', 'shrubs_full', '-v7.3');
% cultivated_full=[cultivated;gap];
% save('./grazingniche/matdata/cultivated_full.mat', 'cultivated_full', '-v7.3');
% builtup_full=[builtup;gap];
% save('./grazingniche/matdata/builtup_full.mat', 'builtup_full', '-v7.3');
% barren_full=[barren;gap];
% save('./grazingniche/matdata/barren_full.mat', 'barren_full', '-v7.3');
%% load the new land use data
% Or else you can just download the ('grassland_env.mat') data from below
% section and ignore this section
% you can just load the data from here and neglect the above section
% this section takes around 30s to run
load('herbaceous_full.mat')
load('shrubs_full.mat')
% in case you need the other land use types
% load('cultivated_full.mat')
%load('builtup_full.mat')
load('barren_full.mat')
%% resize the loaded grassland maps (herbaceous+shrubs)
R_env = georefcells([-90,90],[-180,180],size(herbaceous_full));
[herbaceous_resized,R_resized] = georesize(herbaceous_full,R_env,1/12,"bilinear");
[shrubs_resized,R_resized] = georesize(shrubs_full,R_env,1/12,"bilinear");
[barren_resized,R_resized]=georesize(barren_full,R_env,1/12,"bilinear")
%[builtup_resized,R_resized]=georesize(builtup_full,R_env,1/12,"bilinear")
%% 
imagesc(barren_resized)
colorbar
%% obtaining a grassland map called grassland_env
grassland_env=shrubs_resized+herbaceous_resized;
save('./grazingniche/matdata/grassland_env.mat','grassland_env');

%% now let's decide a niche: the best grazing conditions for grazing.
load('grassland_env.mat')
imagesc(grassland_env)
colorbar
%% nice graph of grassland
% to run this you need to add the command file into directory
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[-180,180],size(grassland_env));
geoshow(flipud(grassland_env), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/grassland.svg');

%% let's see if the niche we decided is suitable.
%% first graph the niche range
% This graph is really nice! You can actually see the grazing belt in the
% first graph. It is a very good presentation. The rest of the graphs cut
% the belts further into grazing zones. "Thermalneutral belt" for grazing.
% I need a really good graphic presentation of this to support my findings 

% below see another section on bi-color graphs, cleaer presentation

% Create a figure with tiled subplots

figure;
t = tiledlayout(2, 2);

% Plot mean annual temperature
nexttile;
imagesc(aggregate_tas);
colorbar;
caxis([5 29]);
title('Mean Annual Temperature');

% Plot precipitation
nexttile;
imagesc(aggregate_pr);
caxis([324 2805]);
colorbar;
title('Precipitation');

% Plot near surface wind speed
nexttile;
imagesc(aggregate_sfcWind);
colorbar;
caxis([1.5 6.0]);
title('Near Surface Wind Speed');

% Plot relative humidity
nexttile;
imagesc(aggregate_hurs);
colorbar;
caxis([47 87]);
title('Relative Humidity');

% Adjust the spacing and title of the tiled layout
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t, 'Climate Data');

% Save the figure as a PNG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_range.png');
%%
% To change the color scheme so that values within a 
% specified range are shown in green and values outside 
% this range are shown in white, you can create a custom 
% colormap and use logical indexing.
figure;
t = tiledlayout(2, 2);

% Custom colormap
colormap_matrix = [1 1 1; 0 1 0]; % White and green

% Plot mean annual temperature
nexttile;
imagesc(aggregate_tas);
colormap(gca, colormap_matrix);
mask = aggregate_tas >= 5 & aggregate_tas <= 29;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Mean Annual Temperature');

% Plot precipitation
nexttile;
imagesc(aggregate_pr);
colormap(gca, colormap_matrix);
mask = aggregate_pr >= 324 & aggregate_pr <= 2805;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Precipitation');

% Plot near surface wind speed
nexttile;
imagesc(aggregate_sfcWind);
colormap(gca, colormap_matrix);
mask = aggregate_sfcWind >= 1.51 & aggregate_sfcWind <= 6;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Near Surface Wind Speed');

% Plot relative humidity
nexttile;
imagesc(aggregate_hurs);
colormap(gca, colormap_matrix);
mask = aggregate_hurs >= 47 & aggregate_hurs <= 87;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Relative Humidity');

% Adjust the spacing and title of the tiled layout
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t, 'Climate Data');

% Save the figure as a PNG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_range.png');

% using this graph I can now see that relative humidity and near surface
% wind area are not playing a role in determining the grazing niche at all.
%% let's couple the niche with the land use data
% in order to run this section you need to add the command file into
% directory in order for addcolorplus to work
% this section produces the main graph of the paper. A map of the grazing
% niche. 
% Here I am not using the nasa land use data but the Earthenv data. NASA is
% a catagorical land use dataset but Earthenv is a precentage map. It shows
% the percentage of shrublands and herbacesous vegetation cover. 

%cond1_pr = (aggregate_pr >= 433 & aggregate_pr <= 2368);
cond1_pr = (aggregate_pr >= 50 & aggregate_pr <= 2627);
cond1_tas = (aggregate_tas >= -3 & aggregate_tas <= 29);
cond1_sfcWind = (aggregate_sfcWind >= 1 & aggregate_sfcWind <= 6);
cond_hurs = (aggregate_hurs >= 39 & aggregate_hurs <= 67);
cond1_landuse = (grassland_env>0);

% Comerbine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% incase you want to see what is holding you back, you can use this code to
% see
%imagesc(cond1_pr)

% %here I am trying to see the numerical values of my resuls
% niche_gridnum=numel(find(niche1(:)==1));
% land_gridnum=numel(find(landuse(:)>0));
% pasture_gridnum=numel(find((landuse(:)>5 & landuse(:)<=10)));
% pasture_prct=pasture_gridnum/land_gridnum;
% niche_prct1=niche_gridnum/pasture_gridnum;

% Create the coupled landuse map
landuse_coupled1 = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,84],[-180,180],size(landuse_coupled1));

figure1 = figure('WindowState','fullscreen');
sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_general.svg');

save('./grazingniche/matdata/niche.mat',"landuse_coupled1")
