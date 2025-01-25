
addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% In this code script I am determining the grassland grazing capacity
load('grassland_env.mat')
% this livestockdensity.mat comes from the below processing of data
% load this into directory and you don't need to process the below n
% sections of code related to livestock. 
% this livestockdensity.mat data has gotten rid of the non-grassland
% livestocks,etc.
load('livestockdensity.mat')
% this AGB.mat data comes from the below section on grassland biomass.
load('AGB.mat')
load('aggregate_tas.mat')
load('niche.mat')% variable name: landuse_coupled1 (in run4)

load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
aggregate_sfcWind(aggregate_sfcWind>10)=NaN
%% determining modern livestock distribution niche
imagesc(aggregate_tas)
imagesc(resizecattle)
imagesc(aggregate_hurs)
colorbar

%%
nexttile;
hold on;
scatter(vector_pr, vector_cattle,'filled');
colormap(jet);
xlabel('Precipitation (mm)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Precipitation');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [300 2700]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

nexttile;
hold on;
scatter(vector_tas, vector_cattle,'filled');
colormap(jet);
xlabel('Mean Annual Temperature (^{o}C)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Temperature');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [5 29]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

nexttile;
hold on;
scatter(vector_hurs, vector_cattle,'filled');
colormap(jet);
xlabel('Relative Humidity (%)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Humidity');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [47 87]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

nexttile;
hold on;
scatter(vector_sfcWind, vector_cattle,'filled');
colormap(jet);
xlabel('Near Surface Windspeed (m/s)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Windspeed');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [1.5 6]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

set(gcf, 'Position',  [478,192,778,687])
%%
imagesc(aggregate_sfcWind)
colorbar
caxis([0,10])
%% read population density data
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");
resizedpop(resizedpop<0)=0
imagesc(resizedpop)
colorbar
%% 导入cattle density 并resize
% this section of code will take 30s to run
[cattle, Rcattle] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/cattle/Glb_Cattle_CC2006_AD.tif');

[allresizecattle,resizedRcattle] = georesize(cattle,Rcattle,1/12,"bilinear");
% 把过于大的cattledensity替换成最大值，否则最后地图的colorbar显示的色差太小（像是没有数据）
%cattle(cattle>250) = 250;
% all nan values are recorded as -1. this interferes with the calculation
allresizecattle(allresizecattle<0)=0

% this is the grassland cattle
%I want resizecattle to equal to allresizecattle where grassland_env==0
resizecattle = allresizecattle; % First, make allresizecattle a copy of resizecattle
resizecattle(resizecattle>250) = 0;
resizecattle(resizedpop >20) = 0; % Then, set values to 0 where grassland_env is not 0
% cattlenum=sum(sum(resizecattle.*worldarea));
% cattlenumall=sum(sum(allresizecattle.*worldarea));



imagesc(allresizecattle)
colorbar

imagesc(resizecattle)
colorbar
%% R is a very important variable for graphing. for each section of code (history/now/future) the R
% shold be the same. Here I am setting R for present evaluation. Note that
% this R can be changed. 
R=georefcells([-90,90],[-180,180],size(resizecattle));
%% better graph: cattle
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland cattle density', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar;
R=georefcells([-90,90],[-180,180],size(resizecattle));
geoshow(flipud(resizecattle), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/cattle_density.png');

%% 导入sheep density 并resize
[sheep, Rsheep] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/sheep/Glb_SHAD_2006.tif');
[allresizesheep,resizedRsheep] = georesize(sheep,Rsheep,1/12,"bilinear");
% 把过于大的sheepdensity替换成最大值，否则最后地图的colorbar显示的色差太小（像是没有数据）
allresizesheep(allresizesheep<0) =0;

% keeping only the grassland sheep
resizesheep=allresizesheep;
resizesheep(resizedpop >20)=0
% sheepsnum=sum(sum(resizesheep.*worldarea));
% sheepsnumall=sum(sum(allresizesheep.*worldarea));
%% graph: sheep
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland sheep density', 'FontSize', 16);
custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(resizesheep), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/sheep_density.png');

%% 导入goats density 并resize
[goats, Rgoats] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/goats/Glb_GTAD_2006.tif');
Rgoats=Rsheep;% somehow Rgoats is a empty array and does not contain georeferencing information
[allresizegoats,resizedRgoats] = georesize(goats,Rgoats,1/12,"bilinear");
% 把过于大的goatsdensity替换成最大值，否则最后地图的colorbar显示的色差太小（像是没有数据）
allresizegoats(allresizegoats<0) =0;
allresizegoats(allresizegoats>250) = 250;
% keeping only the grassland goats!
resizegoats=allresizegoats;
resizegoats(resizedpop >20)=0;
% goatsnum=sum(sum(resizegoats.*worldarea));
% goatsnumall=sum(sum(allresizegoats.*worldarea));

%% graph: goats
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland goats density', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(resizegoats), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/goats_density.png');
%% livestock unit
% Although I calculated this, I did not use this in my analysis
% I individually accounted for the feed of each of the animal types
% but If you want to use standardized animal unit you can use this
% equation. Note the parameters might be different. 
livestock=resizecattle+0.2*resizegoats+0.2*resizesheep
% %% livestock intensity graph
% figure1 = figure('WindowState','fullscreen');
% sgtitle('Grassland livestock density', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1;addcolorplus(341)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% geoshow(flipud(livestock), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/livestock_density.png');
%% put them in one big graph
% figure1 = figure('WindowState','fullscreen');
% sgtitle('Grassland Animal Density', 'FontSize', 16);

% Custom positions for each subplot

% % First map (Top-left)
% subplot(2,3,1);
% title('Grassland cattle density');
% worldmap('world');
% %setm(ax1, 'FFaceColor', [1 1 1]);
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% geoshow(flipud(resizecattle), R, 'DisplayType', 'texturemap');
% load coastlines;
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Second map (Top-middle)
% subplot(2,3,2);
% title('Grassland sheep density');
% worldmap('world');
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% geoshow(flipud(resizesheep), R, 'DisplayType', 'texturemap');
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Third map (Top-right)
% subplot(2,3,3);
% title('Grassland goats density');
% worldmap('world');
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% geoshow(flipud(resizegoats), R, 'DisplayType', 'texturemap');
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Fourth map (Bottom)
% subplot(2,3,[4,6]);
% title('Grassland livestock density');
% worldmap('world');
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% colorbar;
% geoshow(flipud(livestock), R, 'DisplayType', 'texturemap');
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Save the figure
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/manu/livestock_density_subplot.png');

%% The livestock data processing is now complete. 
% save the livestock data
save('./grazingniche/matdata/livestockdensity.mat','resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep','livestock')
load('livestockdensity.mat') %variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep')

%% Here I am calculating the grassland above ground biomass data, 
%% AGB calculation
%% calculate the fANNP by using mean annual temperature to determine what fraction of 
% biomass is aboveground
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
fANPP=0.171+0.0129*aggregate_tas
imagesc(fANPP)
colorbar
%% tree cover data 
% This data is 30arcseconds, around 1km, so takes some
% time to load. Do not run this section of code once the resizing is done.
%data=readgeoraster('gm_ve_v1.tif')
% xlimits = [-90,90];
% ylimits = [-180,180]
% dataR = georefcells(xlimits,ylimits,size(data));
% [treecover,treecoverR] = georesize(data,dataR,1/12,"bilinear");
% save('treecover.mat')
load('treecover.mat')
treecoverR=georefcells([-90,90],[-180,180],size(treecover));
treecover(treecover>100) = 0;

imagesc(treecover)
colorbar
treecover = double(treecover);
treecovermultiplier = 1 ./ exp(4.45521 .* 0.01 .* treecover);
imagesc(treecovermultiplier)
colorbar
%% npp data
% check out the nasa npp data here: https://neo.gsfc.nasa.gov/view.php?datasetId=MOD17A3H_Y_NPP
%npp data
%npp=readgeoraster('MOD17A3H_Y_NPP_2022-01-01_rgb_3600x1800.TIFF')
npp = csvread('MOD17A3H_Y_NPP_2022-01-01_gs_3600x1800.CSV', 0, 0);
%  reading starts from the second row (1 indicates the second row as MATLAB is 0-based indexing) and the first column (0 indicates the first column).
npp(npp>4000)=0
imagesc(npp)
colorbar
%% calculation of grassland AGB by using formula in literature
% calculation
% this section of code is strange in that it always exhibit resizesheep(grassland_env==0)=0
% npp is double when I first loaded into matlab, but somehow it becomes
% unit 8 after I got to this section of code, so I have to load npp again
% conversion factors range from 0.47 to 0.5

AGB=((npp.*fANPP)/0.5).*treecovermultiplier
AGB(AGB>1000)=0
imagesc(AGB)
colorbar
% why is there a line along the coastline?
save('./grazingniche/matdata/AGB.mat','AGB')

%% 最后展示的图：npp、MAT、fANPP、tree cover、AGB
%% 算livestock capacity
% note that the livestock distribution is the number of livestock per
% square kilometer, but the npp data is biomass per square meter, so a
% conversion factor of 0.00001 is needed, also the npp data is grams per
% square meter,so that would mean a conversion factor of 1000 is needed
% here I am assuming the average cattle weight is 1000kg,goat weight is
% 70kg,wheep weight is 100kg, they all consume 2.5% of their body weight
% daily, and graze 365 days, and times a unit conversion factor of 0.001

% so here I am calculating the a percentage of the 
% biomass needed by current grazing animals and the available grassland
% biomass.  ju
% now I am no longer using 1000,70,100. I am using animal unit weight (AU) is 455. Cattle is 1(455kg), goats and sheep are 0.2
%(91kg)
load('livestockdensity.mat') %variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep')
load('grassland_env.mat')
load('AGB.mat')
cattlegrass=resizecattle*455*0.03*365*0.001;
goatsgrass=resizegoats*91*0.03*365*0.001;
sheepgrass=resizesheep*91*0.03*365*0.001;

livestockgrass=cattlegrass+goatsgrass+sheepgrass;

livestockgrass(grassland_env==0)=0

livestockpercentage=(livestockgrass)./AGB

livestockpercentage(grassland_env<3)=0

livestockpercentage(livestockpercentage>1) = 1;
livestockpercentage(livestockpercentage<0) = 0;

save('./grazingniche/matdata/livestockpercentage.mat','livestockpercentage')

imagesc(livestockpercentage);
colorbar

imagesc(AGB);
colorbar
caxis([0,200])

imagesc(grassland_env==0);
colorbar


imagesc(livestockgrass);
colorbar


%% nice graph of the livestock capacity
% figure1 = figure('WindowState','fullscreen');
% sgtitle('grassland capacity exceeded', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(341)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(grassland_env));
% geoshow(flipud(livestockpercentage), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% % set(gcf,'renderer','painters');
% % it seems matlab is unable to produce real svg with this map
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/livestockcapacity.svg');
% % 从这张图可以看出来两级分化现象非常严重，不是没有使用就是使用过度。


%% how many is the exceeded grassland capacity outside of the grazing niche
% get a logic layer of overgraze
overgraze_logic = (livestockpercentage >= 0.65 & livestockpercentage <= 1);
niche_logic=(landuse_coupled1>0);

% I now want to get a 1800*3600 matrix of where the value of overgraze_logic is 1 where as the niche_logic value is 0
overgraze_outsideniche_logic = overgraze_logic & ~niche_logic;
imagesc(overgraze_outsideniche_logic)

count1=sum(sum(overgraze_outsideniche_logic));
count2=sum(sum(niche_logic));
percentage_outniche=count1./count2

overgraze_insideniche_logic = overgraze_logic & niche_logic;
imagesc(overgraze_insideniche_logic)



save('./grazingniche/matdata/overgraze_niche.mat','overgraze_insideniche_logic','overgraze_outsideniche_logic','overgraze_logic')

% %% graph of outside niche+overgraze and inside niche+overgraze
% 
% figure7 = figure('WindowState','fullscreen');
% sgtitle('Overgrazed grassland inside of the niche', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(333)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(grassland_env));
% geoshow(flipud(overgraze_insideniche_logic), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% 
% %saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/insidenicheover_graze.svg');
% %%
% figure8 = figure('WindowState','fullscreen');
% sgtitle('Overgrazed grassland outside of the niche', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(333)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(grassland_env));
% geoshow(flipud(overgraze_outsideniche_logic), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 


%saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/outsideniche_overgraze.svg');

%% grassland that are outside of the niche?
% what are the regions that are covered with herbaceous and shrubs
% vegetation but are unsuitable for grazing?
% so this means grassland_env is >3, but landuse_coupled1 is 0

outniche = (grassland_env > 3) & (landuse_coupled1 == 0); % This creates a logical matrix
outniche = double(outniche); % Converts logicals to doubles: 1 for true, 0 for false
outniche=outniche.*grassland_env
imagesc(outniche)
colorbar

% figure9 = figure('WindowState','fullscreen');
% sgtitle('Grassland outside of the niche', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(341)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(outniche), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/outniche.svg');

%%
save('./grazingniche/matdata/livestock2015.mat','livestockpercentage','outniche','livestock','resizecattle','allresizecattle', 'resizesheep','allresizesheep','resizegoats','allresizegoats','resizedpop');     
%% Here I am dealing with AGBC and BGBC data
% These data are at a resolution of 300m. The input of these two variables
% could take all night (no idea how long exactly). Do not run this section
% of code. The result is two matrixes of 52201*129600 big. I have saved
% the resized version.
% tic
% [AGBC, Ragbc] = readgeoraster('aboveground_biomass_carbon_2010.tif')
% [BGBC, Rbgbc] = readgeoraster('belowground_biomass_carbon_2010.tif')
% toc

% tic
% [resizedAGBC,resizedRagbc] = georesize(AGBC,Ragbc,1/36,"bilinear");
% [resizedBGBC,resizedRbgbc] = georesize(BGBC,Rbgbc,1/36,"bilinear");
% toc
% %%
% placeholder1=zeros(60,3600);
% placeholder2=zeros(290,3600);
% resizedAGBC_=[placeholder1;resizedAGBC;placeholder2]
% resizedBGBC_=[placeholder1;resizedBGBC;placeholder2]
% resizedRagbc_ = georefcells([-90,90],[-180,180],size(resizedAGBC_));
% resizedRbgbc_ = georefcells([-90,90],[-180,180],size(resizedBGBC_));
% %
% save('./grazingniche/matdata/AGBCBGBC.mat','resizedAGBC_','resizedBGBC_','resizedRagbc_','resizedRbgbc_');     
% load('AGBCBGBC.mat')

% AGBC=resizedAGBC_;
% BGBC=resizedBGBC_
% Ragbc=resizedRagbc_;
% Rbgbc=resizedRagbc_
% 
% save('./grazingniche/matdata/AGBCBGBC.mat','AGBC','BGBC','Ragbc','Rbgbc');     

% this is the final result
load('AGBCBGBC.mat')
%%

[SOC0_30_, Rsoc_] = readgeoraster('SOCS_0_30cm_year_2010AD_10km.tif')
SOC0_30_(SOC0_30_<0)=0

[SOC0_200_, Rsoc_] = readgeoraster('SOCS_0_200cm_year_2010AD_10km.tif')
SOC0_200_(SOC0_200_<0)=0

[SOC0_30,Rsoc] = georesize(SOC0_30_,Rsoc_,3600/4320,"bilinear");
[SOC0_200,Rsoc] = georesize(SOC0_200_,Rsoc_,3600/4320,"bilinear");

imagesc(SOC0_200)
caxis([0 550]);
colorbar

save('./grazingniche/matdata/SOC.mat','SOC0_30','SOC0_200','Rsoc');     


%% carbon loss due to grazing
load('AGBCBGBC.mat')
load('SOC.mat')
load("worldarea.mat")
load('overgraze_niche.mat')
SOC0_200=double(SOC0_200);
AGBC=double(AGBC);
BGBC=double(BGBC);

% the unit of SOC, AGBC and BGBC are all MgC/ha, MgC is megagram, which is
% 10^6 gram. Common units used for SOC also include PgC, petagram, which is
% 10^15 grams
% I need to convert MgC/ha to MgC/km^2 first, so this is why I multipled it
% by 100. dividing it by 100 is becasue the grassland_env is supposed to be
% the fraction of grassland per pixle, but the values are not in
% percentages.
% carbonloss_SOC_map=overgraze_logic.*SOC0_200.*0.5.*grassland_env./100.*worldarea.*100;
carbonloss_SOC_map=(overgraze_insideniche_logic.*SOC0_200.*0.1.*grassland_env./100.*worldarea.*100)+((overgraze_outsideniche_logic.*SOC0_200.*0.2.*grassland_env./100.*worldarea.*100));
carbonloss_AGBC_map=overgraze_logic.*AGBC.*grassland_env./100.*worldarea.*100;
carbonloss_BGBC_map=overgraze_logic.*BGBC.*grassland_env./100.*worldarea.*100;

carbonloss_map=carbonloss_SOC_map+carbonloss_AGBC_map+carbonloss_BGBC_map

imagesc(carbonloss_map)
colorbar

sum(sum(carbonloss_map))
sum(sum(carbonloss_SOC_map))
sum(sum(carbonloss_AGBC_map))
sum(sum(carbonloss_BGBC_map))


carbonloss_outniche=sum(carbonloss_map(landuse_coupled1>1))

imagesc(landuse_coupled1)
colorbar

imagesc(overgraze_logic)
colorbar

imagesc(carbonloss_SOC_map)
colorbar

imagesc(carbonloss_AGBC_map)
colorbar

imagesc(carbonloss_BGBC_map)
colorbar


carbonloss_SOC=sum(sum(carbonloss_SOC_map));
carbonloss_AGBC=sum(sum(carbonloss_AGBC_map));
carbonloss_BGBC=sum(sum(carbonloss_BGBC_map));

save('./grazingniche/matdata/carbonloss_map.mat','carbonloss_map','carbonloss_SOC','carbonloss_AGBC','carbonloss_BGBC');     
%% bar figure of carbon loss due to grazing
figure_carbonloss={'SOC loss',carbonloss_SOC;
    'AGBC loss',carbonloss_AGBC;
    'BGBC loss',carbonloss_BGBC}

figure_carbonloss_value=[carbonloss_SOC;carbonloss_AGBC;carbonloss_BGBC]
figure_carbonloss_caption={'SOC loss';'AGBC loss';'BGBC loss'}

bar(figure_carbonloss_value)
xticklabels(figure_carbonloss_caption)
set(gcf, 'Position',  [584,785,378,216])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/bar/carbonlossbar.svg');

%%
% we find carbon loss resulting from grazing amounted to 50PgC, of which
% 26PgC is attributed to above ground carbon loss, 15PgC is attributed to
% below-ground carbon loss, 5PgC is attributed to soil carbon loss. 
% Overgrazing leads to system damage of local vegetation and soil and in
% turn leads to loss of carbon. 
%全球草地生态系统的碳储量约为308PgC，约有92%储存在土壤中，不到10%储存在生物量当中。
%放牧活动影响草草地净初级生产力、土壤呼吸，从而影响植被和土壤的碳储量。过度放牧是导致草地失去碳储量的罪魁祸首。通过分析我们发现：
% Grazing activities affect the net primary productivity of 
% grassland and soil respiration, thereby affecting the carbon 
% storage of vegetation and soil. Overgrazing is a major culprit 
% in grasslands losing carbon storage. Through analysis we found
% we find carbon loss resulting from grazing amounted to 50PgC, of which
% 26PgC is attributed to above ground carbon loss, 15PgC is attributed to
% below-ground carbon loss, 5PgC is attributed to soil carbon loss. 
% Overgrazing leads to system damage of local vegetation and soil and in
% turn leads to loss of carbon. 



% %%
% figure10 = figure('WindowState','fullscreen');
% sgtitle('AGBC', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(302)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(AGBC), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/AGBC.svg');
% 
% figure11 = figure('WindowState','fullscreen');
% sgtitle('BGBC', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(300)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(BGBC), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/BGBC.svg');
% 
% figure12 = figure('WindowState','fullscreen');
% sgtitle('SOC (0-2m)', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(301)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(SOC0_200), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SOC.svg');
% 
% figure13 = figure('WindowState','fullscreen');
% sgtitle('Tree cover map', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(291)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(treecover), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/treecover.svg');
% 
figure14 = figure('WindowState','fullscreen');
sgtitle('Net primary productivity', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(284)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[-180,180],size(npp));
geoshow(flipud(npp), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/NPP.svg');
% 
% figure15 = figure('WindowState','fullscreen');
% sgtitle('Map of carbon loss due to overgrazing', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(288)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(carbonloss_map), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/carbonloss.svg');
% 
% 
