%% This script calculates the turnover impact on population and livestock figure

%% sRead and process 2100 population
% % read ncinfo
% pop_data_2100_info=ncinfo("/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/populationdensity/NetCDF/SSP3_NetCDF/total/NetCDF/ssp3_2100.nc")
% openvar("pop_data_2100_info")
% % read variable
% pop_data_2100=ncread("/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/populationdensity/NetCDF/SSP3_NetCDF/total/NetCDF/ssp3_2100.nc","ssp3_2100")
% openvar("pop_data_2100")
% % process variable so that lat and lon fill -90-90, -180-180
% pop_data_2100_1=pop_data_2100'
% imagesc(pop_data_2100_1)
% NorthZeros=zeros(53,2880)
% SouthZeros=zeros(270,2880)
% pop_data_2100_2=[NorthZeros;pop_data_2100_1;SouthZeros]
% imagesc(pop_data_2100_2)
% % get figure of population in 2100 to check
% figure9 = figure;
% custom_colormap = [1 1 1; addcolorplus(341)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% set(gcf, 'Colormap', custom_colormap);
% c = colorbar;  
% c.Ruler.TickLabelFormat='%g%%';
% R = georefcells([-90,90],[-180,180],size(pop_data_2100_2));
% geoshow(flipud(pop_data_2100_2), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% % process population in 2100, save. 
% pop_data_2100_2(pop_data_2100_2==2147483647)=0
% sum(sum(pop_data_2100_2)) % check sum in all pixel to check
% pop_count_2100=pop_data_2100_2
% save('./grazingniche/matdata/pop_count_2100.mat','pop_count_2100')
%% Read and process population in 2015.

% population count in 2015 (population per pixel)
[pop_count_2015, Rpop] = readgeoraster('gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_15_min.tif')
pop_count_2015(pop_count_2015<0)=0

sum(sum(pop_count_2015)) % check to see if total population match global population in 2015 (7.3 billion)
%% load livestock, grassland, and population data
% step 1: data preparation
load("grassland_env.mat")
load("livestockdensity.mat")
load('turnover_widest.mat')% this mat file includes these variables:compare_turnover_rcp26-compare_turnover_rcp85
load("pop_count_2100.mat")

% population per square kilometer
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
pop(pop<0)=0
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");

%% rescale population in 2100
% Step 2: Calculate the scaling factor to preserve the total population
Rpop_count_2100 = georefcells([-90,90],[-180,180],size(pop_count_2100));
[resizedpop_count_2100, resizedRpop] = georesize(pop_count_2100, Rpop_count_2100, 3600/2880, "bilinear");
original_sum = sum(pop_count_2100(:));  % Sum of the original data
resized_sum = sum(resizedpop_count_2100(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resizedpop_count_2100 = double(resizedpop_count_2100 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resizedpop_count_2100(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

sum(sum(resizedpop_count_2100))
sum(resizedpop_count_2100(:))
%% rescale population in 2015
% Step 2: Calculate the scaling factor to preserve the total population
Rpop_count_2015 = georefcells([-90,90],[-180,180],size(pop_count_2015));
[resizedpop_count_2015, resizedRpop] = georesize(pop_count_2015, Rpop_count_2015, 3600/1440, "bilinear");
original_sum = sum(pop_count_2015(:));  % Sum of the original data
resized_sum = sum(resizedpop_count_2015(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resizedpop_count_2015 = double(resizedpop_count_2015 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resizedpop_count_2015(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

sum(sum(resizedpop_count_2015))
sum(resizedpop_count_2015(:))
%% get the growth rate for population from 2015 to 2100 to scale future livestock
growth_pop=resizedpop_count_2100./resizedpop_count_2015
growth_pop(growth_pop>10)=10 % population projected to grow more than 10 fold per grid cell is unlikely and is due to extremely small values in 2015
% check if pop data and livestock data are aligned
nexttile
imagesc(growth_pop)
colorbar
caxis([0,10])
nexttile
imagesc(resizecattle)
% rescale cattle, sheep, and goats data
resizecattle_2100=resizecattle.*growth_pop
resizesheep_2100=resizesheep.*growth_pop
resizegoats_2100=resizegoats.*growth_pop
% graph future cattle, sheep, and goats
nexttile
imagesc(resizecattle_2100)
colorbar
caxis([0,250])
nexttile
imagesc(resizesheep_2100)
colorbar
caxis([0,250])
nexttile
imagesc(resizegoats_2100)
colorbar
caxis([0,250])
%% flip variables and check if they are aligned
turnover85_resize_widest=imresize(compare_turnover_rcp85_widest,[1800,3600])
imagesc(turnover85_resize_widest)
turnover85_resize_flip = double(turnover85_resize_widest(:, [ceil(end/2+1):end, 1:floor(end/2)]));

% step 2: check if all variables are aligned
nexttile
imagesc(turnover85_resize_flip)
nexttile
imagesc(resizedpop_count_2100)
nexttile
imagesc(resizecattle)
colorbar
%% get map of present population, future population, future cattle, future goats, future sheep
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(resizedpop_count_2100));
geoshow(flipud(resizedpop_count_2100), R, 'DisplayType', 'texturemap');   
caxis([0,10000])
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_2100.svg');
%% present pop
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(resizedpop_count_2100));
geoshow(flipud(resizedpop_count_2015), R, 'DisplayType', 'texturemap');   
caxis([0,10000])
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_2015.svg');
%% pop growth
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(growth_pop));
geoshow(flipud(growth_pop), R, 'DisplayType', 'texturemap');   
caxis([0,10])
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_growth.svg');
%% future cattle
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(resizecattle_2100));
geoshow(flipud(resizecattle_2100), R, 'DisplayType', 'texturemap');   
caxis([0,250])
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/resizedcattle_2100.svg');
%% future sheep
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(resizesheep_2100));
geoshow(flipud(resizesheep_2100), R, 'DisplayType', 'texturemap');   
caxis([0,250])
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/resizedsheep_2100.svg');
%% future goats
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(resizegoats_2100));
geoshow(flipud(resizegoats_2100), R, 'DisplayType', 'texturemap');   
caxis([0,250])
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/resizedgoats_2100.svg');

%% start calculation: population
% This calculates the population that is inhabited on grasslands
grassland_pop=zeros(size(resizedpop_count_2100))
grassland_pop(livestock>0&resizedpop<20&grassland_env>50)=resizedpop_count_2100(livestock>0&resizedpop<20&grassland_env>50)

grassland_pop_sum=sum(sum(grassland_pop))
% there are 0.16 billion grassland-based population
% how to define pastoralists? 
% pastoralists
%168,417,460

imagesc(grassland_pop)
colorbar

% This calculates the population that is positively influenced by GN shift
grassland_pop_positive=zeros(size(resizedpop_count_2100))
grassland_pop_positive(turnover85_resize_flip>0&grassland_env>50)=grassland_pop(turnover85_resize_flip>0&grassland_env>50)
grassland_pop_positive_sum=sum(sum(grassland_pop_positive))

% 22 million population will be positively impacted

grassland_pop_negative=zeros(size(resizedpop_count_2100))
grassland_pop_negative(turnover85_resize_flip<0&grassland_env>50)=grassland_pop(turnover85_resize_flip<0&grassland_env>50)
grassland_pop_negative_sum=sum(sum(grassland_pop_negative))

% 0.14 billion population will be negatively impacted
% 140 million will be negatively impacted

% of the negatively impacted 140 million populations, 85% are in
% middle-to-low income countries, and 50% in low-income countries, in extreme poverty countries 


disp(['grassland_pop_positive_sum:',num2str(grassland_pop_positive_sum)])
disp(['grassland_pop_negative_sum:',num2str(grassland_pop_negative_sum)])

% In total 1.3 billion population will be influenced, of which 1 billion is negatively influecned, and 0.2 billion is positivel influenced.

% In total 150 million population will be influenced, 80 million will be
% negatively influenced, 50 million will be positively
% influenced

% 3 billion population will be influenced. 1.7 will be negatively
% influenced, 1.2 billion will be positively influenced

% 0.48 billion population will be influenced. 0.18 will be positively
% influenced, 0.27 billion population will be negatively influenced

%% produce graph: negatively impacted population
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(grassland_pop_negative));
geoshow(flipud(grassland_pop_negative), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,100])
set(gcf, 'Position',  [361,651,489,258])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_negative_2100.svg');
%% produce graph: positively impacted population
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(287)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(grassland_pop_negative));
geoshow(flipud(grassland_pop_positive), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,100])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_positive_2100.svg');
%% Get country information-Get country mask
% Define path to your file
mapsize=[160,320]
x=ncinfo('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc')
file_path = '/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc';
countries = {x.Variables(3:217).Name};

% Create a structure to hold all country data
country_data = struct();

% Read lat and lon, and create a meshgrid
lon = ncread(file_path, 'lon');
lat = ncread(file_path, 'lat');
[lon, lat] = meshgrid(lon, lat);
Rmask = georefcells([-90,90],[-180,180], size(ncread(file_path, 'm_world')));
 
% Loop through countries, read and resize each mask
for i = 1:length(countries)
    country_name = countries{i};
    
    % Read country data
    country_data_raw = ncread(file_path, country_name);
    
    % Resize the country data
    [resized_data, ~] = georesize(country_data_raw, Rmask, mapsize(2)/720, "bilinear");
    
    % Store in the structure
    country_data.(country_name) = resized_data';
end
%% rescale grassland influenced population to match country mask (country mask is 160*320)
%% Positive rescale
resized_grassland_pop_positive= imresize(grassland_pop_positive,320/3600);
original_sum = sum(grassland_pop_positive(:));  % Sum of the original data
resized_sum = sum(resized_grassland_pop_positive(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_grassland_pop_positive = double(resized_grassland_pop_positive * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_grassland_pop_positive(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale
resized_grassland_pop_negative= imresize(grassland_pop_negative,320/3600);
original_sum = sum(grassland_pop_negative(:));  % Sum of the original data
resized_sum = sum(resized_grassland_pop_negative(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_grassland_pop_negative = double(resized_grassland_pop_negative * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_grassland_pop_negative(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% Get for all countries for negative
% Initialize an empty cell array to store results
pop_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in pop_neg
    pop_neg{i} = sum(sum(resized_grassland_pop_negative .* double(country_data.(fieldName))));
end
pop_neg_values = cell2mat(pop_neg)';

%% Get for all countries for positive
% Initialize an empty cell array to store results
pop_posi = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in pop_posi
    pop_posi{i} = sum(sum(resized_grassland_pop_positive .* double(country_data.(fieldName))));
end
pop_posi_values = cell2mat(pop_posi)';
openvar('pop_posi_values')
openvar('pop_neg_values')

%% %%%%%%%%% cattle

%% start calculation
load('worldarea.mat')
resizecattle_2100(isnan(resizecattle_2100))=0
cattle_grassland_all=sum(sum(resizecattle_2100.*worldarea))

% it is estimated that in total there will be 11,134,682 (11 million
% grassland-based cattle in 2100)

% This calculates the cattle population that is positively influenced by GN shift
cattle_positive_2100=zeros(size(resizecattle_2100))
cattle_positive_2100(turnover85_resize_flip>0)=resizecattle(turnover85_resize_flip>0)
cattle_positive_2100_count=cattle_positive_2100.*worldarea
cattle_positive_sum_2100=sum(sum(cattle_positive_2100_count))

cattle_negative_2100=zeros(size(resizecattle))
cattle_negative_2100(turnover85_resize_flip<0)=resizecattle(turnover85_resize_flip<0)
cattle_negative_2100_count=cattle_negative_2100.*worldarea
cattle_negative_sum_2100=sum(sum(cattle_negative_2100_count))


% In total there are 1.2 billion grassland-based cattle. 420 million will be
% negatively influenced, and 92 million will be positively influenced

%% produce graph: negatively impacted cattle
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(cattle_negative_2100));
geoshow(flipud(cattle_negative_2100), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'head/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,50])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/cattle_negative.svg');
%% produce graph: positively impacted cattle
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(287)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(cattle_positive_2100));
geoshow(flipud(cattle_positive_2100), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'head/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,50])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/cattle_positive.svg');
%%
%% rescale grassland influenced cattle
% Step 2: Calculate the scaling factor to preserve the total population
resized_cattle_positive= imresize(cattle_positive_2100_count,320/3600);
original_sum = sum(cattle_positive_2100_count(:));  % Sum of the original data
resized_sum = sum(resized_cattle_positive(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_cattle_positive = double(resized_cattle_positive * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_cattle_positive(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale cattle
resized_cattle_negative_2100= imresize(cattle_negative_2100_count,320/3600);
original_sum = sum(cattle_negative_2100_count(:));  % Sum of the original data
resized_sum = sum(resized_cattle_negative_2100(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_cattle_negative_2100 = double(resized_cattle_negative_2100 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_cattle_negative_2100(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);
%% Get for all countries for negative cattle
% Initialize an empty cell array to store results
cattle_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in cattle_neg
    cattle_neg{i} = sum(sum(resized_cattle_negative_2100 .* double(country_data.(fieldName))));
end
cattle_neg_values = cell2mat(cattle_neg)';
%% Get for all countries for positive cattle
% Initialize an empty cell array to store results
cattle_posi = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in cattle_neg
    cattle_posi{i} = sum(sum(resized_cattle_positive .* double(country_data.(fieldName))));
end
cattle_posi_values = cell2mat(cattle_posi)';

%% %%%%% sheep
%% start calculation

resizesheep_2100(isnan(resizesheep_2100))=0
sheep_grassland_all=sum(sum(resizesheep_2100.*worldarea))

% This calculates the population that is positively influenced by GN shift
sheep_positive_2100=zeros(size(resizesheep_2100))
sheep_positive_2100(turnover85_resize_flip>0)=resizesheep_2100(turnover85_resize_flip>0)
sheep_positive_2100_count=sheep_positive_2100.*worldarea
sheep_positive_sum_2100=sum(sum(sheep_positive_2100_count))

sheep_negative_2100=zeros(size(resizesheep_2100))
sheep_negative_2100(turnover85_resize_flip<0)=resizesheep_2100(turnover85_resize_flip<0)
sheep_negative_2100_count=sheep_negative_2100.*worldarea
sheep_negative_sum_2100=sum(sum(sheep_negative_2100_count))

% In total there are 1.2 billion grassland-based sheep. 710 million will be
% negatively influenced, and 500 million will be positively influenced


%% produce graph: negatively impacted sheep
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(sheep_negative_2100));
geoshow(flipud(sheep_negative_2100), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'head/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,50])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/sheep_negative.svg');
%% produce graph: positively impacted sheep
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(287)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(sheep_positive_2100));
geoshow(flipud(sheep_positive_2100), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'head/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,50])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/sheep_positive.svg');
%%
%% rescale grassland influenced sheep
% Step 2: Calculate the scaling factor to preserve the total population
resized_sheep_positive_2100= imresize(sheep_positive_2100_count,320/3600);
original_sum = sum(sheep_positive_2100_count(:));  % Sum of the original data
resized_sum = sum(resized_sheep_positive_2100(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_sheep_positive_2100 = double(resized_sheep_positive_2100 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_sheep_positive_2100(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale sheep
resized_sheep_negative_2100= imresize(sheep_negative_2100_count,320/3600);
original_sum = sum(sheep_negative_2100_count(:));  % Sum of the original data
resized_sum = sum(resized_sheep_negative_2100(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_sheep_negative_2100 = double(resized_sheep_negative_2100 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_sheep_negative_2100(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);
%% Get for all countries for negative sheep
% Initialize an empty cell array to store results
sheep_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in sheep_neg
    sheep_neg{i} = sum(sum(resized_sheep_negative_2100 .* double(country_data.(fieldName))));
end
sheep_neg_values = cell2mat(sheep_neg)';
%% Get for all countries for positive sheep
% Initialize an empty cell array to store results
sheep_posi = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in sheep_neg
    sheep_posi{i} = sum(sum(resized_sheep_positive_2100 .* double(country_data.(fieldName))));
end
sheep_posi_values = cell2mat(sheep_posi)';

%% %%%%% goats
%% start calculation



resizegoats_2100(isnan(resizegoats_2100))=0
goats_grassland_all=sum(sum(resizegoats_2100.*worldarea))

% This calculates the population that is positively influenced by GN shift
goats_positive_2100=zeros(size(resizegoats_2100))
goats_positive_2100(turnover85_resize_flip>0)=resizegoats_2100(turnover85_resize_flip>0)
goats_positive_2100_count=goats_positive_2100.*worldarea
goats_positive_sum_2100=sum(sum(goats_positive_2100_count))

goats_negative_2100=zeros(size(resizegoats_2100))
goats_negative_2100(turnover85_resize_flip<0)=resizegoats_2100(turnover85_resize_flip<0)
goats_negative_2100_count=goats_negative_2100.*worldarea
goats_negative_sum_2100=sum(sum(goats_negative_2100_count))

% In total there are 698 million grassland-based sheep. 500 million will be
% negatively influenced, and 170 million will be positively influenced

%% produce graph: negatively impacted goats
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(goats_negative_2100));
geoshow(flipud(goats_negative_2100), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'head/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,50])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/goats_negative.svg');
%% produce graph: positively impacted goats
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(287)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(goats_positive_2100));
geoshow(flipud(goats_positive_2100), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'head/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,50])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/goats_positive.svg');
%%
%% rescale grassland influenced goats
% Step 2: Calculate the scaling factor to preserve the total population
resized_goats_positive_2100= imresize(goats_positive_2100_count,320/3600);
original_sum = sum(goats_positive_2100_count(:));  % Sum of the original data
resized_sum = sum(resized_goats_positive_2100(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_goats_positive_2100 = double(resized_goats_positive_2100 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_goats_positive_2100(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale goats
resized_goats_negative_2100= imresize(goats_negative_2100_count,320/3600);
original_sum = sum(goats_negative_2100_count(:));  % Sum of the original data
resized_sum = sum(resized_goats_negative_2100(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_goats_negative_2100 = double(resized_goats_negative_2100 * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_goats_negative_2100(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);
%% Get for all countries for negative goats
% Initialize an empty cell array to store results
goats_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in goats_neg
    goats_neg{i} = sum(sum(resized_goats_negative_2100 .* double(country_data.(fieldName))));
end
goats_neg_values = cell2mat(goats_neg)';
%% Get for all countries for positive goats
% Initialize an empty cell array to store results
goats_posi = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in goats_neg
    goats_posi{i} = sum(sum(resized_goats_positive_2100 .* double(country_data.(fieldName))));
end
goats_posi_values = cell2mat(goats_posi)';


%%
livestock_neg_values=goats_neg_values+cattle_neg_values+sheep_neg_values;
livestock_posi_values=goats_posi_values+cattle_posi_values+sheep_posi_values;

goats_net=goats_neg_values-goats_posi_values;
cattle_net=cattle_neg_values-cattle_posi_values;
sheep_net=sheep_neg_values-sheep_posi_values;

negative_livestock=goats_negative_sum_2100+cattle_negative_sum_2100+sheep_negative_sum_2100

% 1.6 billion livestock will be impacted?
% as opposed to 160 million population?
% 这个数字要怎么算？
%% create bar chart for net values (how many livestock will be unable to relocate within a country?)
% Create a bar chart with country names as labels
figure;
bar(goats_net);
title('Net negative impact');
xlabel('Country');
ylabel('Goats');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order 
% Sort sheep_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(goats_net, 'descend');
% Sort country names based on the sorting indices
sorted_countries = countryFields(sort_index);

% Create a bar chart with sorted values and country names as labels
figure;
bar(sorted_values);
title('Population in Grassland by Country (Descending Order)');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(sorted_countries), 'XTickLabel', sorted_countries);
xlim([0,10])
% Rotate x-axis labels for better readability
xtickangle(90);


save('./grazingniche/matdata/sensitivity_poplivestock_widest.mat','negative_livestock','goats_negative_sum_2100','cattle_negative_sum_2100','sheep_negative_sum_2100','pop_neg_sum','grassland_pop_negative_sum','*posi_values','*neg_values')
%%
openvar('goats_neg_values')
openvar('cattle_neg_values')
openvar('sheep_neg_values')