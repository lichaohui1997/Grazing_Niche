%% This script calculates the turnover impact on population and livestock figure
%% impact on population and livestock
% step 1: data preparation
load("grassland_env.mat")
load("livestockdensity.mat")
load('turnover_average.mat')% this mat file includes these variables:compare_turnover_rcp26-compare_turnover_rcp85
% population count (population per pixel)
[pop_count, Rpop] = readgeoraster('gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_15_min.tif')
pop_count(pop_count<0)=0
%[resizedpop_count,resizedRpop] = georesize(pop_count,Rpop,1800/720,"bilinear");
% population per square kilometer
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
pop(pop<0)=0
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");

%% rescale population
% Step 2: Calculate the scaling factor to preserve the total population
Rpop_count = georefcells([-90,90],[-180,180],size(pop_count));
[resizedpop_count,resizedRpop] = georesize(pop_count,Rpop_count,1800/720,"bilinear");
original_sum = sum(pop_count(:));  % Sum of the original data
resized_sum = sum(resizedpop_count(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resizedpop_count = double(resizedpop_count * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resizedpop_count(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

sum(sum(resizedpop_count))
sum(resizedpop_count(:))
%% flip variables and check if they are aligned
turnover85_resize_average=imresize(compare_turnover_rcp85_average,[1800,3600])
imagesc(turnover85_resize_average)
turnover85_resize_flip = double(turnover85_resize_average(:, [ceil(end/2+1):end, 1:floor(end/2)]));
grassland_env_flip = grassland_env(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% step 2: check if all variables are aligned
imagesc(turnover85_resize_flip)
imagesc(grassland_env_flip)
imagesc(resizedpop_count)
imagesc(resizecattle)
colorbar

%% start calculation
% This calculates the population that is inhabited near grasslands
grassland_pop=zeros(size(resizedpop_count))
grassland_pop(livestock>0&resizedpop<20)=resizedpop_count(livestock>0&resizedpop<20)
grassland_pop_sum=sum(sum(grassland_pop))

imagesc(grassland_pop)
colorbar

% This calculates the population that is positively influenced by GN shift
grassland_pop_positive=zeros(size(resizedpop_count))
grassland_pop_positive(turnover85_resize_flip>0)=grassland_pop(turnover85_resize_flip>0)
grassland_pop_positive_sum=sum(sum(grassland_pop_positive))

grassland_pop_negative=zeros(size(resizedpop_count))
grassland_pop_negative(turnover85_resize_flip<0)=grassland_pop(turnover85_resize_flip<0)
grassland_pop_negative_sum=sum(sum(grassland_pop_negative))

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
caxis([0,500])
set(gcf, 'Position',  [361,651,489,258])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_negative.svg');
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
caxis([0,500])
set(gcf, 'Position',  [361,651,489,258])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_positive.svg');
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
%% %%%%%%%%% cattle

%% start calculation
% This calculates the population that is inhabited near grasslands
cattle_pop=zeros(size(resizecattle))
%cattle(grassland_env>0)=resizedpop_count(grassland_env>0)
cattle(livestock>0&resizedpop>20)=resizedpop_count(livestock>0&resizedpop>20)
cattle_sum=sum(sum(resizecattle))

% This calculates the population that is positively influenced by GN shift
cattle_positive=zeros(size(resizecattle))
cattle_positive(turnover85_resize_flip>0)=resizecattle(turnover85_resize_flip>0)
cattle_positive_sum=sum(sum(cattle_positive))

cattle_negative=zeros(size(resizecattle))
cattle_negative(turnover85_resize_flip<0)=resizecattle(turnover85_resize_flip<0)
cattle_negative_sum=sum(sum(cattle_negative))

% In total there are 5 million grassland-based cattle. 3 million will be
% negatively influenced, and 2 million will be positively influenced
% negatively influenced, 50 million will be positively
% influenced

%% produce graph: negatively impacted cattle
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(cattle_negative));
geoshow(flipud(cattle_negative), R, 'DisplayType', 'texturemap');   
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
R = georefcells([-90,90],[-180,180],size(cattle_positive));
geoshow(flipud(cattle_positive), R, 'DisplayType', 'texturemap');   
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
resized_cattle_positive= imresize(cattle_positive,320/3600);
original_sum = sum(cattle_positive(:));  % Sum of the original data
resized_sum = sum(resized_cattle_positive(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_cattle_positive = double(resized_cattle_positive * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_cattle_positive(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale cattle
resized_cattle_negative= imresize(cattle_negative,320/3600);
original_sum = sum(cattle_negative(:));  % Sum of the original data
resized_sum = sum(resized_cattle_negative(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_cattle_negative = double(resized_cattle_negative * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_cattle_negative(:));  % Should be close to original_sum
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
    cattle_neg{i} = sum(sum(resized_cattle_negative .* double(country_data.(fieldName))));
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
% This calculates the population that is inhabited near grasslands
sheep_pop=zeros(size(resizesheep))
%sheep(grassland_env>0)=resizedpop_count(grassland_env>0)
sheep(livestock>0&resizedpop>20)=resizedpop_count(livestock>0&resizedpop>20)
sheep_sum=sum(sum(resizesheep))

% This calculates the population that is positively influenced by GN shift
sheep_positive=zeros(size(resizesheep))
sheep_positive(turnover85_resize_flip>0)=resizesheep(turnover85_resize_flip>0)
sheep_positive_sum=sum(sum(sheep_positive))

sheep_negative=zeros(size(resizesheep))
sheep_negative(turnover85_resize_flip<0)=resizesheep(turnover85_resize_flip<0)
sheep_negative_sum=sum(sum(sheep_negative))

% In total there are 5 million grassland-based sheep. 3 million will be
% negatively influenced, and 2 million will be positively influenced
% negatively influenced, 50 million will be positively
% influenced

%% produce graph: negatively impacted sheep
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(sheep_negative));
geoshow(flipud(sheep_negative), R, 'DisplayType', 'texturemap');   
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
R = georefcells([-90,90],[-180,180],size(sheep_positive));
geoshow(flipud(sheep_positive), R, 'DisplayType', 'texturemap');   
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
resized_sheep_positive= imresize(sheep_positive,320/3600);
original_sum = sum(sheep_positive(:));  % Sum of the original data
resized_sum = sum(resized_sheep_positive(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_sheep_positive = double(resized_sheep_positive * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_sheep_positive(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale sheep
resized_sheep_negative= imresize(sheep_negative,320/3600);
original_sum = sum(sheep_negative(:));  % Sum of the original data
resized_sum = sum(resized_sheep_negative(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_sheep_negative = double(resized_sheep_negative * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_sheep_negative(:));  % Should be close to original_sum
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
    sheep_neg{i} = sum(sum(resized_sheep_negative .* double(country_data.(fieldName))));
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
    sheep_posi{i} = sum(sum(resized_sheep_positive .* double(country_data.(fieldName))));
end
sheep_posi_values = cell2mat(sheep_posi)';

%% %%%%% goats
%% start calculation
% This calculates the population that is inhabited near grasslands
goats_pop=zeros(size(resizegoats))
%goats(grassland_env>0)=resizedpop_count(grassland_env>0)
goats(livestock>0&resizedpop>20)=resizedpop_count(livestock>0&resizedpop>20)
goats_sum=sum(sum(resizegoats))

% This calculates the population that is positively influenced by GN shift
goats_positive=zeros(size(resizegoats))
goats_positive(turnover85_resize_flip>0)=resizegoats(turnover85_resize_flip>0)
goats_positive_sum=sum(sum(goats_positive))

goats_negative=zeros(size(resizegoats))
goats_negative(turnover85_resize_flip<0)=resizegoats(turnover85_resize_flip<0)
goats_negative_sum=sum(sum(goats_negative))

% In total there are 5 million grassland-based goats. 3 million will be
% negatively influenced, and 2 million will be positively influenced
% negatively influenced, 50 million will be positively
% influenced

%% produce graph: negatively impacted goats
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(275)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(goats_negative));
geoshow(flipud(goats_negative), R, 'DisplayType', 'texturemap');   
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
R = georefcells([-90,90],[-180,180],size(goats_positive));
geoshow(flipud(goats_positive), R, 'DisplayType', 'texturemap');   
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
resized_goats_positive= imresize(goats_positive,320/3600);
original_sum = sum(goats_positive(:));  % Sum of the original data
resized_sum = sum(resized_goats_positive(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_goats_positive = double(resized_goats_positive * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_goats_positive(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% negative rescale goats
resized_goats_negative= imresize(goats_negative,320/3600);
original_sum = sum(goats_negative(:));  % Sum of the original data
resized_sum = sum(resized_goats_negative(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resized_goats_negative = double(resized_goats_negative * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resized_goats_negative(:));  % Should be close to original_sum
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
    goats_neg{i} = sum(sum(resized_goats_negative .* double(country_data.(fieldName))));
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
    goats_posi{i} = sum(sum(resized_goats_positive .* double(country_data.(fieldName))));
end
goats_posi_values = cell2mat(goats_posi)';


%%
livestock_neg_values=goats_neg_values+cattle_neg_values+sheep_neg_values;
livestock_posi_values=goats_posi_values+cattle_posi_values+sheep_posi_values;

goats_net=goats_neg_values-goats_posi_values;
cattle_net=cattle_neg_values-cattle_posi_values;
sheep_net=sheep_neg_values-sheep_posi_values;

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


