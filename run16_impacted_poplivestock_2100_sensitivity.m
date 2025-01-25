%% This script calculates the uncertainty for the impacted population using the widest and thinnest thresholds
%% Need to alternate (replace thinnest/widest) to get the values

pop_data_2100=ncread("/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/populationdensity/NetCDF/SSP3_NetCDF/total/NetCDF/ssp3_2100.nc","ssp3_2100")
openvar("pop_data_2100")

pop_data_2100_info=ncinfo("/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/populationdensity/NetCDF/SSP3_NetCDF/total/NetCDF/ssp3_2100.nc")
openvar("pop_data_2100_info")

pop_data_2100_1=pop_data_2100'
imagesc(pop_data_2100_1)
NorthZeros=zeros(53,2880)
SouthZeros=zeros(270,2880)
pop_data_2100_2=[NorthZeros;pop_data_2100_1;SouthZeros]

imagesc(pop_data_2100_2)


figure9 = figure;
sgtitle('Grassland outside of the niche', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[-180,180],size(pop_data_2100_2));
geoshow(flipud(pop_data_2100_2), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 


pop_data_2100_2(pop_data_2100_2==2147483647)=0
sum(sum(pop_data_2100_2))
pop_count_2100=pop_data_2100_2
save('./grazingniche/matdata/pop_count_2100.mat','pop_count_2100')
%% This script calculates the turnover impact on population and livestock figure
%% impact on population and livestock
% step 1: data preparation
load("grassland_env.mat")
load("livestockdensity.mat")
load('turnover_thinnest.mat')% this mat file includes these variables:compare_turnover_rcp26-compare_turnover_rcp85
load("pop_count_2100.mat")
pop_count=pop_count_2100;

% population count (population per pixel)
% [pop_count, Rpop] = readgeoraster('gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_15_min.tif')
% pop_count(pop_count<0)=0
%[resizedpop_count,resizedRpop] = georesize(pop_count,Rpop,1800/720,"bilinear");
% population per square kilometer
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
pop(pop<0)=0
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");

%% rescale population
% Step 2: Calculate the scaling factor to preserve the total population
Rpop_count = georefcells([-90,90],[-180,180],size(pop_count));
[resizedpop_count, resizedRpop] = georesize(pop_count, Rpop_count, 3600/2880, "bilinear");
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
turnover85_resize_thinnest=imresize(compare_turnover_rcp85_thinnest,[1800,3600])
imagesc(turnover85_resize_thinnest)
turnover85_resize_flip = double(turnover85_resize_thinnest(:, [ceil(end/2+1):end, 1:floor(end/2)]));

% step 2: check if all variables are aligned
nexttile
imagesc(turnover85_resize_flip)
nexttile
imagesc(resizedpop_count)
nexttile
imagesc(resizecattle)
colorbar

figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[-180,180],size(resizedpop_count));
geoshow(flipud(resizedpop_count), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 


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

