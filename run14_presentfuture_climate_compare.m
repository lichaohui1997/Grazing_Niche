 %% This script compares the differences between the temperature between 2100 and present
% first run script 6. Future data such as data_tas_rcp85 comes from script 6
imagesc(data_tas_rcp85)
imagesc(data_hurs_rcp85)
imagesc(data_pr_rcp85)
 
colormap("turbo")
%colormap("hsv")
caxis([0,2000])
colorbar
%%
load('./grazingniche/matdata/aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('./grazingniche/matdata/aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('./grazingniche/matdata/aggregate_tas.mat')
aggregate_tas=aggregate_data;
load('./grazingniche/matdata/aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
% land mask
landmask=aggregate_sfcWind==66|aggregate_sfcWind==0
landmask=imresize(landmask,[160,320])%在逻辑值的时候就要resize，不要等到double之后再resize
landmask=double(landmask)
landmask(landmask == 1) = NaN;
landmask(landmask == 0) = 1;
imagesc(landmask)
colorbar
% sfcdata
aggregate_sfcWind(aggregate_sfcWind>20)=0
aggregate_sfcWind(aggregate_sfcWind>10)=0

load('grassland_env.mat')

%% check the data really quickly
imagesc(aggregate_hurs)
colormap('turbo')
colorbar

imagesc(aggregate_pr)
colormap('turbo')
colorbar
caxis([0,2000])
%%
data_tas_present=imresize(aggregate_tas,[160,320])
data_hurs_present=imresize(aggregate_hurs,[160,320])
data_sfcWind_present=imresize(aggregate_sfcWind,[160,320])
data_pr_present=imresize(aggregate_pr,[160,320])


data_tas_rcp85 = data_tas_rcp85(:, [ceil(end/2+1):end, 1:floor(end/2)]);
data_pr_rcp85 = data_pr_rcp85(:, [ceil(end/2+1):end, 1:floor(end/2)]);
data_hurs_rcp85 = data_hurs_rcp85(:, [ceil(end/2+1):end, 1:floor(end/2)]);
data_sfcWind_rcp85 = data_sfcWind_rcp85(:, [ceil(end/2+1):end, 1:floor(end/2)]);

%%

imagesc(data_tas_present)
colorbar

imagesc(data_tas_rcp85)
colorbar

data_tas_diff=data_tas_rcp85-data_tas_present
data_pr_diff=data_pr_rcp85-data_pr_present
data_hurs_diff=data_hurs_rcp85-data_hurs_present
data_sfcWind_diff=data_sfcWind_rcp85-data_sfcWind_present


%%
% This is actually a very useful color pallet. Can adjst for different colors
num_colors = 256;  % Total number of colors for smooth transition

% Define specific color points: dark red, white, dark green
color_points = [ 
    0.6, 0.0, 0.0;  % Dark Red
    1.0, 1.0, 1.0;  % White
    0.1, 0.8, 0.3;  % Dark Green
];

% Interpolate to create a smooth transition colormap
reds_to_white = [linspace(color_points(1,1), color_points(2,1), num_colors/2)', ...
                 linspace(color_points(1,2), color_points(2,2), num_colors/2)', ...
                 linspace(color_points(1,3), color_points(2,3), num_colors/2)'];

white_to_greens = [linspace(color_points(2,1), color_points(3,1), num_colors/2)', ...
                   linspace(color_points(2,2), color_points(3,2), num_colors/2)', ...
                   linspace(color_points(2,3), color_points(3,3), num_colors/2)'];

% Combine the two gradients into a single colormap
custom_colormap = [reds_to_white; white_to_greens];

%% tas diff
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
%custom_colormap = [reds_to_white; white_to_greens];
nexttile
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_tas_diff));
geoshow(flipud(data_tas_diff.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
caxis([-6,10])
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/diff_tas.svg');
%% pr diff
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
%custom_colormap = [reds_to_white; white_to_greens];

ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0fmm';
R = georefcells([-90,90],[-180,180],size(data_pr_diff));
geoshow(flipud(data_pr_diff.*landmask), R, 'DisplayType', 'texturemap');   
%caxis([-20,30])
load coastlines
%caxis([-6,10])
caxis([-400,400])
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/diff_pr.svg');

%% hurs diff
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%g%%';
R = georefcells([-90,90],[-180,180],size(data_hurs_diff));
geoshow(flipud(data_hurs_diff.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
caxis([-20,40])
plotm(coastlat, coastlon, 'Color', 'black');
set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/diff_hurs.svg');

%% sfcWind diff
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0fm/s';
R = georefcells([-90,90],[-180,180],size(data_sfcWind_diff));
geoshow(flipud(data_sfcWind_diff.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
caxis([-5,5])
plotm(coastlat, coastlon, 'Color', 'black');
set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/diff_sfcWind.svg');

%%
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
%custom_colormap = [reds_to_white; white_to_greens];

ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_hurs_present));
geoshow(flipud(data_hurs_present.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
%%
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
%custom_colormap = [reds_to_white; white_to_greens];

ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_pr_present));
geoshow(flipud(data_pr_present.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
caxis([0,3000])
plotm(coastlat, coastlon, 'Color', 'black'); 
%%
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
%custom_colormap = [reds_to_white; white_to_greens];

ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_sfcWind_present));
geoshow(flipud(data_sfcWind_present.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
caxis([0,10])
plotm(coastlat, coastlon, 'Color', 'black'); 


%%
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_hurs_diff));
geoshow(flipud(data_hurs_rcp85.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
%caxis([-6,10])
plotm(coastlat, coastlon, 'Color', 'black'); 
%%
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_hurs_diff));
geoshow(flipud(data_tas_rcp85.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
%caxis([-6,10])
plotm(coastlat, coastlon, 'Color', 'black'); 
%%
figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(310)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat = '%.0f°C';
R = georefcells([-90,90],[-180,180],size(data_hurs_diff));
geoshow(flipud(data_pr_rcp85.*landmask), R, 'DisplayType', 'texturemap');   
load coastlines
caxis([0,3000])
plotm(coastlat, coastlon, 'Color', 'black'); 
%%
nexttile
imagesc(data_tas_present)
colorbar

nexttile
imagesc(data_hurs_present)
colorbar

nexttile
imagesc(data_pr_present)
colorbar

nexttile
imagesc(data_sfcWind_present)
colorbar
%%
nexttile
imagesc(data_tas_rcp85)
colorbar

nexttile
imagesc(data_hurs_rcp85)
colorbar

nexttile
imagesc(data_pr_rcp85)
colorbar

nexttile
imagesc(data_sfcWind_rcp85)
colorbar