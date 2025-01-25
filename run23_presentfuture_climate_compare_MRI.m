%% This script compares the differences between the temperature between 2100 and present using the same set of MRI data
%% Get land mask
load('./grazingniche/matdata/aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
% land mask
landmask=aggregate_sfcWind==66|aggregate_sfcWind==0
landmask=imresize(landmask,[160,320])
landmask=double(landmask)
landmask(landmask == 1) = NaN;
landmask(landmask == 0) = 1;
%% run script 6. and script 20. Future data such as data_tas_rcp85 comes from script 6
data_tas_diff=data_tas_rcp85-data_tas_rcp85_present
data_pr_diff=data_pr_rcp85-data_pr_rcp85_present
data_hurs_diff=data_hurs_rcp85-data_hurs_rcp85_present
data_sfcWind_diff=data_sfcWind_rcp85-data_sfcWind_rcp85_present
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
caxis([-10,10])
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
caxis([-2,2])
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