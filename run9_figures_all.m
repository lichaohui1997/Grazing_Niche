clear;clear all;
addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% vector climate and grassland data
load('final_data.mat')
pr_final_vector=sum(pr_final,2)/7;
tas_final_vector=sum(tas_final,2)/7;
sfcWind_final_vector=sum(sfcWind_final,2)/7;
hurs_final_vector=sum(hurs_final,2)/7;
%% 散点图 e.g. tas-pr

subplot(2,3,1)
X=pr_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

% 使用scatter绘图
scatter(X,Y,'filled','CData',H);
xlabel('Precipitation (mm)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Precipitation and Temperature'});
defualtAxes()


subplot(2,3,2)
X=hurs_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

% 使用scatter绘图
scatter(X,Y,'filled','CData',H);
xlabel('Relative humidity (%)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Humidity and Temperature'});
defualtAxes()

subplot(2,3,3)
X=sfcWind_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

% 使用scatter绘图
scatter(X,Y,'filled','CData',H);
xlabel('Near Surface Windspeed (m/s)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Windspeed and Temperature'});
defualtAxes()

% Save the figure as a SVG and specify the size
set(gcf, 'Position',  [751,163,1092,753])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter.svg');
%% 等高线图
subplot(2,3,1)
X=pr_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Precipitation (mm)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Precipitation and Temperature'});

defualtAxes()

subplot(2,3,2)
X=hurs_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Relative Humidity (%)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Humidity and Temperature'});
defualtAxes()

subplot(2,3,3)
X=sfcWind_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Near Surface Windspeed (m/s)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Windspeed and Temperature'});
defualtAxes()

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatterfill.svg');


%% 带直方图的散点图

PntSet=[pr_final_vector,tas_final_vector]
X=pr_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
CData=interp2(XMesh,YMesh,ZMesh,X,Y);

% 主分布图
ax1=axes('Parent',gcf);hold(ax1,'on')
scatter(ax1,PntSet(:,1),PntSet(:,2),'filled','CData',CData);
ax1.Position=[0.1,0.1,0.6,0.6];

% X轴直方图
ax2=axes('Parent',gcf);hold(ax2,'on')
histogram(ax2,PntSet(:,1),'FaceColor',addcolorplus(115),...
    'EdgeColor','none','FaceAlpha',0.7)
ax2.Position=[0.1,0.75,0.6,0.15];
ax2.YColor='none';
ax2.XTickLabel='';
ax2.TickDir='out';
ax2.XLim=ax1.XLim;

% Y轴直方图
ax3=axes('Parent',gcf);hold(ax3,'on')
histogram(ax3,PntSet(:,2),'FaceColor',addcolorplus(115),...
    'EdgeColor','none','FaceAlpha',0.7,'Orientation','horizontal')
ax3.Position=[0.75,0.1,0.15,0.6];
ax3.XColor='none';
ax3.YTickLabel='';
ax3.TickDir='out';
ax3.YLim=ax1.YLim;

set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter_histogram.svg');

%% 带灰色区域的散点图
nexttile;
hold on;
scatter(pr_final_vector, resizedhyde_vector_s,'filled');
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
scatter(tas_final_vector, resizedhyde_vector_s,'filled');
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
scatter(hurs_final_vector, resizedhyde_vector_s,'filled');
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
scatter(sfcWind_final_vector, resizedhyde_vector_s,'filled');
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
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter_grey.svg');



%% aggregate climate data (matrix form) from run4

load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;

%% grassland data
load('grassland_env.mat')
%% Niche map

% cond1_pr = (aggregate_pr >= 324 & aggregate_pr <= 2548);
% cond1_tas = (aggregate_tas >= 2 & aggregate_tas <= 28.24);
% cond1_sfcWind = (aggregate_sfcWind >= 1.13 & aggregate_sfcWind <= 6.50);
% cond_hurs = (aggregate_hurs >= 34 & aggregate_hurs <= 84);
% cond1_landuse = (grassland_env>0);

cond1_pr = (aggregate_pr >= 433 & aggregate_pr <= 2368);
cond1_tas = (aggregate_tas >= 4 & aggregate_tas <= 28.5);
cond1_sfcWind = (aggregate_sfcWind >= 1.4 & aggregate_sfcWind <= 6);
cond_hurs = (aggregate_hurs >= 46 & aggregate_hurs <= 87);
cond1_landuse = (grassland_env>0);

% Combine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% Create the coupled landuse map
landuse_coupled1 = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,90],[-180,180],size(landuse_coupled1));

figure1 = figure;
%sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
%custom_colormap = [1 1 1; addcolorplus(332)];
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
%set(gcf, 'Position',  [584,449,684,537])
set(gcf, 'Position',  [361,614,899,295])


saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/niche_general.svg');
save('./grazingniche/matdata/niche.mat',"landuse_coupled1")

%% nice graph of grassland
% to run this you need to add the command file into directory
figure1 = figure;
sgtitle('Grassland map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[0,360],size(grassland_env));
geoshow(flipud(grassland_env), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])

% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/grassland.svg');

%% grassland that is outside of the niche
outniche = (grassland_env > 1) & (landuse_coupled1 == 0); % This creates a logical matrix
outniche = double(outniche); % Converts logicals to doubles: 1 for true, 0 for false
outniche=outniche.*grassland_env

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
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(outniche), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
set(gcf, 'Position',  [584,449,684,537])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/outniche.svg');

%% 山形图-绘制七组气候数据的误差-为什么降水的没有数据？
subplot(4,1,1);
Data=mat2cell(hurs_final,ones(1,1)*10602,ones(7,1))
JP=joyPlot(Data,'ColorMode','Order');
JP=JP.draw();

subplot(4,1,2);
Data=mat2cell(pr_final,ones(1,1)*10602,ones(7,1))
JP=joyPlot(Data,'ColorMode','Order');
JP=JP.draw();

subplot(4,1,3);
Data=mat2cell(tas_final,ones(1,1)*10602,ones(7,1))
JP=joyPlot(Data,'ColorMode','Order');
JP=JP.draw();

subplot(4,1,4);
Data=mat2cell(sfcWind_final,ones(1,1)*10602,ones(7,1))
JP=joyPlot(Data,'ColorMode','Order');
JP=JP.draw();
%% from script run5
load('livestockdensity.mat')
% this AGB.mat data comes from the below section on grassland biomass.
load('AGB.mat')
load('aggregate_tas.mat')
load('niche.mat')% variable name: landuse_coupled1 (in run4)
load('livestockpercentage.mat')
load('overgraze_niche.mat')%variable name:overgraze_insideniche_logic;overgraze_outsideniche_logic
load('livestock2015.mat');%variable name: 'livestockpercentage','outniche','livestock','resizecattle','allresizecattle', 'resizesheep','allresizesheep','resizegoats','allresizegoats','resizedpop');
% carbon loss due to grazing
load('AGBCBGBC.mat');
load('SOC.mat');
load("worldarea.mat");
load('treecover.mat');
load('carbonloss_map');


%% Figure 3: subfigure 1: livestock density, grassland livestock population distrubution
figure1 = figure;
%sgtitle('Grassland livestock distribution', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(livestock));
geoshow(flipud(livestock), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'heads (livestock units)/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])

set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/livestock_density.svg');
%% cattle density
figure1 = figure;
%sgtitle('Grassland livestock distribution', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
%custom_colormap = [1 1 1; addcolorplus(306)];
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(resizecattle));
geoshow(flipud(resizecattle), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'heads (livestock units)/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])

set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/cattle_density.svg');
%% cattle+sheep density
figure1 = figure;
%sgtitle('Grassland livestock distribution', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(306)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(resizecattle));
geoshow(flipud(resizecattle+0.2.*resizesheep), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'heads (livestock units)/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])

set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/cattlesheep_density.svg');

%% Figure 3, subfigure 2: empty niche
load('empty_niche.mat');
figure15 = figure;
%sgtitle('Empty niche (Niche under grazing capacity)', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(282)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(empty_niche), R, 'DisplayType', 'texturemap');   
load coastlines
hc=colorbar;
title(hc, 'km^2/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/emptyniche.svg');

%% livestock capacity
figure1 = figure;
sgtitle('percentage of consumed biomass', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
% c = colorbar;  
% c.Ruler.TickLabelFormat='%g%%';
colorbar
R = georefcells([-90,90],[-180,180],size(grassland_env));
geoshow(flipud(livestockpercentage), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
set(gcf, 'Position',  [584,449,684,537])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/livestockcapacity.svg');
% 从这张图可以看出来两级分化现象非常严重，不是没有使用就是使用过度。
%% graph of outside niche+overgraze and inside niche+overgraze

%% overgrazed grassland inside niche
% figure7 = figure;
% sgtitle('Overgrazed grassland inside of the niche', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(276)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(overgraze_insideniche_logic), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% hc=colorbar;
% title(hc,'km^2/cell');
% 
% set(gcf, 'Position',  [456,474,465,437])
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/insidenicheover_graze.svg');
%% Figure 3, subfigure 3: Overgrazed non-GN grassland: they shouldn't be grazed!
figure8 = figure;
%sgtitle('Overgrazed grassland outside of the niche', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(333)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(overgraze_outsideniche_logic.*grassland_env), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
hc=colorbar;
title(hc, 'km^2/cell', 'FontSize', 14)

caxis([0,60])

set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/outsideniche_overgraze.svg');
%% Figure 3, subfigure 4: carbon loss 

load('carbonloss_map.mat');
carbonloss_map=carbonloss_map./100; %this is to convert the km^2 back to ha
carbonloss_map(carbonloss_map>5000)=5000;
figure15 = figure;
%sgtitle('Map of carbon loss due to overgrazing', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(288)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc, 'MgC/ha', 'FontSize', 14)

R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(carbonloss_map), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/carbonloss.svg');

%% Figure 3, subfigure 6:bar figure of carbon loss
figure_carbonloss={'SOC loss',carbonloss_SOC;
    'AGBC loss',carbonloss_AGBC;
    'BGBC loss',carbonloss_BGBC}
figure_carbonloss_value=[carbonloss_SOC;carbonloss_AGBC;carbonloss_BGBC]
figure_carbonloss_caption={'SOC loss';'AGBC loss';'BGBC loss'}
bar(figure_carbonloss_value)
xticklabels(figure_carbonloss_caption)
ylabel('Carbon loss due to grazing (MgC)');
set(gcf, 'Position',  [584,545,418,456])
ax = gca
ax.FontSize = 14;
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/carbonlossbar.svg');

%% bar chart of overgraze in and out of niche
%% load these first from run7.m analysis
% 这里提前跑了run7.m的数据，为了和前面的figure3接上。
load('figure_nichelandarea_value');
load('figure_overgrazeGN_value.mat');
%
figure_overgrazedGN_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}
% Create the bar chart
barh(figure_overgrazeGN_value);
% Set x-axis labels
set(gca, 'yticklabel', figure_overgrazedGN_caption);
% Add a title and labels
title('Overgrazed Area in GN vs. Outside GN');
ylabel('Region');
xlabel('Overgrazed Area (km^2)');
legend('Overgrazed Area in GN', 'Overgrazed Area outside GN');%% how much of the unusable grassland are overgrazed?
ax = gca
ax.FontSize = 14;
%get(gcf, 'Position')
set(gcf, 'Position',  [312,381,543,578])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/overgrazeGNbar.svg');

%% carbon loss due to grazing
%% grassland above ground biomass
figure9 = figure;
sgtitle('Grassland Aboveground Biomass', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(333)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'gC/m^2/yr');
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(AGB), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/AGB.svg');

%% grassland above ground biomass Carbon
figure10 = figure;
sgtitle('AGBC', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(302)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'MgC/ha');
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(AGBC), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

% 这个数值看起来可能会非常的小。因为是草地的。大多数AGBC都是在树里面的
set(gcf, 'Position',  [584,449,684,537]);
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/AGBC.svg');
%% BGBC
figure11 = figure;
sgtitle('BGBC', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(300)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'MgC/ha');
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(BGBC), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/BGBC.svg');
%% SOC

SOC0_200(SOC0_200>1000)=1000;

figure12 = figure;
sgtitle('SOC (0-2m)', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(300)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'MgC/ha');
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(SOC0_200), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/SOC.svg');
%% tree cover map
treecover(treecover>100) = 0;

figure13 = figure;
sgtitle('Tree cover map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(291)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(treecover), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/treecover.svg');
%% NPP
npp(npp>249)=0;

figure14 = figure;
sgtitle('Net primary productivity', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(284)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'gC/m^2/yr');
R = georefcells([-90,90],[-180,180],size(outniche));
geoshow(flipud(npp), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

set(gcf, 'Position',  [584,449,684,537])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/NPP.svg');
%% run6 needs to be run entirely. Too many variables. Script has simple logic, produces only 1 important graph. Does not worth summarizing.
%% from run7.m analysis
load('figure_nichelandarea_value');
load('figure_overgrazeGN_value.mat');
%% bar figure of niche area 

figure_nichelandarea_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}

bar(figure_nichelandarea_value);
% Set x-axis labels
set(gca, 'xticklabel', figure_nichelandarea_caption);
% Add a title and labels
title('GN grassland vs Non-GN grassland');
xlabel('Region');
ylabel('Grassland area (km^2)');
legend('GN area', 'Non-GN area');%% how much of the unusable grassland are overgrazed?
ax = gca
ax.FontSize = 14;
set(gcf, 'Position',  [324,437,918,428])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/nichelandareabar.svg');

%% from run8
load('figure_rcp_value.mat')
%% bar chart of future niche area
figure_rcp_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}
bar(figure_rcp_value)
% Set x-axis labels
set(gca, 'xticklabel', figure_rcp_caption);
% Add a title and labels
title('Future GN area of continents');
xlabel('Region');
ylabel('Grassland area (km^2)');
ax = gca
ax.FontSize = 14;
legend('present','rcp2.6', 'rcp4.5','rcp6.0','rcp8.5');
set(gcf, 'Position',  [218,293,1345,433])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_bar.svg');
%% run run6+run6_regional. The figures there are not included in this script.