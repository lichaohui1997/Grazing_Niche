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
%% 等高线图 niche
subplot(2,3,1)
X=pr_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,10000,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Precipitation (mm)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Historical Distribution'; 'Precipitation and Temperature'});
defualtAxes()

subplot(2,3,2)
X=hurs_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,100,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
%F = ksdensity([X, Y], [XMesh(:), YMesh(:)], 'Bandwidth', [10, 3]); 
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Relative Humidity (%)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Historical Distribution'; 'Humidity and Temperature'});
defualtAxes()

subplot(2,3,3)
X=sfcWind_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,10,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Near Surface Windspeed (m/s)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Historical Distribution'; 'Windspeed and Temperature'});
defualtAxes()

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatterfill.svg');

%%
colorbar_handle = colorbar; % Add a colorbar
colorbar_handle.Ticks = [0 1]; % Set ticks at the bottom and top
colorbar_handle.TickLabels = {'Low', 'High'}; % Add "Low" and "High" labels
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_colorbar.svg');

%% Max suggestion: grassland/distribution (Lenton paper Fig1G)
%% final step
% my pr_5 is 12960*4 variable. resizedhyde_vector is 12960*1 variable. I want the rows values in the resizedhyde_vector is smaller than 0, correspondingly, the rows inpr_5 will be deleted
load('grassland_env.mat')
grassland_env_resized = imresize(grassland_env, [90 144]);
imagesc(grassland_env_resized)
colorbar

grassland_env_vector=reshape(grassland_env_resized,[],1)

pr_max_final(grassland_env_vector<=0,:)=[]
tas_max_final(grassland_env_vector<=0,:)=[]
sfcWind_max_final(grassland_env_vector<=0,:)=[]
hurs_max_final(grassland_env_vector<=0,:)=[]

%% scatter图 max
subplot(2,3,1)
X=pr_max_final;
Y=tas_max_final;
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
title({'Grassland Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
defualtAxes()


subplot(2,3,2)
X=hurs_max_final;
Y=tas_max_final;
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
title({'Grassland Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])


subplot(2,3,3)
X=sfcWind_max_final;
Y=tas_max_final;
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
title({'Grassland Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])


% Save the figure as a SVG and specify the size
set(gcf, 'Position',  [751,163,1092,753])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter_max.svg');
%% max：grassland等高线图
subplot(2,3,1)
X=pr_max_final;
Y=tas_max_final;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,10000,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F = ksdensity([X, Y], [XMesh(:), YMesh(:)], 'Bandwidth', [800, 10]); 
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Precipitation (mm)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Grassland Distribution'; 'Precipitation and Temperature'});

defualtAxes()

subplot(2,3,2)
X=hurs_max_final;
Y=tas_max_final;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,100,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F = ksdensity([X, Y], [XMesh(:), YMesh(:)], 'Bandwidth', [30, 5]); 
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Relative Humidity (%)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Grassland Distribution'; 'Humidity and Temperature'});
defualtAxes()

subplot(2,3,3)
X=sfcWind_max_final;
Y=tas_max_final;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,10,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Near Surface Windspeed (m/s)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Grassland Distribution'; 'Windspeed and Temperature'});
defualtAxes()

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayheat_max.svg');

%% scatter figure: present livestock distribution
%% data preparation
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/livestockdensity.mat')
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_tas.mat')
aggregate_tas=aggregate_data;

aggregate_sfcWind(aggregate_sfcWind==66)=0
aggregate_sfcWind(aggregate_sfcWind>20)=0

imagesc(aggregate_sfcWind)
colorbar


imagesc(sfcWind_resize)
colorbar


grassland_env_resized = imresize(grassland_env, [90 144]);

livestock_resize=imresize(livestock, [90 144]);
pr_resize=imresize(aggregate_pr, [90 144]);
tas_resize=imresize(aggregate_tas, [90 144]);
hurs_resize=imresize(aggregate_hurs, [90 144]);
sfcWind_resize=imresize(aggregate_sfcWind, [90 144]);

livestock_vector=reshape(livestock_resize,[],1)
pr_livestock_vector=pr_resize(:)
tas_livestock_vector=tas_resize(:)
hurs_livestock_vector=hurs_resize(:)
sfcWind_livestock_vector=sfcWind_resize(:)

pr_livestock_vector(livestock_vector<=2,:)=[]
tas_livestock_vector(livestock_vector<=2,:)=[]
sfcWind_livestock_vector(livestock_vector<=2,:)=[]
hurs_livestock_vector(livestock_vector<=2,:)=[]
%% start scatter map for all livestock
subplot(2,3,1)
X=pr_livestock_vector;
Y=tas_livestock_vector;
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
title({'Livestock Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
xlim([0,10000])
defualtAxes()


subplot(2,3,2)
X=hurs_livestock_vector;
Y=tas_livestock_vector;
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
title({'Livestock Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,100])

subplot(2,3,3)
X=sfcWind_livestock_vector;
Y=tas_livestock_vector;
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
title({'Livestock Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,10])


% Save the figure as a SVG and specify the size
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter_livestock.svg');

%% livestock等高线图
subplot(2,3,1)
X=pr_livestock_vector;
Y=tas_livestock_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,10000,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Precipitation (mm)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Livestock Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
xlim([0,10000])
defualtAxes()


subplot(2,3,2)
X=hurs_livestock_vector;
Y=tas_livestock_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,100,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);

ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Relative Humidity (%)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Livestock Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,100])

subplot(2,3,3)
X=sfcWind_livestock_vector;
Y=tas_livestock_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(0,10,n);
YList=linspace(-40,40,n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);
% 绘制等高线填充图
hold on
contourf(XMesh,YMesh,ZMesh,15,'EdgeColor','none')
xlabel('Near Surface Windspeed (m/s)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Livestock Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,10])

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayheat_livestock.svg');

%% future grassland in LUH dataset
load('futureland.mat')
imagesc(futureland_rcp85)

figure9 = figure;
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
R = georefcells([-90,90],[-180,180],size(futureland_rcp85));
geoshow(flipud(futureland_rcp85), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
set(gcf, 'Position',  [584,449,684,537])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/future_grassland_85.svg');
%% turnover colormap
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
%% Africa
figure1=figure
custom_colormap = [1 1 1; addcolorplus(309)];
ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
setm(ax, 'FFaceColor', [1 1 1]);    
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90,90],[-180,180],size(livestock));
geoshow(flipud(livestock), R, 'DisplayType', 'texturemap');   
load coastlines
setm(ax, 'MapLatLimit', [-35, 40], 'MapLonLimit', [-20, 55]); % Focus on Africa
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
%set(gcf, 'Position',  [361,614,899,295])
set(figure1, 'Position', [100, 100, 900, 500]); % Wider figure for the map
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Africa_modern.svg');

figure2=figure
ydata = linspace(90, -90, 1800);
latitude_Africa_norm=latitude_Africa./sum(latitude_Africa)
plot(latitude_Africa_norm,ydata)
ylim([-35,40])
%set(gcf, 'Position',  [361,614,189,434])
set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Africa_lat.svg');

%% South America
% First figure: Map of South America
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(309)];
ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
setm(ax, 'FFaceColor', [1 1 1]); % Set background color
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90, 90], [-180, 180], size(livestock)); % Full global extent
geoshow(flipud(livestock), R, 'DisplayType', 'texturemap'); % Display livestock data
setm(ax, 'MapLatLimit', [-60, 15], 'MapLonLimit', [-90, -30]); % Adjusted for South America
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black');
caxis([0, 150]);
set(figure1, 'Position', [100, 100, 900, 500]); % Wider figure for the map

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/SouthAmerica_modern.svg');

% Second figure: Latitude Profile Plot for South America
figure2 = figure;
ydata = linspace(90, -90, 1800); % Adjust latitude range if specific to South America
latitude_SouthAmerica_norm = latitude_SouthAmerica ./ sum(latitude_SouthAmerica); % Assuming `latitude_SouthAmerica` data
plot(latitude_SouthAmerica_norm, ydata, 'LineWidth', 1.5);
ylim([-60, 15]); % Adjusted for South America's latitude range
xlabel('Normalized Value');
ylabel('Latitude');
title('Latitude Profile for South America');
set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

% Save the figure as an SVG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/SouthAmerica_lat.svg');
%% Oceania
% First figure: Map of Oceania
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(309)];
ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
setm(ax, 'FFaceColor', [1 1 1]); % Set background color
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90, 90], [-180, 180], size(livestock)); % Full global extent
geoshow(flipud(livestock), R, 'DisplayType', 'texturemap'); % Display livestock data
setm(ax, 'MapLatLimit', [-50, 10], 'MapLonLimit', [110, 180]); % Adjusted for Oceania
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black');
caxis([0, 150]);
set(figure1, 'Position', [100, 100, 900, 500]); % Wider figure for the map

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Oceania_modern.svg');

% Second figure: Latitude Profile Plot for Oceania
figure2 = figure;
ydata = linspace(90, -90, 1800); % Adjust latitude range if specific to Oceania
latitude_Oceania_norm = latitude_Oceania ./ sum(latitude_Oceania); % Assuming `latitude_Oceania` data
plot(latitude_Oceania_norm, ydata, 'LineWidth', 1.5);
ylim([-50, 10])
xlabel('Normalized Value');
ylabel('Latitude');
title('Latitude Profile for Oceania');
set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

% Save the figure as an SVG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Oceania_lat.svg');
%% Asia
% First figure: Map of Asia
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(309)];
ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
setm(ax, 'FFaceColor', [1 1 1]); % Set background color
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90, 90], [-180, 180], size(livestock)); % Full global extent
geoshow(flipud(livestock), R, 'DisplayType', 'texturemap'); % Display livestock data
setm(ax, 'MapLatLimit', [5, 55], 'MapLonLimit', [60, 150]); % Adjusted for Asia
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black');
caxis([0, 150]);
set(figure1, 'Position', [100, 100, 900, 500]); % Wider figure for the map

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Asia_modern.svg');

% Second figure: Latitude Profile Plot for Asia
figure2 = figure;
ydata = linspace(90, -90, 1800); % Adjust latitude range if specific to Asia
latitude_Asia_norm = latitude_Asia ./ sum(latitude_Asia); % Assuming `latitude_Asia` data
plot(latitude_Asia_norm, ydata, 'LineWidth', 1.5);
ylim([5, 55]);
xlabel('Normalized Value');
ylabel('Latitude');
title('Latitude Profile for Asia');
set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

% Save the figure as an SVG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Asia_lat.svg');
%% Europe
% First figure: Map of Europe
figure1 = figure;
custom_colormap = [1 1 1; addcolorplus(309)];
ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
setm(ax, 'FFaceColor', [1 1 1]); % Set background color
set(gcf, 'Colormap', custom_colormap);
R = georefcells([-90, 90], [-180, 180], size(livestock)); % Full global extent
geoshow(flipud(livestock), R, 'DisplayType', 'texturemap'); % Display livestock data
setm(ax, 'MapLatLimit', [35, 70], 'MapLonLimit', [-25, 45]); % Adjusted for Europe
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black');
caxis([0, 150]);
set(figure1, 'Position', [100, 100, 900, 500]); % Wider figure for the map

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Europe_modern.svg');

% Second figure: Latitude Profile Plot for Europe
figure2 = figure;
ydata = linspace(90, -90, 1800); % Adjust latitude range if specific to Europe
latitude_Europe_norm = latitude_Europe ./ sum(latitude_Europe); % Assuming `latitude_Europe` data
plot(latitude_Europe_norm, ydata, 'LineWidth', 1.5);
ylim([35, 70]);
xlabel('Normalized Value');
ylabel('Latitude');
title('Latitude Profile for Europe');
set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

% Save the figure as an SVG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/Europe_lat.svg');
%% North America
% First figure: Map of North America
figure1 = figure;
%custom_colormap = [1 1 1; addcolorplus(309)];
custom_colormap = [reds_to_white; white_to_greens];
ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
setm(ax, 'FFaceColor', [1 1 1]); % Set background color
set(gcf, 'Colormap', custom_colormap);

% Georeference and display the livestock data
R = georefcells([-90, 90], [0, 360], size(compare_turnover_rcp85)); % Full global extent
geoshow(flipud(compare_turnover_rcp85), R, 'DisplayType', 'texturemap'); % Display livestock data

% Set map limits to focus on North America
setm(ax, 'MapLatLimit', [15, 85], 'MapLonLimit', [-170, -50]); % Adjusted for North America

% Load and plot coastlines
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black');

% Set color scale and colorbar
hc = colorbar;
title(hc, 'heads (livestock units)/cell', 'FontSize', 12);
caxis([-14000, 14000]);

% Adjust the figure position and size
set(figure1, 'Position', [100, 100, 900, 500]); % Width 900, Height 500 for a wider map view

% Save the figure as an SVG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/NorthAmerica_modern.svg');

figure2 = figure;
ydata = linspace(90, -90, 160); % Latitude values covering entire globe for general consistency
latitude_NorthAmerica_norm = latitude_NorthAmerica ./ sum(latitude_NorthAmerica); % Assuming `latitude_NorthAmerica` data

% Plot normalized latitude profile for North America
plot(latitude_NorthAmerica_norm, ydata, 'LineWidth', 1.5);
% Set y-axis limits to focus on North America's latitude range
ylim([15, 85]);
xlabel('Normalized Value');
ylabel('Latitude');
title('Latitude Profile for North America');

% Adjust the figure position and size to have the same height as the map, but narrower width
set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

% Save the second figure as an SVG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/NorthAmerica_lat.svg');

%% For loop for all continents
% Define a struct array with region-specific parameters
regions = struct( ...
    'name', {'Africa', 'SouthAmerica', 'Oceania', 'Asia', 'Europe', 'NorthAmerica'}, ...
    'latLimit', {[-35, 40], [-60, 15], [-50, 10], [5, 55], [35, 70], [15, 85]}, ...
    'lonLimit', {[-20, 55], [-90, -30], [110, 180], [60, 150], [-25, 45], [-170, -50]}, ...
    'data', {'latitude_Africa', 'latitude_SouthAmerica', 'latitude_Oceania', 'latitude_Asia', 'latitude_Europe', 'latitude_NorthAmerica'} ...
);

% Loop through each region
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    latLimit = regions(i).latLimit;
    lonLimit = regions(i).lonLimit;
    latitudeData = eval(regions(i).data);  % Load latitude data specific to each region
    
    figure1 = figure;
    custom_colormap = [1 1 1; addcolorplus(309)];
    ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
    setm(ax, 'FFaceColor', [1 1 1]); % Set background color
    set(gcf, 'Colormap', custom_colormap);
    R = georefcells([-90, 90], [-180, 180], size(livestock)); % Full global extent
    geoshow(flipud(livestock), R, 'DisplayType', 'texturemap'); % Display livestock data

    % Set map limits based on region
    setm(ax, 'MapLatLimit', latLimit, 'MapLonLimit', lonLimit);
    
    % Load and plot coastlines
    load coastlines;
    plotm(coastlat, coastlon, 'Color', 'black');

     caxis([0, 150]);

    % Adjust the figure position and size
    set(figure1, 'Position', [100, 100, 900, 500]); % Width 900, Height 500

    % Save the map figure
    saveas(gcf, sprintf('/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/%s_modern.svg', regionName));

    figure2 = figure;
    ydata = linspace(90, -90, 1800); % Latitude values covering entire globe
    latitude_norm = latitudeData ./ sum(latitudeData); % Normalize latitude data

    % Plot normalized latitude profile
    plot(latitude_norm, ydata, 'LineWidth', 1.5);

    % Set y-axis limits to match the region's latitude range
    ylim(latLimit);
    xlabel('Normalized Value');
    ylabel('Latitude');
    title(sprintf('Latitude Profile for %s', regionName));

    % Adjust the figure position and size
    set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

    % Save the latitude profile figure
    saveas(gcf, sprintf('/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/%s_lat.svg', regionName));
end
%% For loop for all continents for future projections under ssp585
% Define a struct array with region-specific parameters
regions = struct( ...
    'name', {'Africa85', 'SouthAmerica85', 'Oceania85', 'Asia85', 'Europe85', 'NorthAmerica85'}, ...
    'latLimit', {[-35, 40], [-60, 15], [-50, 10], [5, 55], [35, 70], [15, 85]}, ...
    'lonLimit', {[-20, 55], [-90, -30], [110, 180], [60, 150], [-25, 45], [-170, -50]}, ...
    'data', {'latitude_Africa85', 'latitude_SouthAmerica85', 'latitude_Oceania85', 'latitude_Asia85', 'latitude_Europe85', 'latitude_NorthAmerica85'} ...
);

% Loop through each region
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    latLimit = regions(i).latLimit;
    lonLimit = regions(i).lonLimit;
    latitudeData = eval(regions(i).data);  % Load latitude data specific to each region
    
    figure1 = figure;
    custom_colormap = [1 1 1; addcolorplus(309)];
    ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
    setm(ax, 'FFaceColor', [1 1 1]); % Set background color
    set(gcf, 'Colormap', custom_colormap);
    R = georefcells([-90, 90], [0, 360], size(futureniche_struct.rcp85)); % Full global extent
    geoshow(flipud(futureniche_struct.rcp85), R, 'DisplayType', 'texturemap'); % Display livestock data

    % Set map limits based on region
    setm(ax, 'MapLatLimit', latLimit, 'MapLonLimit', lonLimit);
    
    % Load and plot coastlines
    load coastlines;
    plotm(coastlat, coastlon, 'Color', 'black');


    % Adjust the figure position and size
    set(figure1, 'Position', [100, 100, 900, 500]); % Width 900, Height 500

    % Save the map figure
    saveas(gcf, sprintf('/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/%s_modern.svg', regionName));

    figure2 = figure;
    ydata = linspace(90, -90, 160); % Latitude values covering entire globe
    latitude_norm = latitudeData ./ sum(latitudeData); % Normalize latitude data

    % Plot normalized latitude profile
    plot(latitude_norm, ydata, 'LineWidth', 1.5);

    % Set y-axis limits to match the region's latitude range
    ylim(latLimit);
    xlabel('Normalized Value');
    ylabel('Latitude');
    title(sprintf('Latitude Profile for %s', regionName));

    % Adjust the figure position and size
    set(figure2, 'Position', [1100, 100, 300, 500]); % Narrower figure, but same height as the map figure

    % Save the latitude profile figure
    saveas(gcf, sprintf('/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/%s_lat.svg', regionName));
end
%% longitude for all continents for modern
% Define a struct array with region-specific parameters for longitude profiles
regions = struct( ...
    'name', {'Africa', 'SouthAmerica', 'Oceania', 'Asia', 'Europe', 'NorthAmerica'}, ...
    'lonLimit', {[-20, 55], [-90, -30], [110, 180], [60, 150], [-25, 45], [-170, -50]}, ...
    'data', {'longitude_Africa', 'longitude_SouthAmerica', 'longitude_Oceania', 'longitude_Asia', 'longitude_Europe', 'longitude_NorthAmerica'} ...
);

% Loop through each region to generate longitude profile plots
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    lonLimit = regions(i).lonLimit;
    longitudeData = eval(regions(i).data);  % Load longitude data specific to each region

    % Create the longitude profile figure
    figure;
    xdata = linspace(-180, 180, 3600); % Longitude values covering the entire globe
    longitude_norm = longitudeData ./ sum(longitudeData); % Normalize longitude data

    % Plot normalized longitude profile
    plot(xdata, longitude_norm, 'LineWidth', 1.5);

    % Set x-axis limits to match the region's longitude range
    xlim(lonLimit);
    xlabel('Longitude');
    ylabel('Normalized Value');
    title(sprintf('Longitude Profile for %s', regionName));

    % Adjust the figure position and size to match previous specifications
    set(gcf, 'Position', [100, 100, 900, 300]); % Wider figure, but consistent height

    % Save the longitude profile figure as an SVG
    saveas(gcf, sprintf('/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/%s_lon.svg', regionName));
end
%% longitude for all continents future
% Define a struct array with region-specific parameters for longitude profiles
regions = struct( ...
    'name', {'Africa85', 'SouthAmerica85', 'Oceania85', 'Asia85', 'Europe85', 'NorthAmerica85'}, ...
    'lonLimit', {[-20, 55], [-90, -30], [110, 180], [60, 150], [-25, 45], [-170, -50]}, ...
    'data', {'longitude_Africa85', 'longitude_SouthAmerica85', 'longitude_Oceania85', 'longitude_Asia85', 'longitude_Europe85', 'longitude_NorthAmerica85'} ...
);

% Loop through each region to generate longitude profile plots
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    lonLimit = regions(i).lonLimit;
    longitudeData = eval(regions(i).data);  % Load longitude data specific to each region

    % Create the longitude profile figure
    figure;
    xdata = linspace(-180, 180, 320); % Longitude values covering the entire globe
    longitude_norm = longitudeData ./ sum(longitudeData); % Normalize longitude data

    % Plot normalized longitude profile
    plot(xdata, longitude_norm, 'LineWidth', 1.5);

    % Set x-axis limits to match the region's longitude range
    xlim(lonLimit);
    xlabel('Longitude');
    ylabel('Normalized Value');
    title(sprintf('Longitude Profile for %s', regionName));

    % Adjust the figure position and size to match previous specifications
    set(gcf, 'Position', [100, 100, 900, 300]); % Wider figure, but consistent height

    % Save the longitude profile figure as an SVG
    saveas(gcf, sprintf('/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/%s_lon.svg', regionName));
end
%% impact on population and livestock
% step 1: data preparation
load('turnover.mat')% this mat file includes these variables:compare_turnover_rcp26-compare_turnover_rcp85
% population count (population per pixel)
[pop_count, Rpop] = readgeoraster('gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_15_min.tif')
Rpop = georefcells([-90,90],[-180,180],size(pop_count));
pop_count(pop_count<0)=0
[resizedpop_count,resizedRpop] = georesize(pop_count,Rpop,1800/720,"bilinear");
% population per square kilometer
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
pop(pop<0)=0
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");

%% rescale population
% Step 2: Calculate the scaling factor to preserve the total population
[resizedpop_count,resizedRpop] = georesize(pop_count,Rpop,1800/720,"bilinear");
original_sum = sum(pop_count(:));  % Sum of the original data
resized_sum = sum(resizedpop_count(:));  % Sum of the resized data
scaling_factor = original_sum / resized_sum;

% Step 3: Apply the scaling factor to the resized data
resizedpop_count = double(resizedpop_count * scaling_factor);

% Verify if the new resized total matches the original total
new_resized_sum = sum(resizedpop_count(:));  % Should be close to original_sum
disp(['Original total population: ', num2str(original_sum)]);
disp(['Adjusted resized total population: ', num2str(new_resized_sum)]);

%% flip variables and check if they are aligned
turnover85_resize=imresize(compare_turnover_rcp85,[1800,3600])
imagesc(turnover85_resize)
turnover85_resize_flip = double(turnover85_resize(:, [ceil(end/2+1):end, 1:floor(end/2)]));
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
%grassland_pop(grassland_env>0)=resizedpop_count(grassland_env>0)
grassland_pop(livestock>0)=resizedpop_count(livestock>0)
grassland_pop(resizedpop>20)=0
grassland_pop_sum=sum(sum(grassland_pop))

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
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/population_positive.svg');
%%
%% rescale grassland influenced population
% Step 2: Calculate the scaling factor to preserve the total population
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
%% try for one
% Initialize pop_neg as an empty cell array
pop_neg = {};
pop_neg{1} = sum(sum(resized_grassland_pop_negative .* double(country_data.m_AFG)));

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
%% create bar chart for negative
% Convert pop_neg to a numeric vector (assuming all cells contain numeric values)
pop_neg_values = cell2mat(pop_neg);

% Create a bar chart with country names as labels
figure;
bar(pop_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for negative
% Sort pop_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(pop_neg_values, 'descend');
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
%% create bar chart for positive
% Convert pop_posi to a numeric vector (assuming all cells contain numeric values)
pop_posi_values = cell2mat(pop_posi);

% Create a bar chart with country names as labels
figure;
bar(pop_posi_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for positive
% Sort pop_posi_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(pop_posi_values, 'descend');
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
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
%% try for one
% Initialize cattle_neg as an empty cell array
cattle_neg = {};
cattle_neg{1} = sum(sum(resized_cattle_negative .* double(country_data.m_AFG)));

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
%% create bar chart for negative cattle
% Convert cattle_neg to a numeric vector (assuming all cells contain numeric values)
cattle_neg_values = cell2mat(cattle_neg);

% Create a bar chart with country names as labels
figure;
bar(cattle_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for negative cattle
% Sort cattle_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(cattle_neg_values, 'descend');
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
%% Get for all countries for positive cattle
% Initialize an empty cell array to store results
cattle_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in cattle_neg
    cattle_neg{i} = sum(sum(resized_cattle_positive .* double(country_data.(fieldName))));
end
%% create bar chart for positive cattle
% Convert cattle_neg to a numeric vector (assuming all cells contain numeric values)
cattle_neg_values = cell2mat(cattle_neg);

% Create a bar chart with country names as labels
figure;
bar(cattle_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for positive cattle
% Sort cattle_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(cattle_neg_values, 'descend');
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
%% %%%%%%% sheep
%% start calculation
% This calculates the population that is inhabited near grasslands
sheep_pop=zeros(size(resizesheep))
%sheep(grassland_env>0)=resizedpop_count(grassland_env>0)
sheep(livestock>0)=resizedpop_count(livestock>0)
sheep(resizedpop>20)=0
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
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
%% try for one
% Initialize sheep_neg as an empty cell array
sheep_neg = {};
sheep_neg{1} = sum(sum(resized_sheep_negative .* double(country_data.m_AFG)));

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
%% create bar chart for negative sheep
% Convert sheep_neg to a numeric vector (assuming all cells contain numeric values)
sheep_neg_values = cell2mat(sheep_neg);

% Create a bar chart with country names as labels
figure;
bar(sheep_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for negative sheep
% Sort sheep_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(sheep_neg_values, 'descend');
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
%% Get for all countries for positive sheep
% Initialize an empty cell array to store results
sheep_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in sheep_neg
    sheep_neg{i} = sum(sum(resized_sheep_positive .* double(country_data.(fieldName))));
end
%% create bar chart for positive sheep
% Convert sheep_neg to a numeric vector (assuming all cells contain numeric values)
sheep_neg_values = cell2mat(sheep_neg);

% Create a bar chart with country names as labels
figure;
bar(sheep_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for positive sheep
% Sort sheep_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(sheep_neg_values, 'descend');
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
%% %%%%%% goats
%% start calculation
% This calculates the population that is inhabited near grasslands
goats_pop=zeros(size(resizegoats))
%goats(grassland_env>0)=resizedpop_count(grassland_env>0)
goats(livestock>0)=resizedpop_count(livestock>0)
goats(resizedpop>20)=0
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
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
title(hc, 'people/cell', 'FontSize', 14)
plotm(coastlat, coastlon, 'Color', 'black'); 
caxis([0,150])
set(gcf, 'Position',  [361,614,899,295])
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
%% try for one
% Initialize goats_neg as an empty cell array
goats_neg = {};
goats_neg{1} = sum(sum(resized_goats_negative .* double(country_data.m_AFG)));

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
%% create bar chart for negative goats
% Convert goats_neg to a numeric vector (assuming all cells contain numeric values)
goats_neg_values = cell2mat(goats_neg);

% Create a bar chart with country names as labels
figure;
bar(goats_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for negative goats
% Sort goats_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(goats_neg_values, 'descend');
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
%% Get for all countries for positive goats
% Initialize an empty cell array to store results
goats_neg = {};

% Get all field names in the struct `country_data`
countryFields = fieldnames(country_data);

% Loop through each country field
for i = 1:length(countryFields)-1
    % Get the name of the current country field
    fieldName = countryFields{i};
    
    % Calculate the result for the current country and store it in goats_neg
    goats_neg{i} = sum(sum(resized_goats_positive .* double(country_data.(fieldName))));
end
%% create bar chart for positive goats
% Convert goats_neg to a numeric vector (assuming all cells contain numeric values)
goats_neg_values = cell2mat(goats_neg);

% Create a bar chart with country names as labels
figure;
bar(goats_neg_values);
title('Population in Grassland by Country');
xlabel('Country');
ylabel('Population Count');
set(gca, 'XTick', 1:length(countryFields), 'XTickLabel', countryFields);
%% in descending order for positive goats
% Sort goats_neg_values in descending order and get the sorting indices
[sorted_values, sort_index] = sort(goats_neg_values, 'descend');
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