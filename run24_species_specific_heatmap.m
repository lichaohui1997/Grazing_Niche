%% cattle

cattle_resize=imresize(resizecattle, [90 144]);
pr_resize=imresize(aggregate_pr, [90 144]);
tas_resize=imresize(aggregate_tas, [90 144]);
hurs_resize=imresize(aggregate_hurs, [90 144]);
sfcWind_resize=imresize(aggregate_sfcWind, [90 144]);

cattle_vector=reshape(cattle_resize,[],1)
pr_cattle_vector=pr_resize(:)
tas_cattle_vector=tas_resize(:)
hurs_cattle_vector=hurs_resize(:)
sfcWind_cattle_vector=sfcWind_resize(:)

pr_cattle_vector(cattle_vector<=2,:)=[]
tas_cattle_vector(cattle_vector<=2,:)=[]
sfcWind_cattle_vector(cattle_vector<=2,:)=[]
hurs_cattle_vector(cattle_vector<=2,:)=[]

%% cattle等高线图
subplot(2,3,1)
X=pr_cattle_vector;
Y=tas_cattle_vector;
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
title({'Cattle Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
xlim([0,10000])
defualtAxes()

subplot(2,3,2)
X=hurs_cattle_vector;
Y=tas_cattle_vector;
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
title({'Cattle Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,100])

subplot(2,3,3)
X=sfcWind_cattle_vector;
Y=tas_cattle_vector;
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
title({'Cattle Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,10])

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayheat_cattle.svg');
%% sheep
sheep_resize=imresize(resizesheep, [90 144]);
pr_resize=imresize(aggregate_pr, [90 144]);
tas_resize=imresize(aggregate_tas, [90 144]);
hurs_resize=imresize(aggregate_hurs, [90 144]);
sfcWind_resize=imresize(aggregate_sfcWind, [90 144]);

sheep_vector=reshape(sheep_resize,[],1)
pr_sheep_vector=pr_resize(:)
tas_sheep_vector=tas_resize(:)
hurs_sheep_vector=hurs_resize(:)
sfcWind_sheep_vector=sfcWind_resize(:)

pr_sheep_vector(sheep_vector<=2,:)=[]
tas_sheep_vector(sheep_vector<=2,:)=[]
sfcWind_sheep_vector(sheep_vector<=2,:)=[]
hurs_sheep_vector(sheep_vector<=2,:)=[]

%% sheep等高线图
subplot(2,3,1)
X=pr_sheep_vector;
Y=tas_sheep_vector;
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
title({'Sheep Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
xlim([0,10000])
defualtAxes()

subplot(2,3,2)
X=hurs_sheep_vector;
Y=tas_sheep_vector;
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
title({'Sheep Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,100])

subplot(2,3,3)
X=sfcWind_sheep_vector;
Y=tas_sheep_vector;
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
title({'Sheep Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,10])

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayheat_sheep.svg');
%% goats
goats_resize=imresize(resizegoats, [90 144]);
pr_resize=imresize(aggregate_pr, [90 144]);
tas_resize=imresize(aggregate_tas, [90 144]);
hurs_resize=imresize(aggregate_hurs, [90 144]);
sfcWind_resize=imresize(aggregate_sfcWind, [90 144]);

goats_vector=reshape(goats_resize,[],1)
pr_goats_vector=pr_resize(:)
tas_goats_vector=tas_resize(:)
hurs_goats_vector=hurs_resize(:)
sfcWind_goats_vector=sfcWind_resize(:)

pr_goats_vector(goats_vector<=2,:)=[]
tas_goats_vector(goats_vector<=2,:)=[]
sfcWind_goats_vector(goats_vector<=2,:)=[]
hurs_goats_vector(goats_vector<=2,:)=[]

%% goats等高线图
subplot(2,3,1)
X=pr_goats_vector;
Y=tas_goats_vector;
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
title({'Goats Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
xlim([0,10000])
defualtAxes()

subplot(2,3,2)
X=hurs_goats_vector;
Y=tas_goats_vector;
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
title({'Goats Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,100])

subplot(2,3,3)
X=sfcWind_goats_vector;
Y=tas_goats_vector;
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
title({'Goats Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])
xlim([0,10])

% Save the figure as a SVG
set(gcf, 'Position',  [751,163,1092,753])
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayheat_goats.svg');
