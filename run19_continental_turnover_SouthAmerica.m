%%
load('livestockdensity.mat')
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('/Users/lichaohui/Desktop/calculation/grazingniche/matdata/aggregate_tas.mat')
aggregate_tas=aggregate_data;

aggregate_sfcWind(aggregate_sfcWind>10)=0
load("grassland_env.mat")
%% 提前准备：pixel to area calculation [function]
% this will give you a matlab standardized function that allows you to
% calculate any size matrix into earth area size. 
% you only need to insert your "mapsize" into the first line of command and
% then will be able to obtain a pixel to area conversion matrix "areas"
% (latitude only) and "worldareas" variable (or whatever you wish to call
% it) full size conversion matrix.

mapsize=[160,320]

R = 6371; % Earth's radius in km
latitude_diff = 180 / mapsize(1); % in degrees
longitude_diff = 360 / mapsize(2); % in degrees, but note that it remains constant around the globe

% Convert longitude difference to radians
longitude_diff = deg2rad(longitude_diff);

% Initialize the areas matrix
areas = zeros(mapsize(1), 1);

% Compute area for each latitude band
for i = 1:mapsize(1)
    % Calculate the latitude bounds for the current band
    latitude1 = 90 - (i-1) * latitude_diff; % Start from the top and move downwards
    latitude2 = 90 - i * latitude_diff;
    
    % Convert to radians
    latitude1 = deg2rad(latitude1);
    latitude2 = deg2rad(latitude2);
    
    % Compute the area and store in the matrix
    areas(i) = R^2 * longitude_diff * (sin(latitude1) - sin(latitude2));
end

% Display the areas matrix
disp(areas);

worldarea=repmat(areas,[1,mapsize(2)])

sum(sum(worldarea))
%% 提前准备：code for the continent mask
%% now we also need to convert our continent masks into the same matrix size.
% the result is mask_futures

% Define path to your file
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
    [resized_data, ~] = georesize(country_data_raw, Rmask, 320/720, "bilinear");
    
    % Store in the structure
    country_data.(country_name) = resized_data';
end


% mask for continents
% asia
Asia = {'m_AFG', 'm_ARM', 'm_AZE', 'm_BHR', 'm_BGD', 'm_BTN', 'm_BRN', 'm_KHM', 'm_CHN', 'm_CYM',... 
        'm_GEO', 'm_IND', 'm_IDN', 'm_IRN', 'm_IRQ', 'm_ISR', 'm_JPN', 'm_JOR', 'm_KAZ', 'm_KWT', ...
        'm_KGZ', 'm_LAO', 'm_LBN', 'm_MYS', 'm_MNG', 'm_MMR', 'm_NPL', 'm_PRK', 'm_OMN', ...
        'm_PAK', 'm_PHL', 'm_QAT', 'm_RUS', 'm_SAU', 'm_SGP', 'm_KOR', 'm_LKA', 'm_SYR', 'm_TWN', ...
        'm_TJK', 'm_THA', 'm_TUR', 'm_TKM', 'm_ARE', 'm_UZB', 'm_VNM', 'm_YEM'};

Asiamask = country_data.m_AFG;
for i = 2:length(Asia)
    Asiamask = Asiamask + country_data.(Asia{i});
end

imagesc(Asiamask)
colorbar

Asiamask = Asiamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% for europe
Europe = {'m_ALB', 'm_AND', 'm_AUT', 'm_BEL', 'm_BIH', 'm_BGR', 'm_HRV', 'm_CYP', 'm_CZE', 'm_DNK',... 
          'm_EST', 'm_FIN', 'm_FRA', 'm_DEU', 'm_GRC', 'm_HUN', 'm_ISL', 'm_IRL', 'm_ITA', ...
          'm_LVA', 'm_LTU', 'm_LUX', 'm_MLT', 'm_MDA', 'm_MNE', 'm_NLD', 'm_MKD', ...
          'm_NOR', 'm_POL', 'm_PRT', 'm_ROU', 'm_SRB', 'm_SVK', 'm_SVN', 'm_ESP', ...
          'm_SWE', 'm_CHE', 'm_UKR', 'm_GBR'};

Europemask = country_data.m_ALB;
for i = 2:length(Europe)
    Europemask = Europemask + country_data.(Europe{i});
end

imagesc(Europemask)
colorbar

Europemask = Europemask(:, [ceil(end/2+1):end, 1:floor(end/2)]);


% now for south america
SouthAmerica = {'m_ARG', 'm_BOL', 'm_BRA', 'm_CHL', 'm_COL', 'm_ECU', 'm_GUY', 'm_PRY', 'm_PER', 'm_SUR', 'm_URY', 'm_VEN'};
SouthAmericamask = country_data.m_ARG;
for i = 2:length(SouthAmerica)
    SouthAmericamask = SouthAmericamask + country_data.(SouthAmerica{i});
end
imagesc(SouthAmericamask)
colorbar

SouthAmericamask = SouthAmericamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% north america
NorthAmerica = {'m_ATG', 'm_BHS', 'm_BRB', 'm_BEL', 'm_CAN', 'm_CYM', 'm_CRI', 'm_CUB', 'm_DMA', 'm_DOM', 'm_SLV', 'm_GRL', 'm_GRD', 'm_GLP', 'm_GTM', 'm_HND', 'm_JAM', 'm_MEX', 'm_MTQ', 'm_NIC', 'm_PAN', 'm_PRI', 'm_SPM', 'm_LCA', 'm_VCT', 'm_TTO', 'm_USA', 'm_VIR'};
NorthAmericamask = country_data.m_ATG;
for i = 2:length(NorthAmerica)
    NorthAmericamask = NorthAmericamask + country_data.(NorthAmerica{i});
end
imagesc(NorthAmericamask)
colorbar
% SouthAmerica
SouthAmerica = {'m_DZA', 'm_COD','m_SDN', 'm_AGO', 'm_BEN', 'm_BWA', 'm_BFA', 'm_BDI', 'm_CMR', 'm_CPV', 'm_CAF', 'm_TCD', 'm_COM', 'm_COG', 'm_CIV', 'm_DJI', 'm_EGY', 'm_GNQ', 'm_ERI', 'm_ETH', 'm_GAB', 'm_GMB', 'm_GHA', 'm_GIN', 'm_GNB', 'm_KEN', 'm_LSO', 'm_LBR', 'm_LBY', 'm_MDG', 'm_MWI', 'm_MLI', 'm_MRT', 'm_MUS', 'm_MAR', 'm_MOZ', 'm_NAM', 'm_NER', 'm_NGA', 'm_REU', 'm_RWA', 'm_STP', 'm_SEN', 'm_SLE', 'm_SOM', 'm_ZAF', 'm_SSD', 'm_ESH', 'm_SWZ', 'm_TZA', 'm_TGO', 'm_TUN', 'm_UGA', 'm_ZMB', 'm_ZWE'};
Africamask = country_data.m_DZA;
for i = 2:length(SouthAmerica)
    Africamask = Africamask + country_data.(SouthAmerica{i});
end
imagesc(Africamask)
colorbar

Africamask = Africamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% Oceania
Oceania = {'m_AUS', 'm_FJI', 'm_KIR', 'm_FSM', 'm_NCL', 'm_NZL', 'm_NIU', 'm_PLW', 'm_PNG', 'm_WSM', 'm_SLB', 'm_TON', 'm_VUT'};

Oceaniamask = country_data.m_AUS;
for i = 2:length(Oceania)
    Oceaniamask = Oceaniamask + country_data.(Oceania{i});
end
imagesc(Oceaniamask)
colorbar

Oceaniamask = Oceaniamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% convert the masks to double so that it's easier to manipulate them later
% then not only can they be used as logical layers we can directly multiply
% them with matricies, so we don't have to resort to too many conditions
% because the 0 is able to cancle any variable we want

Asiamask=double(Asiamask);
Europemask=double(Europemask);
Oceaniamask=double(Oceaniamask);
Africamask=double(Africamask);
NorthAmericamask=double(NorthAmericamask);
SouthAmericamask=double(SouthAmericamask);

%% 决定SouthAmerica Niche
cond1_pr = (aggregate_pr >= 228 & aggregate_pr <= 1889);
cond1_tas = (aggregate_tas >= 3.1 & aggregate_tas <= 26);
cond1_sfcWind = (aggregate_sfcWind >= 1.2 & aggregate_sfcWind <= 5.3);
cond_hurs = (aggregate_hurs >= 43 & aggregate_hurs <= 63);
cond1_landuse = (grassland_env>0);


% Combine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% Create the coupled landuse map
landuse_coupled1_average = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,90],[-180,180],size(landuse_coupled1_average));

figure1 = figure('WindowState','fullscreen');
sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap([-60, 15],[-90, -30]);
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(landuse_coupled1_average), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_general_average.svg');
close all;

%save('./grazingniche/matdata/niche_average.mat',"landuse_coupled1_average")

%% 决定futuer niche: average
load("future_climate.mat")
load("futureland.mat")
% Define the list of scenarios
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};

% Initialize a structure to store the niche maps
futureniche_struct_average = struct();

% Loop through each scenario
for s = 1:length(scenarios)
    % Get the scenario
    scenario = scenarios{s};

    % Get the data for the scenario
    data_pr = eval(['data_pr_' scenario]);
    data_tas = eval(['data_tas_' scenario]);
    data_sfcWind = eval(['data_sfcWind_' scenario]);
    data_hurs = eval(['data_hurs_' scenario]);
    data_landuse = eval(['futureland_' scenario]);
%     imagesc(data_landuse) %发现数据不一致
%     imagesc(data_tas)
%     imagesc(data_landuse_new)


% 调整数据
    halfSize = size(data_landuse, 2) / 2;
    data_landuse_new = circshift(data_landuse, [0, halfSize]);
 
        % Define the conditions

    cond_pr = (data_pr >= 228 & data_pr <= 1889);
    cond_tas = (data_tas >= 3.1 & data_tas <= 26);
    cond_sfcWind = (data_sfcWind >= 1.2 & data_sfcWind <= 5.3);
    cond_hurs = (data_hurs >= 43 & data_hurs <= 63);
    cond_landuse = (data_landuse_new >= 0.01);


    % Combine the conditions to create the niche map
    all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;


    futureniche = data_landuse_new .* all_conditions;

    % Store the niche map in the structure
    futureniche_struct_average.(scenario) = futureniche;

    imagesc(futureniche_struct_average.rcp26)
    imagesc(data_landuse_new)
    imagesc(data_landuse)
    imagesc(all_conditions)

    % Display the niche map
subplot(2, 2, s);
R=georefcells([-90,90],[0,360],size(futureniche));    
sgtitle('Future livestock niche', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(332)];
ax = worldmap([-60, 15],[-90, -30]);
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
    title(['GN distribution under ' scenario]);
    sgtitle('Future grazing niche under scenarios' )
end

set(gcf, 'Position',  [259,455,1140,528])

%saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_average.svg');

%save('./grazingniche/matdata/futureniche_average.mat','futureniche_struct_average')
%% turnover rate suggested by Max
%% Step 1.Find the relevant variables
% 要找到nicheland_world和nicheland_world_rcp26的根变量，这两个已经是有面积的了，这样再去插值升尺度或者降尺度就麻烦了。
%nicheland_world的根变量是landuse_coupled1
%nicheland_world_rcp26的根变量是futureniche_struct.rcp26

% imagesc(landuse_coupled1_average) %美国在左边，值在100以内，1800*3600
% colorbar
% imagesc(futureniche_struct.rcp26) %美国在右边，值在1以内，160*320
% colorbar

%% Step 2.Unifying the resolution, unit, also unifying the coordinates 
% Unifying coordination system
halfSize = size(landuse_coupled1_average, 2) / 2;
landuse_coupled1_resize1_average = circshift(landuse_coupled1_average, [0, halfSize]);
% Unifying unit (into percentage of grassland per pixel, 0-1)
landuse_coupled1_resize2_average=landuse_coupled1_resize1_average/100
% Unifying resolution
landuse_coupled1_resize3_average = imresize(landuse_coupled1_resize2_average, [160,320], 'bilinear');
%% Test if they are unified
% imagesc(landuse_coupled1_resize3_average) %美国在右边，值0-1，160*320
% colorbar
% 
% imagesc(futureniche_struct.rcp26) %美国在右边，值0-1，160*320
% colorbar

%% Step 3.Turn percentage of area per pixel into area per pixel

compare_futureniche_rcp26_average=futureniche_struct_average.rcp26.*worldarea;
compare_futureniche_rcp45_average=futureniche_struct_average.rcp45.*worldarea;
compare_futureniche_rcp60_average=futureniche_struct_average.rcp60.*worldarea;
compare_futureniche_rcp85_average=futureniche_struct_average.rcp85.*worldarea;

compare_niche=landuse_coupled1_resize3_average.*worldarea

sum(sum(compare_niche))

%% Step 4. Compare them
compare_turnover_rcp26_average=compare_futureniche_rcp26_average-compare_niche
compare_turnover_rcp45_average=compare_futureniche_rcp45_average-compare_niche
compare_turnover_rcp60_average=compare_futureniche_rcp60_average-compare_niche
compare_turnover_rcp85_average=compare_futureniche_rcp85_average-compare_niche

save('./grazingniche/matdata/turnover_average.mat',"compare_turnover_rcp26_average","compare_turnover_rcp45_average","compare_turnover_rcp60_average","compare_turnover_rcp85_average")
load('turnover_average.mat')% this .mat contains compare_turnover_rcp85 and other scenarios variable
imagesc(compare_niche) 
colorbar

imagesc(compare_futureniche_rcp26_average)
colorbar

sum(sum(compare_niche))%结果是7.5e6
sum(sum(compare_futureniche_rcp26_average))%结果是9.1e6
sum(sum(compare_futureniche_rcp85_average))%结果是8.6e6

% 所以最后的结果是未来的niche会增加？？？？还增加了挺多？？
% 是不是算错了？？？？

%% Step 5: Calculate turnover rate for each continent

%% for each continent
% This is the net value. Net turnover rate
turnover_total_rcp26_average=sum(sum(compare_turnover_rcp26_average));

% How much will be gained in new areas?
turnover_increase_rcp26_average=compare_turnover_rcp26_average;
turnover_increase_rcp26_average(turnover_increase_rcp26_average<0)=0
turnover_increase_value_rcp26_average=sum(sum(turnover_increase_rcp26_average))

% How much will be lost?
turnover_decrease_rcp26_average=compare_turnover_rcp26_average;
turnover_decrease_rcp26_average(turnover_decrease_rcp26_average>0)=0
turnover_decrease_value_rcp26_average=sum(sum(turnover_decrease_rcp26_average))

% How much in each continent?
turnover_SouthAmerica_rcp26_average=sum(sum(compare_turnover_rcp26_average(Africamask==1)));
turnover_increase_SouthAmerica_rcp26_average=sum(sum(turnover_increase_rcp26_average((Africamask==1))))
turnover_decrease_SouthAmerica_rcp26_average=sum(sum(turnover_decrease_rcp26_average((Africamask==1))))

%% Batch making for all continents and all rcps
% Define RCPs and Continents
rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
continents = {'Asia', 'Europe', 'Oceania', 'Africa', 'NorthAmerica', 'SouthAmerica'};

% Initialize structures to store results
turnover_total = struct();
turnover_increase_value = struct();
turnover_decrease_value = struct();
turnover_continents_average = struct();

% Loop over each RCP scenario
for i = 1:length(rcps)
    rcp = rcps{i};
    compare_turnover = eval(['compare_turnover_' rcp,'_average']); % Dynamically load the variable

    % Calculate total turnover for the current RCP
    turnover_total.(rcp) = sum(sum(compare_turnover));

    % Calculate turnover increase and decrease
    turnover_increase = compare_turnover;
    turnover_increase(turnover_increase < 0) = 0;
    turnover_increase_value.(rcp) = sum(sum(turnover_increase));

    turnover_decrease = compare_turnover;
    turnover_decrease(turnover_decrease > 0) = 0;
    turnover_decrease_value.(rcp) = sum(sum(turnover_decrease));

    % Calculate turnover for each continent
    for j = 1:length(continents)
        continent = continents{j};
        mask = eval([continent 'mask']); % Dynamically use the mask variable

        turnover_continents_average.(continent).(rcp).total = sum(sum(compare_turnover(mask == 1)));
        turnover_continents_average.(continent).(rcp).increase = sum(sum(turnover_increase(mask == 1)));
        turnover_continents_average.(continent).(rcp).decrease = sum(sum(turnover_decrease(mask == 1)));
    end
end

% Display results (optional)
disp(turnover_total);
disp(turnover_increase_value);
disp(turnover_decrease_value);
disp(turnover_continents_average);

%save('./grazingniche/matdata/turnover_continents_average.mat','turnover_decrease_SouthAmerica_rcp26_average')
%% Creat figures
bar([turnover_continents_average.SouthAmerica.rcp26.increase,turnover_continents_average.SouthAmerica.rcp26.decrease,turnover_continents_average.SouthAmerica.rcp45.increase,turnover_continents_average.SouthAmerica.rcp45.decrease,turnover_continents_average.SouthAmerica.rcp60.increase,turnover_continents_average.SouthAmerica.rcp60.decrease,turnover_continents_average.SouthAmerica.rcp85.increase,turnover_continents_average.SouthAmerica.rcp85.decrease])
%%
% Data for turnover
data = [
    turnover_continents_average.SouthAmerica.rcp26.increase, turnover_continents_average.SouthAmerica.rcp26.decrease;
    turnover_continents_average.SouthAmerica.rcp45.increase, turnover_continents_average.SouthAmerica.rcp45.decrease;
    turnover_continents_average.SouthAmerica.rcp60.increase, turnover_continents_average.SouthAmerica.rcp60.decrease;
    turnover_continents_average.SouthAmerica.rcp85.increase, turnover_continents_average.SouthAmerica.rcp85.decrease
];

% Bar chart
figure;
bar(data, 'grouped');

% X-axis labels
set(gca, 'XTickLabel', {'RCP 2.6', 'RCP 4.5', 'RCP 6.0', 'RCP 8.5'});

% Labels and title
xlabel('RCP Scenarios');
ylabel('Turnover for CNG (m^2)');
title('SouthAmerica');

% Add legend
legend({'Increase', 'Decrease'}, 'Location', 'northwest');

% Enhance appearance (optional)
grid on;
box on;

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_allrcp_bar_SouthAmerica.svg');

%%
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

%% Four rcps one plot
% Define the RCP scenarios
rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
titles = {'Turnover for RCP 2.6', 'Turnover for RCP 4.5', 'Turnover for RCP 6.0', 'Turnover for RCP 8.5'};

% Create a figure for the subplots
figure;
sgtitle('Turnover Across Different RCP Scenarios', 'FontSize', 16);

% Loop over each RCP scenario to create subplots
for i = 1:length(rcps)
    rcp = rcps{i};
    title_text = titles{i};
    
    % Create a subplot for each RCP scenario
    subplot(2, 2, i);
    
    % Set up the world map for the current subplot
ax = worldmap([-60, 15],[-90, -30]);
    setm(ax, 'FFaceColor', [1 1 1]);  % Set the map's background color
    
    % Apply the custom colormap
    set(gcf, 'Colormap', custom_colormap);
    
    % Set color axis limits to ensure zero aligns with white
    caxis([-14000 14000]);  % Adjust according to your data's range
    
    % Load and display the data for the current RCP
    compare_turnover = eval(['compare_turnover_' rcp '_average']); % Dynamically load the variable
    R = georefcells([-90, 90], [0, 360], size(compare_turnover));
    geoshow(flipud(compare_turnover), R, 'DisplayType', 'texturemap');
    
    % Load and plot coastlines
    load coastlines;
    plotm(coastlat, coastlon, 'Color', 'black');
    
    % Add a title to each subplot
    title(title_text, 'FontSize', 12);
    
    % Adjust subplot for tight and compact style
end

% Adjust the figure's size for a tight, compact display
set(gcf, 'Position', [100, 100, 1000, 600]);

% Display a colorbar that applies to all subplots
% Display a single colorbar that applies to all subplots
hc = colorbar('Position', [0.93, 0.1, 0.02, 0.8]);  % Custom position for a single colorbar
hc.Label.String = 'Decrease/Increase CNG (m^2)';
hc.Label.FontSize = 14;
hc.TicksMode = 'auto';  % Allow automatic ticks based on caxis range
hc.TickDirection = 'out';  % Ensure ticks are outward-facing

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_allrcp_SouthAmerica.svg');
%close all;
save('./grazingniche/matdata/turnover_average_SouthAmerica.mat','compare_turnover_rcp26_average','compare_turnover_rcp45_average','compare_turnover_rcp60_average','compare_turnover_rcp85_average');     


