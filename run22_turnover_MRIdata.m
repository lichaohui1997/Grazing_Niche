%% This script calculates the turnover rate using the new niche
%% first decide the niche region now
%% couple the niche with the land use data
% in order to run this section you need to add the command file into
% directory in order for addcolorplus to work
% this section produces the main graph of the paper. A map of the grazing
% niche. 
% Here I am not using the nasa land use data but the Earthenv data. NASA is
% a catagorical land use dataset but Earthenv is a precentage map. It shows
% the percentage of shrublands and herbacesous vegetation cover. 
load('grassland_env.mat')

imagesc(grassland_env)
grassland_env = grassland_env(:, [ceil(end/2+1):end, 1:floor(end/2)]);
grassland_env=imresize(grassland_env,[160,320])

imagesc(data_sfcWind_rcp85_present)

%cond1_pr = (aggregate_pr >= 433 & aggregate_pr <= 2368);
cond1_pr = (data_pr_rcp85_present >= 100 & data_pr_rcp85_present <= 2300);
cond1_tas = (data_tas_rcp85_present >= -2 & data_tas_rcp85_present <= 28);
cond1_sfcWind = (data_sfcWind_rcp85_present >= 1.4 & data_sfcWind_rcp85_present <= 6);
cond_hurs = (data_hurs_rcp85_present >= 30 & data_hurs_rcp85_present <= 87);
cond1_landuse = (grassland_env>0);

% Combine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% incase you want to see what is holding you back, you can use this code to
% see
%imagesc(cond1_pr)

% %here I am trying to see the numerical values of my resuls
% niche_gridnum=numel(find(niche1(:)==1));
% land_gridnum=numel(find(landuse(:)>0));
% pasture_gridnum=numel(find((landuse(:)>5 & landuse(:)<=10)));
% pasture_prct=pasture_gridnum/land_gridnum;
% niche_prct1=niche_gridnum/pasture_gridnum;

% Create the coupled landuse map
landuse_coupled1 = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,90],[0,360],size(landuse_coupled1));

figure1 = figure('WindowState','fullscreen');
sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_general.svg');

save('./grazingniche/matdata/niche.mat',"landuse_coupled1")
%% Next decide the future niche
%% 决定niche 并graph niche
% Define the list of scenarios
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};

% Initialize a structure to store the niche maps
futureniche_struct = struct();

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

%     imagesc(data_landuse) 
%     imagesc(data_tas)
%     imagesc(data_landuse_new)

% 调整数据
    halfSize = size(data_landuse, 2) / 2;
    data_landuse_new = circshift(data_landuse, [0, halfSize]);
 
        % Define the conditions
    %cond_pr = (data_pr >= 433 & data_pr <= 2368);
    cond_pr = (data_pr >= 100 & data_pr <= 2300);
    cond_tas = (data_tas >= -2 & data_tas <= 28);
    cond_sfcWind = (data_sfcWind >= 1.4 & data_sfcWind <= 6);
    cond_hurs = (data_hurs >= 30 & data_hurs <= 80);
    cond_landuse = (data_landuse_new >= 0.01);


    % Combine the conditions to create the niche map
    all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;


    futureniche = data_landuse_new .* all_conditions;

    % Store the niche map in the structure
    futureniche_struct.(scenario) = futureniche;

    imagesc(futureniche_struct.rcp26)
    imagesc(data_landuse_new)
    imagesc(data_landuse)
    imagesc(all_conditions)

    % Display the niche map
subplot(2, 2, s);
R=georefcells([-90,90],[0,360],size(futureniche));    
sgtitle('Future livestock niche', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(332)];
ax = worldmap('world')      
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

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche.svg');

save('./grazingniche/matdata/futureniche.mat','futureniche_struct')

%% turnover rate suggested by Max
%% Step 1.Find the relevant variables
% 找到nicheland_world和nicheland_world_rcp26的根变量
%nicheland_world的根变量是landuse_coupled1
%nicheland_world_rcp26的根变量是futureniche_struct.rcp26
load('./grazingniche/matdata/niche.mat')
load('./grazingniche/matdata/futureniche.mat')
nexttile
imagesc(landuse_coupled1) 
colorbar
nexttile
imagesc(futureniche_struct.rcp26) 
colorbar


%% Step 3.Turn percentage of area per pixel into area per pixel

compare_futureniche_rcp26=futureniche_struct.rcp26.*worldarea;
compare_futureniche_rcp45=futureniche_struct.rcp45.*worldarea;
compare_futureniche_rcp60=futureniche_struct.rcp60.*worldarea;
compare_futureniche_rcp85=futureniche_struct.rcp85.*worldarea;

compare_niche=landuse_coupled1.*worldarea/100

%% Step 4. Compare them
compare_turnover_rcp26=compare_futureniche_rcp26-compare_niche
compare_turnover_rcp45=compare_futureniche_rcp45-compare_niche
compare_turnover_rcp60=compare_futureniche_rcp60-compare_niche
compare_turnover_rcp85=compare_futureniche_rcp85-compare_niche

save('./grazingniche/matdata/turnover.mat',"compare_turnover_rcp26","compare_turnover_rcp45","compare_turnover_rcp60","compare_turnover_rcp85")
load('turnover.mat')% this .mat contains compare_turnover_rcp85 and other scenarios variable
imagesc(compare_niche) 
colorbar

imagesc(compare_futureniche_rcp26)
colorbar

sum(sum(compare_niche))%结果是7.5e6
sum(sum(compare_futureniche_rcp26))%结果是9.1e6
sum(sum(compare_futureniche_rcp85))%结果是8.6e6



%% Step 5: Calculate turnover rate for each continent

% This is the net value. Net turnover rate
turnover_total_rcp26=sum(sum(compare_turnover_rcp26));

% How much will be gained in new areas?
turnover_increase_rcp26=compare_turnover_rcp26;
turnover_increase_rcp26(turnover_increase_rcp26<0)=0
turnover_increase_value_rcp26=sum(sum(turnover_increase_rcp26))

% How much will be lost?
turnover_decrease_rcp26=compare_turnover_rcp26;
turnover_decrease_rcp26(turnover_decrease_rcp26>0)=0
turnover_decrease_value_rcp26=sum(sum(turnover_decrease_rcp26))

% How much in each continent?
turnover_africa_rcp26=sum(sum(compare_turnover_rcp26(Africamask==1)));
turnover_increase_africa_rcp26=sum(sum(turnover_increase_rcp26((Africamask==1))))
turnover_decrease_africa_rcp26=sum(sum(turnover_decrease_rcp26((Africamask==1))))

%% Batch making for all continents and all rcps
% Define RCPs and Continents
rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
continents = {'Asia', 'Europe', 'Oceania', 'Africa', 'NorthAmerica', 'SouthAmerica'};

% Initialize structures to store results
turnover_total = struct();
turnover_increase_value = struct();
turnover_decrease_value = struct();
turnover_continents = struct();

% Loop over each RCP scenario
for i = 1:length(rcps)
    rcp = rcps{i};
    compare_turnover = eval(['compare_turnover_' rcp]); % Dynamically load the variable

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

        turnover_continents.(continent).(rcp).total = sum(sum(compare_turnover(mask == 1)));
        turnover_continents.(continent).(rcp).increase = sum(sum(turnover_increase(mask == 1)));
        turnover_continents.(continent).(rcp).decrease = sum(sum(turnover_decrease(mask == 1)));
    end
end

% Display results (optional)
disp(turnover_total);
disp(turnover_increase_value);
disp(turnover_decrease_value);
disp(turnover_continents);

%% Creat bar chart
%% all rcps

bar([turnover_continents.Africa.rcp26.increase,turnover_continents.Africa.rcp26.decrease,turnover_continents.Africa.rcp45.increase,turnover_continents.Africa.rcp45.decrease,turnover_continents.Africa.rcp60.increase,turnover_continents.Africa.rcp60.decrease,turnover_continents.Africa.rcp85.increase,turnover_continents.Africa.rcp85.decrease])
bar([turnover_continents.Asia.rcp26.increase,turnover_continents.Asia.rcp26.decrease,turnover_continents.Asia.rcp45.increase,turnover_continents.Asia.rcp45.decrease,turnover_continents.Asia.rcp60.increase,turnover_continents.Asia.rcp60.decrease,turnover_continents.Asia.rcp85.increase,turnover_continents.Asia.rcp85.decrease])
bar([turnover_continents.SouthAmerica.rcp26.increase,turnover_continents.SouthAmerica.rcp26.decrease,turnover_continents.SouthAmerica.rcp45.increase,turnover_continents.SouthAmerica.rcp45.decrease,turnover_continents.SouthAmerica.rcp60.increase,turnover_continents.SouthAmerica.rcp60.decrease,turnover_continents.SouthAmerica.rcp85.increase,turnover_continents.SouthAmerica.rcp85.decrease])
bar([turnover_continents.NorthAmerica.rcp26.increase,turnover_continents.NorthAmerica.rcp26.decrease,turnover_continents.NorthAmerica.rcp45.increase,turnover_continents.NorthAmerica.rcp45.decrease,turnover_continents.NorthAmerica.rcp60.increase,turnover_continents.NorthAmerica.rcp60.decrease,turnover_continents.NorthAmerica.rcp85.increase,turnover_continents.NorthAmerica.rcp85.decrease])
bar([turnover_continents.Europe.rcp26.increase,turnover_continents.Europe.rcp26.decrease,turnover_continents.Europe.rcp45.increase,turnover_continents.Europe.rcp45.decrease,turnover_continents.Europe.rcp60.increase,turnover_continents.Europe.rcp60.decrease,turnover_continents.Europe.rcp85.increase,turnover_continents.Europe.rcp85.decrease])
bar([turnover_continents.Oceania.rcp26.increase,turnover_continents.Oceania.rcp26.decrease,turnover_continents.Oceania.rcp45.increase,turnover_continents.Oceania.rcp45.decrease,turnover_continents.Oceania.rcp60.increase,turnover_continents.Oceania.rcp60.decrease,turnover_continents.Oceania.rcp85.increase,turnover_continents.Oceania.rcp85.decrease])

%% Only rcp85

% Set up the tiled layout for the subplots
figure;
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % 2x3 layout

% Define the colors for gain and loss
colors = [0 1 0; 1 0 0]; % Green for gain, Red for loss

% Plot each continent in a separate subplot
nexttile
bar([turnover_continents.Africa.rcp85.increase, turnover_continents.Africa.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Africa")

nexttile
bar([turnover_continents.Asia.rcp85.increase, turnover_continents.Asia.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Asia")

nexttile
bar([turnover_continents.SouthAmerica.rcp85.increase, turnover_continents.SouthAmerica.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("South America")

nexttile
bar([turnover_continents.NorthAmerica.rcp85.increase, turnover_continents.NorthAmerica.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("North America")

nexttile
bar([turnover_continents.Oceania.rcp85.increase, turnover_continents.Oceania.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Oceania")

nexttile
bar([turnover_continents.Europe.rcp85.increase, turnover_continents.Europe.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Europe")


saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_rcp85_allcontinents.svg');

%% Create maps
%% Custom colormap just for grazing turnover rate (central white, two spectrum red and green.)
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
%% Try one plot first
figure2 = figure;
sgtitle('Turnover', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
%custom_colormap = [1 1 1; addcolorplus(289)];
custom_colormap = [reds_to_white; white_to_greens];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
caxis([-14000 14000]);
c = colorbar;  
R = georefcells([-90,90],[0,360],size(compare_turnover_rcp26));
geoshow(flipud(compare_turnover_rcp85), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])

% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover85.svg');

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
    ax = worldmap('world');
    setm(ax, 'FFaceColor', [1 1 1]);  % Set the map's background color
    
    % Apply the custom colormap
    set(gcf, 'Colormap', custom_colormap);
    
    % Set color axis limits to ensure zero aligns with white
    caxis([-14000 14000]);  % Adjust according to your data's range
    
    % Load and display the data for the current RCP
    compare_turnover = eval(['compare_turnover_' rcp]); % Dynamically load the variable
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
colorbar('Position', [0.93, 0.1, 0.02, 0.8]);  % Custom position for a single colorbar
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_allrcp.svg');

save('./grazingniche/matdata/turnover.mat','compare_turnover_rcp26','compare_turnover_rcp45','compare_turnover_rcp60','compare_turnover_rcp85');     
