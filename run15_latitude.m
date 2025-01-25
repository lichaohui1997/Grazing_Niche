%% This script calculates the latitude dispersion of grassland for present and future

% load the niche turnover from run8_analysis_future
load('turnover.mat')% this .mat contains compare_turnover_rcp85 and other scenarios variable
load('latitude_modern.mat') %this .mat contains latitude data for modern: variables: latitude_Africa, latitude_Asia, etc.
load('latitude_future.mat')% this .mat file contains future latitude of GN, including variables latitude_Oceania85, latitude_Africa85, etc.
load('longitude_modern.mat') %this .mat contains longitude data for modern: variables: longitude_Africa, longitude_Asia, etc.
load('longitude_future.mat')% this .mat file contains future longitude of GN, including variables longitude_Africa85, etc.

% then run run12 for modern latitude data
%% For loop for all continents drawing turnover_continent specific maps

% Define a struct array with region-specific parameters
regions = struct( ...
    'name', {'Africa', 'SouthAmerica', 'Oceania', 'Asia', 'Europe', 'NorthAmerica'}, ...
    'latLimit', {[-35, 40], [-60, 15], [-50, 10], [5, 55], [35, 70], [15, 85]}, ...
    'lonLimit', {[-20, 55], [-90, -30], [110, 180], [60, 150], [-25, 45], [-170, -50]}, ...
    'data', {'latitude_Africa', 'latitude_SouthAmerica', 'latitude_Oceania', 'latitude_Asia', 'latitude_Europe', 'latitude_NorthAmerica'} ...
);

% Create a single figure with a tiled layout
figure;
tiledlayout(2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % 2x3 grid for six regions

% Loop through each region and create a subplot in each tile
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    latLimit = regions(i).latLimit;
    lonLimit = regions(i).lonLimit;
    latitudeData = eval(regions(i).data); % Load latitude data specific to each region

    % Select the next tile for the current region
    nexttile;
    
    % Custom colormap for display
    custom_colormap = [reds_to_white; white_to_greens];
    ax = axesm('MapProjection', 'miller'); % Miller projection for a rectangular map
    setm(ax, 'FFaceColor', [1 1 1]); % Set background color to white
    set(gcf, 'Colormap', custom_colormap);
    
    % Georeference and display the turnover data for the region
    R = georefcells([-90, 90], [0, 360], size(compare_turnover_rcp85)); % Full global extent
    geoshow(flipud(compare_turnover_rcp85), R, 'DisplayType', 'texturemap'); % Display data
    
    % Set map limits based on the region
    setm(ax, 'MapLatLimit', latLimit, 'MapLonLimit', lonLimit);
    
    % Load and plot coastlines
    load coastlines;
    plotm(coastlat, coastlon, 'Color', 'black');
    
    % Set color axis limits
    caxis([-14000, 14000]);
    
    % Make all axes the same scale
    axis equal;
    
    % Add a title to each subplot
    title(regionName, 'FontSize', 12);
end

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_continents.svg');

%% Continent latitude statistics in subfigure (compare present and future)
% Define a struct array with region-specific parameters
regions = struct( ...
    'name', {'Africa', 'Africa85', 'SouthAmerica', 'SouthAmerica85', 'Oceania', 'Oceania85', ...
             'Asia', 'Asia85', 'Europe', 'Europe85', 'NorthAmerica', 'NorthAmerica85'}, ...
    'latLimit', {[-35, 40], [-35, 40], [-60, 15], [-60, 15], [-50, 10], [-50, 10], ...
                 [5, 55], [5, 55], [35, 70], [35, 70], [15, 85], [15, 85]}, ...
    'lonLimit', {[-20, 55], [-20, 55], [-90, -30], [-90, -30], [110, 180], [110, 180], ...
                 [60, 150], [60, 150], [-25, 45], [-25, 45], [-170, -50], [-170, -50]}, ...
    'data', {'latitude_Africa', 'latitude_Africa85', 'latitude_SouthAmerica', 'latitude_SouthAmerica85', ...
             'latitude_Oceania', 'latitude_Oceania85', 'latitude_Asia', 'latitude_Asia85', ...
             'latitude_Europe', 'latitude_Europe85', 'latitude_NorthAmerica', 'latitude_NorthAmerica85'} ...
);

% Create a single figure with a tiled layout for all subplots
figure;
tiledlayout(3, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % 4x3 grid for twelve regions

% Loop through each region and create a subplot for each
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    latLimit = regions(i).latLimit;
    latitudeData = eval(regions(i).data);  % Load latitude data specific to each region

    % Select the next tile for the current region
    nexttile;
    
    % Generate latitude values for plotting
    ydata = linspace(90, -90, length(latitudeData)); % Latitude values covering the globe
    latitude_norm = latitudeData ./ sum(latitudeData); % Normalize latitude data

    % Plot normalized latitude profile
    plot(latitude_norm, ydata, 'LineWidth', 1.5);

    % Set y-axis limits to match the region's latitude range
    ylim(latLimit);
    %xlabel('Normalized Value');
    ylabel('Latitude');
    xticks([])
    title(sprintf(regionName));
end
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/latitude_all.svg');
%% Continent longitude statistics in subfigure (compare present and future)
% Define a struct array with region-specific parameters
regions = struct( ...
    'name', {'Africa', 'Africa85', 'SouthAmerica', 'SouthAmerica85', 'Oceania', 'Oceania85', ...
             'Asia', 'Asia85', 'Europe', 'Europe85', 'NorthAmerica', 'NorthAmerica85'}, ...
    'latLimit', {[-35, 40], [-35, 40], [-60, 15], [-60, 15], [-50, 10], [-50, 10], ...
                 [5, 55], [5, 55], [35, 70], [35, 70], [15, 85], [15, 85]}, ...
    'lonLimit', {[-20, 55], [-20, 55], [-90, -30], [-90, -30], [110, 180], [110, 180], ...
                 [60, 150], [60, 150], [-25, 45], [-25, 45], [-170, -50], [-170, -50]}, ...    
    'data', {'longitude_Africa', 'longitude_Africa85', 'longitude_SouthAmerica', 'longitude_SouthAmerica85', ...
             'longitude_Oceania', 'longitude_Oceania85', 'longitude_Asia', 'longitude_Asia85', ...
             'longitude_Europe', 'longitude_Europe85', 'longitude_NorthAmerica', 'longitude_NorthAmerica85'} ...
);

% Create a single figure with a tiled layout for all subplots
figure;
tiledlayout(3, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % 4x3 grid for twelve regions

% Loop through each region and create a subplot for each
for i = 1:length(regions)
    % Load specific region parameters
    regionName = regions(i).name;
    lonLimit = regions(i).lonLimit;
    longitudeData = eval(regions(i).data);  % Load longitude data specific to each region

    % Select the next tile for the current region
    nexttile;
    
    % Generate longitude values for plotting
    ydata = linspace(-180, 180, length(longitudeData)); % Latitude values covering the globe
    longitude_norm = longitudeData ./ sum(longitudeData); % Normalize longitude data

    % Plot normalized longitude profile
    plot(ydata,longitude_norm, 'LineWidth', 1.5);

    % Set y-axis limits to match the region's longitude range
    xlim(lonLimit);
    yticks([])
    title(sprintf(regionName));
end

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/longitude_all.svg');
