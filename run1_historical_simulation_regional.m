
variables = {'pr', 'tas', 'hurs', 'sfcWind'};
%variables = {'huss'};

% start and end year
start_year = 05;
end_year = 80;

% base path
base_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/climate';

% create weboptions with a larger timeout
options = weboptions('Timeout', 60);

% loop over the variables
for var = variables
    base_url = ['https://ds.nccs.nasa.gov/thredds/fileServer/CMIP5/NASA/GISS/past1000/E2-R_past1000_r1i1p126/' var{:} '_Amon_GISS-E2-R_past1000_r1i1p126_'];
    base_filename = [var{:} '_Amon_GISS-E2-R_past1000_r1i1p126_'];
    for year = start_year:5:(end_year-5)
        % generate the period
        period = sprintf('1%02d101-1%02d012', year, year+5);

        % construct the url and filename
        url = [base_url period '.nc'];
        filename = [base_filename period '.nc'];

        % full path
        full_path = fullfile(base_path, filename);

        % download the file with increased timeout
        %websave(full_path, url, options);
    end

    % Additional loop for the last file with different naming convention
    final_period = sprintf('1%02d101-1%02d012', start_year, start_year+5);
    url = [base_url final_period '.nc'];
    filename = [base_filename final_period '.nc'];
    full_path = fullfile(base_path, filename);
    %websave(full_path, url, options);
end

%% load the variables into directory, graph them, establish geospatial information
% 

for idx = 1:15
    for var = variables
    filename = ['mean_', var{:},'_',num2str(idx), '.mat'];
    load(filename);
    end
end

% %% Do not run this section because the above section just loaded the needed variable into the directory
% 
% % This section will create a lot of graphs. And also create a lot of data
% % named "mean_pr_1.mat" in the directory. The last section justs loads all
% % of this data into working space.
% base_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate';
% 
% % Path to save figures
% fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI';
% 
% % List of variables to process
% variables = {'pr', 'tas', 'sfcWind', 'hurs'};
% 
% % Loop over all variables
% for var = variables
%     % Get a list of all .nc files for the current variable
%     files = dir(fullfile(base_path, [var{:} '_*.nc']));
%     
%     % Loop over all files
%     for idx = 1:length(files)
%         % Full path to the current file
%         filename = fullfile(files(idx).folder, files(idx).name);
%         
%         % Load the data
%         data = ncread(filename, var{:}); 
% % in this section of code I am reshaping the data into 4 dimensions, so as
% % to aggregate the monthly data into yearly data, then calculate the mean
% % of first of every 50 years. so the result is pr, hurs, yearly value in year 1500,1550,1600, etc.
% 
% % here I am making the above 100 values of relative humidity 100
%         if strcmp(var, 'hurs')
%             data(data > 100) = 100;
%         end
% 
% 
%         % Compute the mean for the first time point and convert to Celsius if temperature data
%         if strcmp(var, 'tas') % If the variable is 'tas', convert from Kelvin to Celsius
%             data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
%             mean_data = rot90(data(:,:,1)) - 273.15;
%         elseif strcmp(var, 'pr') % If the variable is 'pr', sum up monthly precipitations
%             data = squeeze(sum(reshape(data, 144, 90, 12, []), 3)); % Turn monthly aggregate to yearly aggregate
%             mean_data = rot90(data(:,:,1));
%             mean_data=86400*30*mean_data
%         else
%             data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
%             mean_data = rot90(data(:,:,1));
%         end
% 
%         % Store the mean_data to a variable
%         eval(['mean_' var{:} '_' num2str(idx) ' = mean_data;']);
% 
%         % Plot the data
% %         figure;
% %         imagesc(mean_data);
% %         colorbar;
% %         title([var{:} ' File ' num2str(idx)]);
% 
%         % Save the figure to a PNG file
%         saveas(gcf, fullfile(fig_path, ['mean_' var{:} '_' num2str(idx) '.png']), 'png');
% 
%         % Save the variable to a .mat file
%         save(['mean_' var{:} '_' num2str(idx) '.mat'], ['mean_' var{:} '_' num2str(idx)]);
%     end
% end
% 
% % Extract the latitude and longitude variables
% lat = ncread(filename, 'lat');
% lon = ncread(filename, 'lon');
% 
% % Extract the latitude and longitude bounds
% lat_bnds = ncread(filename, 'lat_bnds');
% lon_bnds = ncread(filename, 'lon_bnds');
% 
% % Create the geospatial referencing object
% R_history = georasterref();
% R_history.RasterSize = [size(data, 2), size(data, 1)];
% R_history.LatitudeLimits = [min(lat_bnds(:)), max(lat_bnds(:))];
% R_history.LongitudeLimits = [min(lon_bnds(:)), max(lon_bnds(:))];
% R_history.CellExtentInLatitude = abs(mean(diff(lat)));
% R_history.CellExtentInLongitude = abs(mean(diff(lon)));
% 
% % Display the geospatial referencing object
% disp(R_history);

%% load the ancient grazing data
% I think this bit is a bit tricky because the total number of asc files cell size
% always doesn't seem to equal the HeaderInfo nclo*nrow. I solved this
% temporarily by manipulating x = xllcorner + (0:nrows-1) * cellsize; % x-coordinates
% where it's supposed to be x = xllcorner + (0:ncols-1) * cellsize; %
% x-coordinates, otherwise the established coordinates will not be
% consistent with gridMatrix

%% load the downloaded ancient grazing data into directory and graph the figures
% Base path where the data files are stored
% you will need to run this section
% This section of code will take around 1 min to run

% Path to save figures
fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI';

% Create an array of years
years = 1000:100:1800;

% Loop through the array of years
for i = 1:length(years)
    % Load the data file for the current year
    filename = ['pasture', num2str(years(i)), 'AD.asc'];
    data = importdata(filename);
    
    % Extract the grid values and header information
    gridData = data.data;
    headerInfo = data.textdata;
    
    % Extract header information
    ncols = gridData(1);
    nrows = gridData(2);
    xllcorner = gridData(3);
    yllcorner = gridData(4);
    cellsize = gridData(5);
    
    % Create coordinate vectors
    x = xllcorner + (0:nrows-1) * cellsize;
    y = yllcorner + (0:nrows-1) * cellsize;
    
    % Reshape the grid data into a 2D matrix
    gridMatrix = reshape(gridData(7:size(gridData)), [], nrows)';

    % Display the image
    gridMatrix(gridMatrix < 0) = 0;
    imagesc(gridMatrix);
    colormap([1, 1, 1; jet(255)]); % White color at index 1, followed by the jet colormap
    colorbar;
    
    % Add a title
    title(['Historical Grazing Distribution: Year ', num2str(years(i))]);

    % Save the figure to a PNG file
    saveas(gcf, fullfile(fig_path, ['historical_grazing_distribution_' num2str(years(i)) '.png']), 'png');
    
    % Save the gridMatrix to a variable named gridMatrix and the current year
    assignin('base', ['gridMatrix', num2str(years(i))], gridMatrix);
end

R_hyde = georasterref('LatitudeLimits', [min(y), max(y)], 'LongitudeLimits', [min(x), max(x)], ...
        'RasterSize', [nrows, ncols], 'RasterInterpretation', 'cells');

disp(R_hyde);


%% changing the coordinates of the hyde data from -180-180 to 0-360, 
%% so that the layout of the world map is consistent with the historical climatology maps
% I did this in a smart way (or a dumb way) by swaping the left and right side of the map

% Create an array of years
years = 1000:100:1800;

% Loop through the array of years
for i = 1:length(years)
    % Get the matrix for the current year
    gridMatrix = eval(['gridMatrix', num2str(years(i))]);
    
    % Calculate the half size
    halfSize = size(gridMatrix, 2) / 2;
    
    % Shift the matrix
    shiftedMatrix = circshift(gridMatrix, [0, halfSize]);
    
    % Save the shifted matrix to a variable named shiftedMatrix and the current year
    assignin('base', ['shiftedMatrix', num2str(years(i))], shiftedMatrix);
end


imagesc(shiftedMatrix1800);
colormap([1, 1, 1; jet(255)]); % White color at index 1, followed by the jet colormap
colorbar;  
% alright, the graph shows the conversion is successful.



R_hyde = georasterref('LatitudeLimits', [min(y), max(y)], 'LongitudeLimits', [0,360], ...
        'RasterSize', [nrows, ncols], 'RasterInterpretation', 'cells');

disp(R_hyde)

%% Rescale all of the data into 1degree 
numfile=15
resizedpr = cell(numfile,1);
resizedtas = cell(numfile,1);
resizedhurs = cell(numfile,1);
resizedsfcWind = cell(numfile,1);
resizedhyde = cell(numfile,1);

for i = 1:15
    resizedpr{i} = imresize(eval(['mean_pr_', num2str(i)]), [90 180], 'nearest');
    resizedtas{i} = imresize(eval(['mean_tas_', num2str(i)]), [90 180], 'nearest');
    resizedhurs{i} = imresize(eval(['mean_hurs_', num2str(i)]), [90 180], 'nearest');
    resizedsfcWind{i} = imresize(eval(['mean_sfcWind_', num2str(i)]), [90 180], 'nearest');

end

for i = 1:8
    year = 1000 + (i-1) * 100;
    resizedhyde{i} = imresize(eval(['shiftedMatrix', num2str(year)]), [90, 180]);
end

% %% a nice picture of historical grazing activity
% sgtitle('Historical grazing', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(309)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% R = georefcells([-90,90],[0,360],size(resizedhyde{8,1}));
% geoshow(flipud(resizedhyde{8,1}), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/history_grazing.svg');
%% the above code are the same. now start the regional niche coding
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
Rmask = georefcells([-90,90],[-180,180], size(ncread(file_path, 'm_world')));%这一行不需要变。只需要变下面的144/720，将144便成所需要的即可

% Loop through countries, read and resize each mask
for i = 1:length(countries)
    country_name = countries{i};
    
    % Read country data
    country_data_raw = ncread(file_path, country_name);
    
    % Resize the country data
    [resized_data, ~] = georesize(country_data_raw, Rmask, 180/720, 90/360, "bilinear");% 把这一行便成需要的目标大小。注意因为原始mask数据的经纬度是反着的矩阵，所以这里的经纬度变换数据也需要反着。
    
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
% africa
Africa = {'m_DZA', 'm_COD','m_SDN', 'm_AGO', 'm_BEN', 'm_BWA', 'm_BFA', 'm_BDI', 'm_CMR', 'm_CPV', 'm_CAF', 'm_TCD', 'm_COM', 'm_COG', 'm_CIV', 'm_DJI', 'm_EGY', 'm_GNQ', 'm_ERI', 'm_ETH', 'm_GAB', 'm_GMB', 'm_GHA', 'm_GIN', 'm_GNB', 'm_KEN', 'm_LSO', 'm_LBR', 'm_LBY', 'm_MDG', 'm_MWI', 'm_MLI', 'm_MRT', 'm_MUS', 'm_MAR', 'm_MOZ', 'm_NAM', 'm_NER', 'm_NGA', 'm_REU', 'm_RWA', 'm_STP', 'm_SEN', 'm_SLE', 'm_SOM', 'm_ZAF', 'm_SSD', 'm_ESH', 'm_SWZ', 'm_TZA', 'm_TGO', 'm_TUN', 'm_UGA', 'm_ZMB', 'm_ZWE'};
Africamask = country_data.m_DZA;
for i = 2:length(Africa)
    Africamask = Africamask + country_data.(Africa{i});
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

Asiamask(Asiamask > 1) = 1;
Europemask(Europemask > 1) = 1;
Oceaniamask(Oceaniamask > 1) = 1;
Africamask(Africamask > 1) = 1;
NorthAmericamask(NorthAmericamask > 1) = 1;
SouthAmericamask(SouthAmericamask > 1) = 1;

%% make the climate data and hyde data into regional vectors 
%
resizedpr_Asia = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedpr, 'UniformOutput', false));
resizedpr_Europe = cell2mat(cellfun(@(x) x(Europemask == 1), resizedpr, 'UniformOutput', false));
resizedpr_Africa = cell2mat(cellfun(@(x) x(Africamask == 1), resizedpr, 'UniformOutput', false));
resizedpr_Oceania = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedpr, 'UniformOutput', false));
resizedpr_NorthAmerica = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedpr, 'UniformOutput', false));
resizedpr_SouthAmerica = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedpr, 'UniformOutput', false));
%
resizedhurs_Asia = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedhurs, 'UniformOutput', false));
resizedhurs_Europe = cell2mat(cellfun(@(x) x(Europemask == 1), resizedhurs, 'UniformOutput', false));
resizedhurs_Africa = cell2mat(cellfun(@(x) x(Africamask == 1), resizedhurs, 'UniformOutput', false));
resizedhurs_Oceania = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedhurs, 'UniformOutput', false));
resizedhurs_NorthAmerica = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedhurs, 'UniformOutput', false));
resizedhurs_SouthAmerica = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedhurs, 'UniformOutput', false));
%
resizedtas_Asia = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedtas, 'UniformOutput', false));
resizedtas_Europe = cell2mat(cellfun(@(x) x(Europemask == 1), resizedtas, 'UniformOutput', false));
resizedtas_Africa = cell2mat(cellfun(@(x) x(Africamask == 1), resizedtas, 'UniformOutput', false));
resizedtas_Oceania = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedtas, 'UniformOutput', false));
resizedtas_NorthAmerica = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedtas, 'UniformOutput', false));
resizedtas_SouthAmerica = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedtas, 'UniformOutput', false));
%
resizedsfcWind_Asia = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedsfcWind, 'UniformOutput', false));
resizedsfcWind_Europe = cell2mat(cellfun(@(x) x(Europemask == 1), resizedsfcWind, 'UniformOutput', false));
resizedsfcWind_Africa = cell2mat(cellfun(@(x) x(Africamask == 1), resizedsfcWind, 'UniformOutput', false));
resizedsfcWind_Oceania = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedsfcWind, 'UniformOutput', false));
resizedsfcWind_NorthAmerica = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedsfcWind, 'UniformOutput', false));
resizedsfcWind_SouthAmerica = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedsfcWind, 'UniformOutput', false));

%
resizedhyde_Asia = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedhyde(1:8), 'UniformOutput', false));% becasue hyde does not have data after 9. only 1000-1800
resizedhyde_Europe = cell2mat(cellfun(@(x) x(Europemask == 1), resizedhyde(1:8), 'UniformOutput', false));
resizedhyde_Oceania = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedhyde(1:8), 'UniformOutput', false));
resizedhyde_Africa = cell2mat(cellfun(@(x) x(Africamask == 1), resizedhyde(1:8), 'UniformOutput', false));
resizedhyde_SouthAmerica = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedhyde(1:8), 'UniformOutput', false));
resizedhyde_NorthAmerica = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedhyde(1:8), 'UniformOutput', false));

%% here I only want to take the odd layers, 
% so that the hyde data and climate data are yr 1000, 1100, 1300,...1700. 
% If I don't do this then the resizedpr_vector and resizedhyde_vector won't be the same
%
resizedpr_Asia_odd = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedpr(1:2:end), 'UniformOutput', false));
resizedpr_Africa_odd = cell2mat(cellfun(@(x) x(Africamask == 1), resizedpr(1:2:end), 'UniformOutput', false));
resizedpr_Europe_odd = cell2mat(cellfun(@(x) x(Europemask == 1), resizedpr(1:2:end), 'UniformOutput', false));
resizedpr_Oceania_odd = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedpr(1:2:end), 'UniformOutput', false));
resizedpr_NorthAmerica_odd = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedpr(1:2:end), 'UniformOutput', false));
resizedpr_SouthAmerica_odd = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedpr(1:2:end), 'UniformOutput', false));
%
resizedsfcWind_Asia_odd = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedsfcWind(1:2:end), 'UniformOutput', false));
resizedsfcWind_Africa_odd = cell2mat(cellfun(@(x) x(Africamask == 1), resizedsfcWind(1:2:end), 'UniformOutput', false));
resizedsfcWind_Europe_odd = cell2mat(cellfun(@(x) x(Europemask == 1), resizedsfcWind(1:2:end), 'UniformOutput', false));
resizedsfcWind_Oceania_odd = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedsfcWind(1:2:end), 'UniformOutput', false));
resizedsfcWind_NorthAmerica_odd = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedsfcWind(1:2:end), 'UniformOutput', false));
resizedsfcWind_SouthAmerica_odd = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedsfcWind(1:2:end), 'UniformOutput', false));
%
resizedtas_Asia_odd = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedtas(1:2:end), 'UniformOutput', false));
resizedtas_Africa_odd = cell2mat(cellfun(@(x) x(Africamask == 1), resizedtas(1:2:end), 'UniformOutput', false));
resizedtas_Europe_odd = cell2mat(cellfun(@(x) x(Europemask == 1), resizedtas(1:2:end), 'UniformOutput', false));
resizedtas_Oceania_odd = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedtas(1:2:end), 'UniformOutput', false));
resizedtas_NorthAmerica_odd = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedtas(1:2:end), 'UniformOutput', false));
resizedtas_SouthAmerica_odd = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedtas(1:2:end), 'UniformOutput', false));
%
resizedhurs_Asia_odd = cell2mat(cellfun(@(x) x(Asiamask == 1), resizedhurs(1:2:end), 'UniformOutput', false));
resizedhurs_Africa_odd = cell2mat(cellfun(@(x) x(Africamask == 1), resizedhurs(1:2:end), 'UniformOutput', false));
resizedhurs_Europe_odd = cell2mat(cellfun(@(x) x(Europemask == 1), resizedhurs(1:2:end), 'UniformOutput', false));
resizedhurs_Oceania_odd = cell2mat(cellfun(@(x) x(Oceaniamask == 1), resizedhurs(1:2:end), 'UniformOutput', false));
resizedhurs_NorthAmerica_odd = cell2mat(cellfun(@(x) x(NorthAmericamask == 1), resizedhurs(1:2:end), 'UniformOutput', false));
resizedhurs_SouthAmerica_odd = cell2mat(cellfun(@(x) x(SouthAmericamask == 1), resizedhurs(1:2:end), 'UniformOutput', false));

%% finding the niche for Africa
% Filter out data points where resizedhyde_Africa < 10
middletas = resizedtas_Africa_odd(resizedhyde_Africa >= 2);
middlepr = resizedpr_Africa_odd(resizedhyde_Africa >= 2);
middlesfcWind = resizedsfcWind_Africa_odd(resizedhyde_Africa >= 2);
middlehurs = resizedhurs_Africa_odd(resizedhyde_Africa >= 2);
%
lower_bound = prctile(middletas, 2.75);
upper_bound = prctile(middletas, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlepr, 2.75);
upper_bound = prctile(middlepr, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlesfcWind, 2.75);
upper_bound = prctile(middlesfcWind, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlehurs, 2.75);
upper_bound = prctile(middlehurs, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);

%% finding the niche for Asia
% Filter out data points where resizedhyde_Asia < 10
middletas = resizedtas_Asia_odd(resizedhyde_Asia >= 2);
middlepr = resizedpr_Asia_odd(resizedhyde_Asia >= 2);
middlesfcWind = resizedsfcWind_Asia_odd(resizedhyde_Asia >= 2);
middlehurs = resizedhurs_Asia_odd(resizedhyde_Asia >= 2);
%
lower_bound = prctile(middletas, 2.75);
upper_bound = prctile(middletas, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlepr, 2.75);
upper_bound = prctile(middlepr, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlesfcWind, 2.75);
upper_bound = prctile(middlesfcWind, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlehurs, 2.75);
upper_bound = prctile(middlehurs, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%% finding the niche for NorthAmerica
% Filter out data points where resizedhyde_NorthAmerica < 10
middletas = resizedtas_NorthAmerica_odd(resizedhyde_NorthAmerica >= 2);
middlepr = resizedpr_NorthAmerica_odd(resizedhyde_NorthAmerica >= 2);
middlesfcWind = resizedsfcWind_NorthAmerica_odd(resizedhyde_NorthAmerica >= 2);
middlehurs = resizedhurs_NorthAmerica_odd(resizedhyde_NorthAmerica >= 2);
%
lower_bound = prctile(middletas, 2.75);
upper_bound = prctile(middletas, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlepr, 2.75);
upper_bound = prctile(middlepr, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlesfcWind, 2.75);
upper_bound = prctile(middlesfcWind, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlehurs, 2.75);
upper_bound = prctile(middlehurs, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%% finding the niche for Europe
% Filter out data points where resizedhyde_Europe < 10
middletas = resizedtas_Europe_odd(resizedhyde_Europe >= 2);
middlepr = resizedpr_Europe_odd(resizedhyde_Europe >= 2);
middlesfcWind = resizedsfcWind_Europe_odd(resizedhyde_Europe >= 2);
middlehurs = resizedhurs_Europe_odd(resizedhyde_Europe >= 2);
%
lower_bound = prctile(middletas, 2.75);
upper_bound = prctile(middletas, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlepr, 2.75);
upper_bound = prctile(middlepr, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlesfcWind, 2.75);
upper_bound = prctile(middlesfcWind, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlehurs, 2.75);
upper_bound = prctile(middlehurs, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
