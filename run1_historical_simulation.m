%% For the code to run smoothly, the first step is to click on the 当前文件夹
% （desktop/calculation/grazingniche)
% 再把command也添加到路径当中
%  note that this actually renders this line of repeated code useless:
% "dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate/'"
%%
% in case you want to access the cmip data in the future but forget how,
% a short tutorial can be found in the pasturecarbon_run.m file. 

%% if you only want to see the niche values and see if you can produce graphs 
%% you only need to run the last section of "run3_historical_climate.m"
%% downloading the historical data from esgf,cmip. 
% url: https://ds.nccs.nasa.gov/thredds/idd/cmip5/cmip5.output1.NASA-GISS.GISS-E2-R.past1000.mon.atmos.Amon.r1i1p126.v20160511.xml#cmip5.output1.NASA-GISS.GISS-E2-R.past1000.mon.atmos.Amon.r1i1p126.v20160511%7Capplication/xml+thredds%7CTHREDDS
% downloading the .nc data into the directory. this section of code 
% mission complete after running once.
% do not run this section of code
% list of variables to download

% 
%http://cmip.bcc.cma.cn/thredds/fileServer
%%
% how to download the data from cmip. Do not use safari. Use chrome. Safari
% cannot open esgf website.
% use the websave command in matlab to directly download the data from
% esgf. The url consists of three parts.
% Use chrome to open this website. https://aims2.llnl.gov/search/cmip6/
% choose from the list. Go to the Metadata-url. The https://... threadds
% will be your first part of the websave url.
% copy paste this entire url into another tab. Open this url. 
% press command+find and find the variables you want to download, such as
% search for "pr"
% near the variable, you can find urlPath. This urlPath will be your last
% part of the url for websave command download. Important: only copy the
% things within the quotation mark, this url should end with .nc, do not
% copy the stuff after the .nc
% websave url=first part+ "/fileServer/"+last part
% websave(filename, url). For filename you can just copy the dataset name.
% There are lots of places to find this dataset name. 
% this is the end result:
% websave('pr_Amon_BCC-CSM2-MR_esm-hist_r1i1p1f1_gn_185001-201412.nc','http://cmip.bcc.cma.cn/thredds/fileServer/cmip6_data/CMIP/BCC/BCC-CSM2-MR/esm-hist/r1i1p1f1/Amon/pr/gn/v20181122/pr_Amon_BCC-CSM2-MR_esm-hist_r1i1p1f1_gn_185001-201412.nc')
% also best to save the raw url so you can re-visit the data when you want
% to. 
% like this: http://cmip.bcc.cma.cn/thredds/catalog/esgcet/2/CMIP6.CMIP.BCC.BCC-CSM2-MR.esm-hist.r1i1p1f1.Amon.pr.gn.v20181122.xml#CMIP6.CMIP.BCC.BCC-CSM2-MR.esm-hist.r1i1p1f1.Amon.pr.gn.v20181122|application/xml+thredds|THREDDS,http://cmip.bcc.cma.cn/las/getUI.do?

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

%% Do not run this section because the above section just loaded the needed variable into the directory

% This section will create a lot of graphs. And also create a lot of data
% named "mean_pr_1.mat" in the directory. The last section justs loads all
% of this data into working space.
base_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate';

% Path to save figures
fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI';

% List of variables to process
variables = {'pr', 'tas', 'sfcWind', 'hurs'};

% Loop over all variables
for var = variables
    % Get a list of all .nc files for the current variable
    files = dir(fullfile(base_path, [var{:} '_*.nc']));
    
    % Loop over all files
    for idx = 1:length(files)
        % Full path to the current file
        filename = fullfile(files(idx).folder, files(idx).name);
        
        % Load the data
        data = ncread(filename, var{:}); 
% in this section of code I am reshaping the data into 4 dimensions, so as
% to aggregate the monthly data into yearly data, then calculate the mean
% of first of every 50 years. so the result is pr, hurs, yearly value in year 1500,1550,1600, etc.

% here I am making the above 100 values of relative humidity 100
        if strcmp(var, 'hurs')
            data(data > 100) = 100;
        end


        % Compute the mean for the first time point and convert to Celsius if temperature data
        if strcmp(var, 'tas') % If the variable is 'tas', convert from Kelvin to Celsius
            data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
            mean_data = rot90(data(:,:,1)) - 273.15;
        elseif strcmp(var, 'pr') % If the variable is 'pr', sum up monthly precipitations
            data = squeeze(sum(reshape(data, 144, 90, 12, []), 3)); % Turn monthly aggregate to yearly aggregate
            mean_data = rot90(data(:,:,1));
            mean_data=86400*30*mean_data
        else
            data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
            mean_data = rot90(data(:,:,1));
        end

        % Store the mean_data to a variable
        eval(['mean_' var{:} '_' num2str(idx) ' = mean_data;']);

        % Plot the data
%         figure;
%         imagesc(mean_data);
%         colorbar;
%         title([var{:} ' File ' num2str(idx)]);

        % Save the figure to a PNG file
        saveas(gcf, fullfile(fig_path, ['mean_' var{:} '_' num2str(idx) '.png']), 'png');

        % Save the variable to a .mat file
        save(['mean_' var{:} '_' num2str(idx) '.mat'], ['mean_' var{:} '_' num2str(idx)]);
    end
end

% Extract the latitude and longitude variables
% lat = ncread(filename, 'lat');
% lon = ncread(filename, 'lon');

% Extract the latitude and longitude bounds
% lat_bnds = ncread(filename, 'lat_bnds');
% lon_bnds = ncread(filename, 'lon_bnds');

% Create the geospatial referencing object
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

%% a nice picture of historical grazing activity
sgtitle('Historical grazing', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resizedhyde{8,1}));
geoshow(flipud(resizedhyde{8,1}), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/history_grazing.svg');

%% now their coordinates are harmonized, put them in a row vector and draw heat map
resizedpr_vector = reshape(cell2mat(resizedpr(:).'), [],1);
resizedtas_vector = reshape(cell2mat(resizedtas(:).'), [],1);
resizedhurs_vector = reshape(cell2mat(resizedhurs(:).'), [],1);
resizedsfcWind_vector = reshape(cell2mat(resizedsfcWind(:).'), [],1);
resizedhyde_vector = reshape(cell2mat(resizedhyde(:).'), [],1);
%% here I only want to take the odd layers, 
% so that the hyde data and climate data are yr 1000, 1100, 1300,...1700. 
% If I don't do this then the resizedpr_vector and resizedhyde_vector won't be the same
resizedpr_vector = reshape(cell2mat(resizedpr(1:2:end).'), [], 1);
resizedtas_vector = reshape(cell2mat(resizedtas(1:2:end).'), [], 1);
resizedhurs_vector = reshape(cell2mat(resizedhurs(1:2:end).'), [], 1);
resizedsfcWind_vector = reshape(cell2mat(resizedsfcWind(1:2:end).'), [], 1);
resizedhyde_vector = reshape(cell2mat(resizedhyde(:).'), [],1);
%% give them standard names so I can juxtapose them with the other climate data
resize_hurs_giss_vector=resizedhurs_vector;
resize_tas_giss_vector=resizedtas_vector;
resize_pr_giss_vector=resizedpr_vector;
resize_sfcWind_giss_vector=resizedsfcWind_vector;
%% save data for max's analysis
resize_hurs_giss_max=resize_hurs_giss_vector
resize_pr_giss_max=resize_pr_giss_vector
resize_tas_giss_max=resize_tas_giss_vector
resize_sfcWind_giss_max=resize_sfcWind_giss_vector

%% get rid of all of the non livestock located land climate points
% Delete the rows where nan is, so that I can produce fits
% becasue distribution fitting does not allow for nan values

% I am doing this now in run3_5historical_climate.m

% resizedpr_vector(resizedhyde_vector<=0) = []
% resizedsfcWind_vector(resizedhyde_vector<=0) = []
% resizedtas_vector(resizedhyde_vector<=0) = []
% resizedhurs_vector(resizedhyde_vector<=0) = []

resizedhyde_vector_s=resizedhyde_vector;
resizedhyde_vector_s(resizedhyde_vector_s<=0) = []

save('./grazingniche/matdata/resizedhyde_vector_s.mat','resizedhyde_vector_s')

%% finding the niche
%% finding the niche
% Filter out data points where resizedhyde_vector < 10
middletas = resizedtas_vector(resizedhyde_vector >= 2);
middlepr = resizedpr_vector(resizedhyde_vector >= 2);
middlesfcWind = resizedsfcWind_vector(resizedhyde_vector >= 2);
middlehurs = resizedhurs_vector(resizedhyde_vector >= 2);
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
%%
subplot(2,3,1)
X=resizedpr_vector;
Y=resizedtas_vector;
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
X=resizedhurs_vector;
Y=resizedtas_vector;
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
X=resizedsfcWind_vector;
Y=resizedtas_vector;
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
%%
nexttile;
hold on;
scatter(resizedpr_vector, resizedhyde_vector,'filled');
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
scatter(resizedtas_vector, resizedhyde_vector,'filled');
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
scatter(resizedhurs_vector, resizedhyde_vector,'filled');
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
scatter(resizedsfcWind_vector, resizedhyde_vector,'filled');
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
% %% Niche map
% 
% % cond1_pr = (aggregate_pr >= 324 & aggregate_pr <= 2548);
% % cond1_tas = (aggregate_tas >= 2 & aggregate_tas <= 28.24);
% % cond1_sfcWind = (aggregate_sfcWind >= 1.13 & aggregate_sfcWind <= 6.50);
% % cond_hurs = (aggregate_hurs >= 34 & aggregate_hurs <= 84);
% % cond1_landuse = (grassland_env>0);
% 
% cond1_pr = (aggregate_pr >= 433 & aggregate_pr <= 2368);
% cond1_tas = (aggregate_tas >= 4 & aggregate_tas <= 28.5);
% cond1_sfcWind = (aggregate_sfcWind >= 1.4 & aggregate_sfcWind <= 6);
% cond_hurs = (aggregate_hurs >= 46 & aggregate_hurs <= 87);
% cond1_landuse = (grassland_env>0);
% 
% % Combine the conditions to create the niche map
% niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;
% 
% % Create the coupled landuse map
% landuse_coupled1 = double(niche1) .* grassland_env;
% 
% % producing the graph of grazing niche
% R = georefcells([-90,90],[-180,180],size(landuse_coupled1));
% 
% figure1 = figure;
% sgtitle('Niche map', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(332)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% c = colorbar;  
% c.Ruler.TickLabelFormat='%g%%';
% geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% set(gcf, 'Position',  [584,449,684,537])
% %% grassland that is outside of the niche
% outniche = (grassland_env > 1) & (landuse_coupled1 == 0); % This creates a logical matrix
% outniche = double(outniche); % Converts logicals to doubles: 1 for true, 0 for false
% outniche=outniche.*grassland_env
% 
% figure9 = figure;
% sgtitle('Grassland outside of the niche', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1; addcolorplus(341)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% c = colorbar;  
% c.Ruler.TickLabelFormat='%g%%';
% R = georefcells([-90,90],[-180,180],size(outniche));
% geoshow(flipud(outniche), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% 
% % set(gcf,'renderer','painters');
% % it seems matlab is unable to produce real svg with this map
% set(gcf, 'Position',  [584,449,684,537])
