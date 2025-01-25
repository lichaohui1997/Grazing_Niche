close all;
%% In this code I am using two new climate models for historical data
%% here are the data: I found 3 usable threads
%% 更多的历史数据气候模式下载
% % This is the thredd: http://esgf-data01.diasjp.net/thredds/catalog/esgcet/2/cmip5.output1.MIROC.MIROC-ESM.past1000.mon.atmos.Amon.r1i1p1.v20120710.xml#cmip5.output1.MIROC.MIROC-ESM.past1000.mon.atmos.Amon.r1i1p1.v20120710|application/xml+thredds|THREDDS
% url='http://esgf-data01.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/past1000/mon/atmos/Amon/r1i1p1/v20120710/hurs/hurs_Amon_MIROC-ESM_past1000_r1i1p1_085001-184912.nc';
% filename='hurs_Amon_MIROC-ESM_past1000_r1i1p1_085001-184912.nc'
% urlwrite(url, filename);
% data = ncread(filename, 'hurs');
% ncdisp(filename)
% imagesc(data(:,:,7))
% colorbar
% %% Another data
% % http://esgf-data01.diasjp.net/thredds/catalog/esgcet/14/cmip5.output1.MRI.MRI-CGCM3.past1000.mon.atmos.Amon.r1i1p1.v20140306.xml#cmip5.output1.MRI.MRI-CGCM3.past1000.mon.atmos.Amon.r1i1p1.v20140306|application/xml+thredds|THREDDS
% url='http://esgf-data01.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-CGCM3/past1000/mon/atmos/Amon/r1i1p1/v20140306/hurs/hurs_Amon_MRI-CGCM3_past1000_r1i1p1_185001-185012.nc'
% filename='hurs_Amon_MRI-CGCM3_past1000_r1i1p1_185001-185012.nc';
% urlwrite(url, filename);
% data = ncread(filename, 'hurs');
% ncdisp(filename)
% imagesc(data(:,:,7))
% colorbar
% %% Another data 和之前下载的非常像，一个是122一个是126.从21到28都有
% %thread: http://esgf.nci.org.au/thredds/catalog/esgcet/cmip5/output1/NASA-GISS/GISS-E2-R/past1000/mon/atmos/Amon/r1i1p122/cmip5.output1.NASA-GISS.GISS-E2-R.past1000.mon.atmos.Amon.r1i1p122.v20160510.xml#cmip5.output1.NASA-GISS.GISS-E2-R.past1000.mon.atmos.Amon.r1i1p122.v20160510|application/xml+thredds|THREDDS
% url='http://esgf.nci.org.au/thredds/fileServer/replica/CMIP5/output1/NASA-GISS/GISS-E2-R/past1000/mon/atmos/Amon/r1i1p122/v20160510/hurs/hurs_Amon_GISS-E2-R_past1000_r1i1p122_180101-185012.nc'
% filename='hurs_Amon_GISS-E2-R_past1000_r1i1p122_175101-180012.nc';
% %urlwrite(url, filename);
% data = ncread(filename, 'hurs');
% ncdisp(filename)
% imagesc(data(:,:,7))
% colorbar

%% Here I am starting the analysis
%% List of variables to download. 
%% This section of code does not need to be run twice
% directory = '/Users/lichaohui/Desktop/grazingniche/rawdata/climatenewnew/MIROC';
% 
% variables = {'hurs', 'pr', 'tas', 'sfcWind'};
% 
% % Base URL
% base_url = 'http://esgf-data01.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/past1000/mon/atmos/Amon/r1i1p1/v20120710/';
% 
% % Loop through each variable to download the file and read it
% for i = 1:length(variables)
%     var = variables{i};
%     url = [base_url var '/' var '_Amon_MIROC-ESM_past1000_r1i1p1_085001-184912.nc'];
%     filename = [var '_Amon_MIROC-ESM_past1000_r1i1p1_085001-184912.nc'];
%     
%     % Download the file
%     %urlwrite(url, filename);
%     
%     % Read the netCDF file
%     data = ncread(filename, var);
%     
%     % Store the data into the workspace with a variable name as defined in 'variables'
%     assignin('base', [var '_data_miroc'], data);
% end
% 
% if ~exist(directory, 'dir')
%     mkdir(directory);
% end
% 
% directory = '/Users/lichaohui/Desktop/grazingniche/climatenewnew/MIROC';
% 
% save(fullfile(directory, 'pr_data_miroc.mat'),'pr_data_miroc');
% save(fullfile(directory, 'sfcWind_data_miroc.mat'),'sfcWind_data_miroc')
% save(fullfile(directory, 'tas_data_miroc.mat'),'tas_data_miroc');
% save(fullfile(directory, 'hurs_data_miroc.mat'),'hurs_data_miroc');

%% Load data and process the raw data, by doing the following:
%% this section also does not need to be run twice
% 1. taking only years 1000,1100...1800; the main output of this section is
% yearly data
% 2. clip over 100 hurs data to 100, turn pr, tas into right value
% this section also does not need to be run twice
load('pr_data_miroc.mat');
load('sfcWind_data_miroc.mat');
load('tas_data_miroc.mat');
load('hurs_data_miroc.mat');

% Years of interest
years_interest = 1000:100:1800;

% Indices for years of interest (subtract 849 because data starts from 850)
start_index = (years_interest - 849) * 12 + 1; % Start from January of each year
end_index = start_index + 11; % Up to December of each year

% Initialize the matrices to store results
tas_yearly_miroc = zeros(128, 64, numel(years_interest));
sfcWind_yearly_miroc = zeros(128, 64, numel(years_interest));
hurs_yearly_miroc = zeros(128, 64, numel(years_interest));
pr_yearly_miroc = zeros(128, 64, numel(years_interest));

% Calculate yearly mean for tas, sfcWind, and hurs
% Calculate yearly sum for pr
for i = 1:numel(years_interest)
    tas_yearly_miroc(:,:,i) = mean(tas_data_miroc(:,:,start_index(i):end_index(i)), 3)-273.15;
    sfcWind_yearly_miroc(:,:,i) = mean(sfcWind_data_miroc(:,:,start_index(i):end_index(i)), 3);
    %hurs_yearly_miroc(:,:,i) = mean(hurs_data_miroc(:,:,start_index(i):end_index(i)), 3);
    capped_hurs_data = min(hurs_data_miroc(:,:,start_index(i):end_index(i)), 100);
    hurs_yearly_miroc(:,:,i) = mean(capped_hurs_data, 3);
    pr_yearly_miroc(:,:,i) = sum(pr_data_miroc(:,:,start_index(i):end_index(i)), 3)*86400*30;
end

directory = '/Users/lichaohui/Desktop/grazingniche/rawdata/climatenewnew/MIROC';
% Save the yearly data
if ~exist(directory, 'dir')
    mkdir(directory);
end

directory = '/Users/lichaohui/Desktop/grazingniche/rawdata/climatenewnew/MIROC';

save(fullfile(directory, 'tas_yearly_miroc.mat'), 'tas_yearly_miroc');
save(fullfile(directory, 'sfcWind_yearly_miroc.mat'), 'sfcWind_yearly_miroc');
save(fullfile(directory, 'hurs_yearly_miroc.mat'), 'hurs_yearly_miroc');
save(fullfile(directory, 'pr_yearly_miroc.mat'), 'pr_yearly_miroc');
%% load the above data. for this section of code, you can start from this section
load('tas_yearly_miroc.mat');
load('sfcWind_yearly_miroc.mat');
load('hurs_yearly_miroc.mat');
load('pr_yearly_miroc.mat');
%% resize and rotate
% rotate
pr_yearly_miroc=rot90(pr_yearly_miroc)
tas_yearly_miroc=rot90(tas_yearly_miroc)
sfcWind_yearly_miroc=rot90(sfcWind_yearly_miroc)
hurs_yearly_miroc=rot90(hurs_yearly_miroc)
% resize
resize_pr_miroc=imresize(pr_yearly_miroc,[90,180],'nearest')
resize_tas_miroc=imresize(tas_yearly_miroc,[90,180],'nearest')
resize_sfcWind_miroc=imresize(sfcWind_yearly_miroc,[90,180],'nearest')
resize_hurs_miroc=imresize(hurs_yearly_miroc,[90,180],'nearest')

%% get a nice graph of it.

figure14 = figure('WindowState','fullscreen');

subplot(2,2,1)
sgtitle('Near surface wind speed', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(342)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_sfcWind_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

subplot(2,2,2)
sgtitle('Mean annual temperature', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(277)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_tas_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

subplot(2,2,3)
sgtitle('Precipitation', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(284)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_pr_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

subplot(2,2,4)
sgtitle('Relative humidity', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_hurs_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/history_climate.svg');

%% another data
% this section of data does not need to be run twice
% I have % out the key line
% directory = '/Users/lichaohui/Desktop/grazingniche/rawdata/climatenewnew/MRI';
% 
% % Base URLs and time periods
% time_periods = {'135001-184912', '085001-134912'};
% base_urls = 'http://esgf-data01.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-CGCM3/past1000/mon/atmos/Amon/r1i1p1/v20140306/';
% 
% % Variables of interest
% variables = {'pr', 'tas', 'hurs', 'sfcWind'};
% 
% % Directory to save data
% save_dir = '/Users/lichaohui/Desktop/calculation/grazingniche/matdata';  % Update this to your specific directory
% 
% % Loop through time periods and variables
% for t = 1:length(time_periods)
%     for v = 1:length(variables)
%         % Construct the full URL and filename
%         url = [base_urls variables{v} '/' variables{v} '_Amon_MRI-CGCM3_past1000_r1i1p1_' time_periods{t} '.nc'];
%         filename = [variables{v} '_Amon_MRI-CGCM3_past1000_r1i1p1_' time_periods{t} '.nc'];
%         
%         % Full path to save the file
%         full_path = fullfile(save_dir, filename);
%         
%         % Download the file
%         fprintf('Downloading %s...\n', url);
%         %urlwrite(url, full_path);
%         
%         % Read the variable into MATLAB
%         fprintf('Reading %s into MATLAB...\n', filename);
%         var_data = ncread(full_path, variables{v});
%         
%         % Save to a .mat file
%         mat_filename = fullfile(save_dir, [variables{v} '_' time_periods{t} '.mat']);
%         save(mat_filename, 'var_data', '-v7.3');
%     end
% end 

%% this code takes some time to run
%% do not need to be run twice
% Define years of interest
% years_of_interest = 1000:100:1800;
% 
% % Initialize empty arrays to hold yearly values
% pr_yearly_mri = zeros(320, 160, numel(years_of_interest));
% hurs_yearly_mri = zeros(320, 160, numel(years_of_interest));
% sfcWind_yearly_mri = zeros(320, 160, numel(years_of_interest));
% tas_yearly_mri = zeros(320, 160, numel(years_of_interest));
% 
% % Loop through each year of interest
% for year = years_of_interest
%     year_idx = (year - 1000)/100 + 1;
%     
%     % Check which dataset to load
%     if year >= 850 && year <= 1349  % data from 085001-134912
%         offset = year - 850;
%         load(['pr_085001-134912.mat']);
%         pr_data_mri = var_data;
%         load(['hurs_085001-134912.mat']);
%         hurs_data_mri = var_data;
%         load(['sfcWind_085001-134912.mat']);
%         sfcWind_data_mri = var_data;
%         load(['tas_085001-134912.mat']);
%         tas_data_mri = var_data;
%         
%     elseif year >= 1350 && year <= 1849  % data from 135001-184912
%         offset = year - 1350;
%         load(['pr_135001-184912.mat']);
%         pr_data_mri = var_data;
%         load(['hurs_135001-184912.mat']);
%         hurs_data_mri = var_data;
%         load(['sfcWind_135001-184912.mat']);
%         sfcWind_data_mri = var_data;
%         load(['tas_135001-184912.mat']);
%         tas_data_mri = var_data;
%     end
%     
%     % Calculate index range for the year in the 3rd dimension
%     idx_start = offset * 12 + 1;  % 12 months per year
%     idx_end = (offset + 1) * 12;
%     
%     % Aggregate/mean calculations
%     pr_yearly_mri(:,:,year_idx) = sum(pr_data_mri(:,:,idx_start:idx_end), 3)*86400*30;
%     %hurs_yearly_mri(:,:,year_idx) = mean(hurs_data_mri(:,:,idx_start:idx_end), 3);
%     % Clip values above 100 to be 100
%     clipped_hurs_data_mri = min(hurs_data_mri, 100);
%     % Then proceed to calculate the mean
%     hurs_yearly_mri(:,:,year_idx) = mean(clipped_hurs_data_mri(:,:,idx_start:idx_end), 3);
% 
%     sfcWind_yearly_mri(:,:,year_idx) = mean(sfcWind_data_mri(:,:,idx_start:idx_end), 3);
%     tas_yearly_mri(:,:,year_idx) = mean(tas_data_mri(:,:,idx_start:idx_end), 3)-273.15;
%     
% end
% 
% 
% % You can now save these aggregated datasets if you want
% 
% directory = '/Users/lichaohui/Desktop/grazingniche/matdata';
% 
% save(fullfile(directory, 'tas_yearly_mri.mat'), 'tas_yearly_mri');
% save(fullfile(directory, 'sfcWind_yearly_mri.mat'), 'sfcWind_yearly_mri');
% save(fullfile(directory, 'hurs_yearly_mri.mat'), 'hurs_yearly_mri');
% save(fullfile(directory, 'pr_yearly_mri.mat'), 'pr_yearly_mri');
%% run this code from here 
load('tas_yearly_mri.mat');
load('sfcWind_yearly_mri.mat');
load('hurs_yearly_mri.mat');
load('pr_yearly_mri.mat');
%% rotate and resize
% rotate
tas_yearly_mri=rot90(tas_yearly_mri);
pr_yearly_mri=rot90(pr_yearly_mri);
sfcWind_yearly_mri=rot90(sfcWind_yearly_mri);
hurs_yearly_mri=rot90(hurs_yearly_mri);
% resize
resize_pr_mri=imresize(pr_yearly_mri,[90,180],'nearest')
resize_tas_mri=imresize(tas_yearly_mri,[90,180],'nearest')
resize_sfcWind_mri=imresize(sfcWind_yearly_mri,[90,180],'nearest')
resize_hurs_mri=imresize(hurs_yearly_mri,[90,180],'nearest')
%% turn them into row vector 
%now my rersizd_hurs_miroc variable is 90*180*8. I want it to be a 12960*1 vector

resize_hurs_miroc_vector = reshape(resize_hurs_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_pr_miroc_vector = reshape(resize_pr_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_tas_miroc_vector = reshape(resize_tas_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_sfcWind_miroc_vector = reshape(resize_sfcWind_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector

resize_hurs_mri_vector = reshape(resize_hurs_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_pr_mri_vector = reshape(resize_pr_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_tas_mri_vector = reshape(resize_tas_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_sfcWind_mri_vector = reshape(resize_sfcWind_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector

%% save these for max's analysis for grassland distribution
resize_hurs_mri_max=resize_hurs_mri_vector;
resize_pr_mri_max=resize_pr_mri_vector;
resize_tas_mri_max=resize_tas_mri_vector;
resize_sfcWind_mri_max=resize_sfcWind_mri_vector;

resize_hurs_miroc_max=resize_hurs_miroc_vector;
resize_pr_miroc_max=resize_pr_miroc_vector;
resize_tas_miroc_max=resize_tas_miroc_vector;
resize_sfcWind_miroc_max=resize_sfcWind_miroc_vector;

resize_hurs_giss_max=resize_hurs_giss_vector
resize_pr_giss_max=resize_pr_giss_vector
resize_tas_giss_max=resize_tas_giss_vector
resize_sfcWind_giss_max=resize_sfcWind_giss_vector

%% get rid of all of the non livestock located land climate points
% here I am doing this for my 3 climate models so far
% Delete the rows where nan is, so that I can produce fits
% becasue distribution fitting does not allow for nan values

resize_hurs_mri_vector(resizedhyde_vector<=0) = []
resize_pr_mri_vector(resizedhyde_vector<=0) = []
resize_sfcWind_mri_vector(resizedhyde_vector<=0) = []
resize_tas_mri_vector(resizedhyde_vector<=0) = []

resize_hurs_miroc_vector(resizedhyde_vector<=0) = []
resize_pr_miroc_vector(resizedhyde_vector<=0) = []
resize_sfcWind_miroc_vector(resizedhyde_vector<=0) = []
resize_tas_miroc_vector(resizedhyde_vector<=0) = []

resize_hurs_giss_vector(resizedhyde_vector<=0) = []
resize_pr_giss_vector(resizedhyde_vector<=0) = []
resize_sfcWind_giss_vector(resizedhyde_vector<=0) = []
resize_tas_giss_vector(resizedhyde_vector<=0) = []

%% calculate the final mean of my three models%%
pr=[resize_pr_miroc_vector,resize_pr_mri_vector,resize_pr_giss_vector]
hurs=[resize_hurs_miroc_vector,resize_hurs_mri_vector,resize_hurs_giss_vector]
sfcWind=[resize_sfcWind_miroc_vector,resize_sfcWind_mri_vector,resize_sfcWind_giss_vector]
tas=[resize_tas_miroc_vector,resize_tas_mri_vector,resize_tas_giss_vector]

pr_max=[resize_pr_miroc_max,resize_pr_mri_max,resize_pr_giss_max]
hurs_max=[resize_hurs_miroc_max,resize_hurs_mri_max,resize_hurs_giss_max]
sfcWind_max=[resize_sfcWind_miroc_max,resize_sfcWind_mri_max,resize_sfcWind_giss_max]
tas_max=[resize_tas_miroc_max,resize_tas_mri_max,resize_tas_giss_max]


