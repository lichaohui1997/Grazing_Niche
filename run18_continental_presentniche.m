load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
aggregate_sfcWind(aggregate_sfcWind==66)=0
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
imagesc(aggregate_tas)
%% start the regional niche coding
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
    [resized_data, ~] = georesize(country_data_raw, Rmask, 3600/720, 1800/360, "bilinear");% 把这一行便成需要的目标大小。注意因为原始mask数据的经纬度是反着的矩阵，所以这里的经纬度变换数据也需要反着。
    
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



% now for south america
SouthAmerica = {'m_ARG', 'm_BOL', 'm_BRA', 'm_CHL', 'm_COL', 'm_ECU', 'm_GUY', 'm_PRY', 'm_PER', 'm_SUR', 'm_URY', 'm_VEN'};
SouthAmericamask = country_data.m_ARG;
for i = 2:length(SouthAmerica)
    SouthAmericamask = SouthAmericamask + country_data.(SouthAmerica{i});
end
imagesc(SouthAmericamask)
colorbar


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


% Oceania
Oceania = {'m_AUS', 'm_FJI', 'm_KIR', 'm_FSM', 'm_NCL', 'm_NZL', 'm_NIU', 'm_PLW', 'm_PNG', 'm_WSM', 'm_SLB', 'm_TON', 'm_VUT'};

Oceaniamask = country_data.m_AUS;
for i = 2:length(Oceania)
    Oceaniamask = Oceaniamask + country_data.(Oceania{i});
end


imagesc(Oceaniamask)
colorbar
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

%% make the climate data and livestock data into regional vectors 
%
% check if they are alinged
imagesc(aggregate_pr)
imagesc(Asiamask)
% For the 'pr' variable
resizedpr_Asia = aggregate_pr(Asiamask == 1)
resizedpr_Europe = aggregate_pr(Europemask == 1)
resizedpr_Africa = aggregate_pr(Africamask == 1)
resizedpr_NorthAmerica = aggregate_pr(NorthAmericamask == 1)
resizedpr_SouthAmerica = aggregate_pr(SouthAmericamask == 1)
resizedpr_Oceania = aggregate_pr(Oceaniamask == 1)

% For the 'hurs' variable
resizedhurs_Asia = aggregate_hurs(Asiamask == 1)
resizedhurs_Europe = aggregate_hurs(Europemask == 1)
resizedhurs_Africa = aggregate_hurs(Africamask == 1)
resizedhurs_NorthAmerica = aggregate_hurs(NorthAmericamask == 1)
resizedhurs_SouthAmerica = aggregate_hurs(SouthAmericamask == 1)
resizedhurs_Oceania = aggregate_hurs(Oceaniamask == 1)

% For the 'tas' variable
resizedtas_Asia = aggregate_tas(Asiamask == 1)
resizedtas_Europe = aggregate_tas(Europemask == 1)
resizedtas_Africa = aggregate_tas(Africamask == 1)
resizedtas_NorthAmerica = aggregate_tas(NorthAmericamask == 1)
resizedtas_SouthAmerica = aggregate_tas(SouthAmericamask == 1)
resizedtas_Oceania = aggregate_tas(Oceaniamask == 1)

% For the 'sfcWind' variable
resizedsfcWind_Asia = aggregate_sfcWind(Asiamask == 1)
resizedsfcWind_Europe = aggregate_sfcWind(Europemask == 1)
resizedsfcWind_Africa = aggregate_sfcWind(Africamask == 1)
resizedsfcWind_NorthAmerica = aggregate_sfcWind(NorthAmericamask == 1)
resizedsfcWind_SouthAmerica = aggregate_sfcWind(SouthAmericamask == 1)
resizedsfcWind_Oceania = aggregate_sfcWind(Oceaniamask == 1)

%
% For each continent, using the respective mask to extract cattle data

resizedcattle_Asia = resizecattle(Asiamask == 1);
resizedcattle_Europe = resizecattle(Europemask == 1);
resizedcattle_Africa = resizecattle(Africamask == 1);
resizedcattle_NorthAmerica = resizecattle(NorthAmericamask == 1);
resizedcattle_SouthAmerica = resizecattle(SouthAmericamask == 1);
resizedcattle_Oceania = resizecattle(Oceaniamask == 1);

% For each continent, using the respective mask to extract goat data
resizedgoats_Asia = resizegoats(Asiamask == 1);
resizedgoats_Europe = resizegoats(Europemask == 1);
resizedgoats_Africa = resizegoats(Africamask == 1);
resizedgoats_NorthAmerica = resizegoats(NorthAmericamask == 1);
resizedgoats_SouthAmerica = resizegoats(SouthAmericamask == 1);
resizedgoats_Oceania = resizegoats(Oceaniamask == 1);

% For each continent, using the respective mask to extract sheep data
resizedsheep_Asia = resizesheep(Asiamask == 1);
resizedsheep_Europe = resizesheep(Europemask == 1);
resizedsheep_Africa = resizesheep(Africamask == 1);
resizedsheep_NorthAmerica = resizesheep(NorthAmericamask == 1);
resizedsheep_SouthAmerica = resizesheep(SouthAmericamask == 1);
resizedsheep_Oceania = resizesheep(Oceaniamask == 1);
%% finding the niche for Africa
% Filter out data points where resizedcattle_Africa < 10
middletas = resizedtas_Africa(resizedcattle_Africa >= 2);
middlepr = resizedpr_Africa(resizedcattle_Africa >= 2);
middlesfcWind = resizedsfcWind_Africa(resizedcattle_Africa >= 2);
middlehurs = resizedhurs_Africa(resizedcattle_Africa >= 2);
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
% Filter out data points where resizedcattle_Asia < 10
middletas = resizedtas_Asia(resizedcattle_Asia >= 2);
middlepr = resizedpr_Asia(resizedcattle_Asia >= 2);
middlesfcWind = resizedsfcWind_Asia(resizedcattle_Asia >= 2);
middlehurs = resizedhurs_Asia(resizedcattle_Asia >= 2);
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
% Filter out data points where resizedcattle_NorthAmerica < 10
middletas = resizedtas_NorthAmerica(resizedcattle_NorthAmerica >= 2);
middlepr = resizedpr_NorthAmerica(resizedcattle_NorthAmerica >= 2);
middlesfcWind = resizedsfcWind_NorthAmerica(resizedcattle_NorthAmerica >= 2);
middlehurs = resizedhurs_NorthAmerica(resizedcattle_NorthAmerica >= 2);
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
% Filter out data points where resizedcattle_Europe < 10
middletas = resizedtas_Europe(resizedcattle_Europe >= 2);
middlepr = resizedpr_Europe(resizedcattle_Europe >= 2);
middlesfcWind = resizedsfcWind_Europe(resizedcattle_Europe >= 2);
middlehurs = resizedhurs_Europe(resizedcattle_Europe >= 2);
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
%% finding the niche for SouthAmerica
% Filter out data points where resizedcattle_SouthAmerica < 10
middletas = resizedtas_SouthAmerica(resizedcattle_SouthAmerica >= 2);
middlepr = resizedpr_SouthAmerica(resizedcattle_SouthAmerica >= 2);
middlesfcWind = resizedsfcWind_SouthAmerica(resizedcattle_SouthAmerica >= 2);
middlehurs = resizedhurs_SouthAmerica(resizedcattle_SouthAmerica >= 2);
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
%% finding the niche for Oceania
% Filter out data points where resizedcattle_Oceania < 10
middletas = resizedtas_Oceania(resizedcattle_Oceania >= 2);
middlepr = resizedpr_Oceania(resizedcattle_Oceania >= 2);
middlesfcWind = resizedsfcWind_Oceania(resizedcattle_Oceania >= 2);
middlehurs = resizedhurs_Oceania(resizedcattle_Oceania >= 2);
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
