addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% import data
% in order for this .m file to work, you would need to import:
load('grassland_env.mat')
load('livestockdensity.mat')
load('AGB.mat')
load('niche.mat')
load('futureniche.mat')
load('futureland.mat')
load('livestock2015.mat')

imagesc(grassland_env)
imagesc(futureland_rcp26)
imagesc(futureniche_struct.rcp26)
% 土地的数据和气候的数据之间需要进行坐标转换，气候的数据一律是非洲在左边（0-360），土地的数据是非洲在中间（-180-180）
% 气候的数据太多不好调整，所以一律把土地的数据调整成和气候数据一致，即非洲在最左边（0-360）。
grassland_env = grassland_env(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp26 = futureland_rcp26(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp45 = futureland_rcp45(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp60 = futureland_rcp60(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp85 = futureland_rcp85(:, [ceil(end/2+1):end, 1:floor(end/2)]);

%% how will the figures change in the future? 
% will the situation be still the same? I am expecting that the grazing has
% moved to europe, and the niche area in asia will greatly decrease. This
% is because the precipitation in asia has moved outside of the niche
% range. I will need this analysis in the main text, and figure in SI. 

%% first present a world area pixel map that corresponds to the 2100 scenario
%% pixel to area calculation [function]
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


%% land surface area
landsurfacearea=sum(sum(worldarea(country_data.m_world==1)));
Asiaarea=sum(sum(worldarea(Asiamask==1)));
Africaarea=sum(sum(worldarea(Africamask==1)));
SouthAmericAarea=sum(sum(worldarea(SouthAmericamask==1)));
NorthAmericaArea=sum(sum(worldarea(NorthAmericamask==1)));
Oceaniaarea=sum(sum(worldarea(Oceaniamask==1)));
Europearea=sum(sum(worldarea(Europemask==1)));

%% grassland area
% this section of code take 15s to run
% how much grassland are there in the world?
grasslandarea_rcp26=futureland_rcp26.*worldarea
grasslandarea_world_rcp26=sum(sum(grasslandarea_rcp26))

grasslandarea_rcp85=futureland_rcp85.*worldarea
grasslandarea_world_rcp85=sum(sum(grasslandarea_rcp85))

grasslandarea_rcp45 = futureland_rcp45 .* worldarea;
grasslandarea_world_rcp45 = sum(sum(grasslandarea_rcp45));

grasslandarea_rcp60 = futureland_rcp60 .* worldarea;
grasslandarea_world_rcp60 = sum(sum(grasslandarea_rcp60));

imagesc(grasslandarea_rcp26)

% 2.8523e+07
% 未来rcp26值是3.27+07(LUH1)
% 换了一个数据算，结果变成了1.84e+07(under rcp2.6)(gfdl-esm)
% 又换了一个数据，变成了2.4+07(LUH2)-2.8+07

%how much grassland are there in each continent?
grasslandarea_Asia_rcp26=sum(sum(grasslandarea_rcp26(Asiamask==1)))
grasslandarea_Africa_rcp26=sum(sum(grasslandarea_rcp26(Africamask==1)))
grasslandarea_Oceania_rcp26=sum(sum(grasslandarea_rcp26(Oceaniamask==1)))
grasslandarea_SouthAmerica_rcp26=sum(sum(grasslandarea_rcp26(SouthAmericamask==1)))
grasslandarea_NorthAmerica_rcp26=sum(sum(grasslandarea_rcp26(NorthAmericamask==1)))
grasslandarea_Europe_rcp26=sum(sum(grasslandarea_rcp26(Europemask==1)))

grasslandarea_Asia_rcp85=sum(sum(grasslandarea_rcp85(Asiamask==1)))
grasslandarea_Africa_rcp85=sum(sum(grasslandarea_rcp85(Africamask==1)))
grasslandarea_Oceania_rcp85=sum(sum(grasslandarea_rcp85(Oceaniamask==1)))
grasslandarea_SouthAmerica_rcp85=sum(sum(grasslandarea_rcp85(SouthAmericamask==1)))
grasslandarea_NorthAmerica_rcp85=sum(sum(grasslandarea_rcp85(NorthAmericamask==1)))
grasslandarea_Europe_rcp85=sum(sum(grasslandarea_rcp85(Europemask==1)))

grasslandarea_Asia_rcp45 = sum(sum(grasslandarea_rcp45(Asiamask == 1)));
grasslandarea_Africa_rcp45 = sum(sum(grasslandarea_rcp45(Africamask == 1)));
grasslandarea_Oceania_rcp45 = sum(sum(grasslandarea_rcp45(Oceaniamask == 1)));
grasslandarea_SouthAmerica_rcp45 = sum(sum(grasslandarea_rcp45(SouthAmericamask == 1)));
grasslandarea_NorthAmerica_rcp45 = sum(sum(grasslandarea_rcp45(NorthAmericamask == 1)));
grasslandarea_Europe_rcp45 = sum(sum(grasslandarea_rcp45(Europemask == 1)));

grasslandarea_Asia_rcp60 = sum(sum(grasslandarea_rcp60(Asiamask == 1)));
grasslandarea_Africa_rcp60 = sum(sum(grasslandarea_rcp60(Africamask == 1)));
grasslandarea_Oceania_rcp60 = sum(sum(grasslandarea_rcp60(Oceaniamask == 1)));
grasslandarea_SouthAmerica_rcp60 = sum(sum(grasslandarea_rcp60(SouthAmericamask == 1)));
grasslandarea_NorthAmerica_rcp60 = sum(sum(grasslandarea_rcp60(NorthAmericamask == 1)));
grasslandarea_Europe_rcp60 = sum(sum(grasslandarea_rcp60(Europemask == 1)));

% how much percent does that take?
grasslandprc_rcp26=grasslandarea_world_rcp26/landsurfacearea;


% In 2100, under rcp2.6 scenario, there will be 22081754 km2 of grassland 
% and 23797086 km2 of grassland under rcp8.5 scenario. This value does not
% deviate far from the current situation. 

%% how much "usable" grassland area are there in the world?
% for understanding: landuse_coupled1 basically is a map of usable
% grassland for raising livestock. The values in each grid is the percent
% of land area in grassland. 

nicheland_world_rcp26=futureniche_struct.rcp26.*worldarea;
nichelandarea_world_rcp26=sum(sum(nicheland_world_rcp26));

imagesc(futureniche_struct.rcp26)
colorbar

nicheland_world_rcp85=futureniche_struct.rcp85.*worldarea;
nichelandarea_world_rcp85=sum(sum(nicheland_world_rcp85));

nicheland_world_rcp45=futureniche_struct.rcp45.*worldarea;
nichelandarea_world_rcp45=sum(sum(nicheland_world_rcp45));

nicheland_world_rcp60=futureniche_struct.rcp60.*worldarea;
nichelandarea_world_rcp60=sum(sum(nicheland_world_rcp60));
% 1.08*10^7
% 用LUH1数据计算是反而未来还多了
% rcp26: 1.73+07
% rcp85: 1.42+07
% 这是一个很不好的结果：未来不仅总草地的面积会增加，而且未来可用的草地的面积也会增加
% 这如何汇报？
% 我觉得这和未来草地的面积预测很有关系。不同的资料草地面积差别实在是太大了。

% 现在用LUH2计算，面积变成了7.3*10^6(rcp26)-7.8*10^6(rcp85)
% rcp85反而比rcp26的面积还多，因为青藏高原上可利用的面积变多了。

% However, the usable grassland will significantly decline in the future,
% with only 7380705 km2 of grassland still inside the climatic zone and 7.8609e+06 km2
% under the rcp8.5 scenario. We can see from the map that usable grassland
% has significantly reduced in all continents, and have relocated to
% several concentrated regions in the world, most notably the Tibetan
% Plateau. The Tibetan Plateau was outside of the niche in the 2015
% analysis due to its temperature was outside of the niche, but in the year
% 2100, it is the most important source of niche grassland. Other regions
% still within the niche is Australia and the northen Great Plains in
% America, and a bit on Pampas. There will be no usable grassland left in
% Africa in 2100 under both scenarios. 
% The dispersion of usable grassland in 2100 will be significantly
% different from that of the current situation. Usable grassland area in
% Oceania will reduce but the reduction is only slight, from 1.6+10^6 km2 
% to 1.3+10^6 km2 (rcp2.6) and 1.5 km2 (rcp8.5), usable grassland area in 
% North America will greatly decrease,from 2.6+10^6 to 7.5+10^5 km2
% (rcp2.6) and 8.0+10^5 km2 (rcp8.5), respectively. The greatest decrease
% will be that of Europe and Africa, almost all of its usable grassland
% will be outside the climate niche in 2100. As for Asia, its usable
% grassland area will in fact will more than double in 2100, from 2.3+10^6 km2 to
% 5.0+10^6 km2(rcp2.6) and 5.3+10^6 km2 (rcp8.5). This is due to that due
% to temperature rise, the Tibetan Plateau will then be within the suitable
% grazing temperature zones. 
% 
% respectively. Usable grassland will decrease from 

% how big is the area of usable grassland are there in each continent?
nichelandarea_Asia_rcp26=sum(sum(nicheland_world_rcp26(Asiamask==1)));
nichelandarea_Africa_rcp26=sum(sum(nicheland_world_rcp26(Africamask==1)));
nichelandarea_SouthAmerica_rcp26=sum(sum(nicheland_world_rcp26(SouthAmericamask==1)));
nichelandarea_NorthAmerica_rcp26=sum(sum(nicheland_world_rcp26(NorthAmericamask==1)));
nichelandarea_Oceania_rcp26=sum(sum(nicheland_world_rcp26(Oceaniamask==1)));
nichelandarea_Europe_rcp26=sum(sum(nicheland_world_rcp26(Europemask==1)));

nichelandarea_Asia_rcp85=sum(sum(nicheland_world_rcp85(Asiamask==1)));
nichelandarea_Africa_rcp85=sum(sum(nicheland_world_rcp85(Africamask==1)));
nichelandarea_SouthAmerica_rcp85=sum(sum(nicheland_world_rcp85(SouthAmericamask==1)));
nichelandarea_NorthAmerica_rcp85=sum(sum(nicheland_world_rcp85(NorthAmericamask==1)));
nichelandarea_Oceania_rcp85=sum(sum(nicheland_world_rcp85(Oceaniamask==1)));
nichelandarea_Europe_rcp85=sum(sum(nicheland_world_rcp85(Europemask==1)));

nichelandarea_Asia_rcp45 = sum(sum(nicheland_world_rcp45(Asiamask == 1)));
nichelandarea_Africa_rcp45 = sum(sum(nicheland_world_rcp45(Africamask == 1)));
nichelandarea_SouthAmerica_rcp45 = sum(sum(nicheland_world_rcp45(SouthAmericamask == 1)));
nichelandarea_NorthAmerica_rcp45 = sum(sum(nicheland_world_rcp45(NorthAmericamask == 1)));
nichelandarea_Oceania_rcp45 = sum(sum(nicheland_world_rcp45(Oceaniamask == 1)));
nichelandarea_Europe_rcp45 = sum(sum(nicheland_world_rcp45(Europemask == 1)));

nichelandarea_Asia_rcp60 = sum(sum(nicheland_world_rcp60(Asiamask == 1)));
nichelandarea_Africa_rcp60 = sum(sum(nicheland_world_rcp60(Africamask == 1)));
nichelandarea_SouthAmerica_rcp60 = sum(sum(nicheland_world_rcp60(SouthAmericamask == 1)));
nichelandarea_NorthAmerica_rcp60 = sum(sum(nicheland_world_rcp60(NorthAmericamask == 1)));
nichelandarea_Oceania_rcp60 = sum(sum(nicheland_world_rcp60(Oceaniamask == 1)));
nichelandarea_Europe_rcp60 = sum(sum(nicheland_world_rcp60(Europemask == 1)));

% what percentage grassland are usable in each continent?

nicheprc_Asia_rcp26=nichelandarea_Asia_rcp26/grasslandarea_Asia_rcp26;
nicheprc_Africa_rcp26=nichelandarea_Africa_rcp26/grasslandarea_Africa_rcp26;
nicheprc_SouthAmerica_rcp26=nichelandarea_SouthAmerica_rcp26/grasslandarea_SouthAmerica_rcp26;
nicheprc_NorthAmerica_rcp26=nichelandarea_NorthAmerica_rcp26/grasslandarea_NorthAmerica_rcp26;
nicheprc_Oceania_rcp26=nichelandarea_Oceania_rcp26/grasslandarea_Oceania_rcp26;
nicheprc_Europe_rcp26=nichelandarea_Europe_rcp26/grasslandarea_Europe_rcp26;
nicheprc_world_rcp26=nichelandarea_world_rcp26./grasslandarea_world_rcp26;

nicheprc_Asia_rcp85=nichelandarea_Asia_rcp85/grasslandarea_Asia_rcp85;
nicheprc_Africa_rcp85=nichelandarea_Africa_rcp85/grasslandarea_Africa_rcp85;
nicheprc_SouthAmerica_rcp85=nichelandarea_SouthAmerica_rcp85/grasslandarea_SouthAmerica_rcp85;
nicheprc_NorthAmerica_rcp85=nichelandarea_NorthAmerica_rcp85/grasslandarea_NorthAmerica_rcp85;
nicheprc_Oceania_rcp85=nichelandarea_Oceania_rcp85/grasslandarea_Oceania_rcp85;
nicheprc_Europe_rcp85=nichelandarea_Europe_rcp85/grasslandarea_Europe_rcp85;
nicheprc_world_rcp85=nichelandarea_world_rcp85./grasslandarea_world_rcp85;

nicheprc_Asia_rcp45 = nichelandarea_Asia_rcp45 / grasslandarea_Asia_rcp45;
nicheprc_Africa_rcp45 = nichelandarea_Africa_rcp45 / grasslandarea_Africa_rcp45;
nicheprc_SouthAmerica_rcp45 = nichelandarea_SouthAmerica_rcp45 / grasslandarea_SouthAmerica_rcp45;
nicheprc_NorthAmerica_rcp45 = nichelandarea_NorthAmerica_rcp45 / grasslandarea_NorthAmerica_rcp45;
nicheprc_Oceania_rcp45 = nichelandarea_Oceania_rcp45 / grasslandarea_Oceania_rcp45;
nicheprc_Europe_rcp45 = nichelandarea_Europe_rcp45 / grasslandarea_Europe_rcp45;
nicheprc_world_rcp45 = nichelandarea_world_rcp45 / grasslandarea_world_rcp45;

nicheprc_Asia_rcp60 = nichelandarea_Asia_rcp60 / grasslandarea_Asia_rcp60;
nicheprc_Africa_rcp60 = nichelandarea_Africa_rcp60 / grasslandarea_Africa_rcp60;
nicheprc_SouthAmerica_rcp60 = nichelandarea_SouthAmerica_rcp60 / grasslandarea_SouthAmerica_rcp60;
nicheprc_NorthAmerica_rcp60 = nichelandarea_NorthAmerica_rcp60 / grasslandarea_NorthAmerica_rcp60;
nicheprc_Oceania_rcp60 = nichelandarea_Oceania_rcp60 / grasslandarea_Oceania_rcp60;
nicheprc_Europe_rcp60 = nichelandarea_Europe_rcp60 / grasslandarea_Europe_rcp60;
nicheprc_world_rcp60 = nichelandarea_world_rcp60 / grasslandarea_world_rcp60;

%% 加上present（2015年）情况的
load('nichelandarea.mat')
figure_rcp_value=[
    nichelandarea_Asia,nichelandarea_Asia_rcp26,nichelandarea_Asia_rcp45,nichelandarea_Asia_rcp60,nichelandarea_Asia_rcp85;
    nichelandarea_Africa,nichelandarea_Africa_rcp26,nichelandarea_Africa_rcp45,nichelandarea_Africa_rcp60,nichelandarea_Africa_rcp85;
    nichelandarea_SouthAmerica,nichelandarea_SouthAmerica_rcp26,nichelandarea_SouthAmerica_rcp45,nichelandarea_SouthAmerica_rcp60,nichelandarea_SouthAmerica_rcp85;
    nichelandarea_NorthAmerica,nichelandarea_NorthAmerica_rcp26,nichelandarea_NorthAmerica_rcp45,nichelandarea_NorthAmerica_rcp60,nichelandarea_NorthAmerica_rcp85;
    nichelandarea_Oceania,nichelandarea_Oceania_rcp26,nichelandarea_Oceania_rcp45,nichelandarea_Oceania_rcp60,nichelandarea_Oceania_rcp85;
    nichelandarea_Europe,nichelandarea_Europe_rcp26,nichelandarea_Europe_rcp45,nichelandarea_Europe_rcp60,nichelandarea_Europe_rcp85]

% 不加present的
figure_rcp_value_1=[
    nichelandarea_Asia_rcp26,nichelandarea_Asia_rcp45,nichelandarea_Asia_rcp60,nichelandarea_Asia_rcp85;
    nichelandarea_Africa_rcp26,nichelandarea_Africa_rcp45,nichelandarea_Africa_rcp60,nichelandarea_Africa_rcp85;
    nichelandarea_SouthAmerica_rcp26,nichelandarea_SouthAmerica_rcp45,nichelandarea_SouthAmerica_rcp60,nichelandarea_SouthAmerica_rcp85;
    nichelandarea_NorthAmerica_rcp26,nichelandarea_NorthAmerica_rcp45,nichelandarea_NorthAmerica_rcp60,nichelandarea_NorthAmerica_rcp85;
    nichelandarea_Oceania_rcp26,nichelandarea_Oceania_rcp45,nichelandarea_Oceania_rcp60,nichelandarea_Oceania_rcp85;
    nichelandarea_Europe_rcp26,nichelandarea_Europe_rcp45,nichelandarea_Europe_rcp60,nichelandarea_Europe_rcp85]
%% save the data
save('./grazingniche/matdata/figure_rcp_value.mat','figure_rcp_value');     

%% bar with present
figure_rcp_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}

bar(figure_rcp_value)
% Set x-axis labels
set(gca, 'xticklabel', figure_rcp_caption);
% Add a title and labels
title('Future GN area of continents');
xlabel('Region');
ylabel('Grassland area');
legend('present','rcp2.6', 'rcp4.5','rcp6.0','rcp8.5');
%% bar without present
figure_rcp_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}

bar(figure_rcp_value_1)
% Set x-axis labels
set(gca, 'xticklabel', figure_rcp_caption);
% Add a title and labels
title('Future GN area of continents');
xlabel('Region');
ylabel('Grassland area');
legend('rcp2.6', 'rcp4.5','rcp6.0','rcp8.5');
% %% turnover rate suggested by Max

% %% turnover rate suggested by Max
% %% Step 1.Find the relevant variables
% % 要找到nicheland_world和nicheland_world_rcp26的根变量，这两个已经是有面积的了，这样再去插值升尺度或者降尺度就麻烦了。
% %nicheland_world的根变量是landuse_coupled1
% %nicheland_world_rcp26的根变量是futureniche_struct.rcp26
% 
% imagesc(landuse_coupled1) %美国在左边，值在100以内，1800*3600
% colorbar
% imagesc(futureniche_struct.rcp26) %美国在右边，值在1以内，160*320
% colorbar
% 
% %% Step 2.Unifying the resolution, unit, also unifying the coordinates 
% % Unifying coordination system
% halfSize = size(landuse_coupled1, 2) / 2;
% landuse_coupled1_resize1 = circshift(landuse_coupled1, [0, halfSize]);
% % Unifying unit (into percentage of grassland per pixel, 0-1)
% landuse_coupled1_resize2=landuse_coupled1_resize1/100
% % Unifying resolution
% landuse_coupled1_resize3 = imresize(landuse_coupled1_resize2, [160,320], 'bilinear');
% %% Test if they are unified
% imagesc(landuse_coupled1_resize3) %美国在右边，值0-1，160*320
% colorbar
% 
% imagesc(futureniche_struct.rcp26) %美国在右边，值0-1，160*320
% colorbar
% 
% %% Step 3.Turn percentage of area per pixel into area per pixel
% 
% compare_futureniche_rcp26=futureniche_struct.rcp26.*worldarea;
% compare_futureniche_rcp45=futureniche_struct.rcp45.*worldarea;
% compare_futureniche_rcp60=futureniche_struct.rcp60.*worldarea;
% compare_futureniche_rcp85=futureniche_struct.rcp85.*worldarea;
% 
% compare_niche=landuse_coupled1_resize3.*worldarea
% 
% 
% 
% %% Step 4. Compare them
% compare_turnover_rcp26=compare_futureniche_rcp26-compare_niche
% compare_turnover_rcp45=compare_futureniche_rcp45-compare_niche
% compare_turnover_rcp60=compare_futureniche_rcp60-compare_niche
% compare_turnover_rcp85=compare_futureniche_rcp85-compare_niche
% 
% save('./grazingniche/matdata/turnover.mat',"compare_turnover_rcp26","compare_turnover_rcp45","compare_turnover_rcp60","compare_turnover_rcp85")
% load('turnover.mat')% this .mat contains compare_turnover_rcp85 and other scenarios variable
% imagesc(compare_niche) 
% colorbar
% 
% imagesc(compare_futureniche_rcp26)
% colorbar
% 
% sum(sum(compare_niche))%结果是7.5e6
% sum(sum(compare_futureniche_rcp26))%结果是9.1e6
% sum(sum(compare_futureniche_rcp85))%结果是8.6e6
% 
% % 所以最后的结果是未来的niche会增加？？？？还增加了挺多？？
% % 是不是算错了？？？？
% 
% %% Step 5: Calculate turnover rate for each continent
% 
% % This is the net value. Net turnover rate
% turnover_total_rcp26=sum(sum(compare_turnover_rcp26));
% 
% % How much will be gained in new areas?
% turnover_increase_rcp26=compare_turnover_rcp26;
% turnover_increase_rcp26(turnover_increase_rcp26<0)=0
% turnover_increase_value_rcp26=sum(sum(turnover_increase_rcp26))
% 
% % How much will be lost?
% turnover_decrease_rcp26=compare_turnover_rcp26;
% turnover_decrease_rcp26(turnover_decrease_rcp26>0)=0
% turnover_decrease_value_rcp26=sum(sum(turnover_decrease_rcp26))
% 
% % How much in each continent?
% turnover_africa_rcp26=sum(sum(compare_turnover_rcp26(Africamask==1)));
% turnover_increase_africa_rcp26=sum(sum(turnover_increase_rcp26((Africamask==1))))
% turnover_decrease_africa_rcp26=sum(sum(turnover_decrease_rcp26((Africamask==1))))
% 
% %% Batch making for all continents and all rcps
% % Define RCPs and Continents
% rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
% continents = {'Asia', 'Europe', 'Oceania', 'Africa', 'NorthAmerica', 'SouthAmerica'};
% 
% % Initialize structures to store results
% turnover_total = struct();
% turnover_increase_value = struct();
% turnover_decrease_value = struct();
% turnover_continents = struct();
% 
% % Loop over each RCP scenario
% for i = 1:length(rcps)
%     rcp = rcps{i};
%     compare_turnover = eval(['compare_turnover_' rcp]); % Dynamically load the variable
% 
%     % Calculate total turnover for the current RCP
%     turnover_total.(rcp) = sum(sum(compare_turnover));
% 
%     % Calculate turnover increase and decrease
%     turnover_increase = compare_turnover;
%     turnover_increase(turnover_increase < 0) = 0;
%     turnover_increase_value.(rcp) = sum(sum(turnover_increase));
% 
%     turnover_decrease = compare_turnover;
%     turnover_decrease(turnover_decrease > 0) = 0;
%     turnover_decrease_value.(rcp) = sum(sum(turnover_decrease));
% 
%     % Calculate turnover for each continent
%     for j = 1:length(continents)
%         continent = continents{j};
%         mask = eval([continent 'mask']); % Dynamically use the mask variable
% 
%         turnover_continents.(continent).(rcp).total = sum(sum(compare_turnover(mask == 1)));
%         turnover_continents.(continent).(rcp).increase = sum(sum(turnover_increase(mask == 1)));
%         turnover_continents.(continent).(rcp).decrease = sum(sum(turnover_decrease(mask == 1)));
%     end
% end
% 
% % Display results (optional)
% disp(turnover_total);
% disp(turnover_increase_value);
% disp(turnover_decrease_value);
% disp(turnover_continents);
% 
% %% Creat bar chart
% %% all rcps
% 
% bar([turnover_continents.Africa.rcp26.increase,turnover_continents.Africa.rcp26.decrease,turnover_continents.Africa.rcp45.increase,turnover_continents.Africa.rcp45.decrease,turnover_continents.Africa.rcp60.increase,turnover_continents.Africa.rcp60.decrease,turnover_continents.Africa.rcp85.increase,turnover_continents.Africa.rcp85.decrease])
% bar([turnover_continents.Asia.rcp26.increase,turnover_continents.Asia.rcp26.decrease,turnover_continents.Asia.rcp45.increase,turnover_continents.Asia.rcp45.decrease,turnover_continents.Asia.rcp60.increase,turnover_continents.Asia.rcp60.decrease,turnover_continents.Asia.rcp85.increase,turnover_continents.Asia.rcp85.decrease])
% bar([turnover_continents.SouthAmerica.rcp26.increase,turnover_continents.SouthAmerica.rcp26.decrease,turnover_continents.SouthAmerica.rcp45.increase,turnover_continents.SouthAmerica.rcp45.decrease,turnover_continents.SouthAmerica.rcp60.increase,turnover_continents.SouthAmerica.rcp60.decrease,turnover_continents.SouthAmerica.rcp85.increase,turnover_continents.SouthAmerica.rcp85.decrease])
% bar([turnover_continents.NorthAmerica.rcp26.increase,turnover_continents.NorthAmerica.rcp26.decrease,turnover_continents.NorthAmerica.rcp45.increase,turnover_continents.NorthAmerica.rcp45.decrease,turnover_continents.NorthAmerica.rcp60.increase,turnover_continents.NorthAmerica.rcp60.decrease,turnover_continents.NorthAmerica.rcp85.increase,turnover_continents.NorthAmerica.rcp85.decrease])
% bar([turnover_continents.Europe.rcp26.increase,turnover_continents.Europe.rcp26.decrease,turnover_continents.Europe.rcp45.increase,turnover_continents.Europe.rcp45.decrease,turnover_continents.Europe.rcp60.increase,turnover_continents.Europe.rcp60.decrease,turnover_continents.Europe.rcp85.increase,turnover_continents.Europe.rcp85.decrease])
% bar([turnover_continents.Oceania.rcp26.increase,turnover_continents.Oceania.rcp26.decrease,turnover_continents.Oceania.rcp45.increase,turnover_continents.Oceania.rcp45.decrease,turnover_continents.Oceania.rcp60.increase,turnover_continents.Oceania.rcp60.decrease,turnover_continents.Oceania.rcp85.increase,turnover_continents.Oceania.rcp85.decrease])
% 
% %% Only rcp85
% 
% % Set up the tiled layout for the subplots
% figure;
% tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % 2x3 layout
% 
% % Define the colors for gain and loss
% colors = [0 1 0; 1 0 0]; % Green for gain, Red for loss
% 
% % Plot each continent in a separate subplot
% nexttile
% bar([turnover_continents.Africa.rcp85.increase, turnover_continents.Africa.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
% title("Africa")
% 
% nexttile
% bar([turnover_continents.Asia.rcp85.increase, turnover_continents.Asia.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
% title("Asia")
% 
% nexttile
% bar([turnover_continents.SouthAmerica.rcp85.increase, turnover_continents.SouthAmerica.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
% title("South America")
% 
% nexttile
% bar([turnover_continents.NorthAmerica.rcp85.increase, turnover_continents.NorthAmerica.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
% title("North America")
% 
% nexttile
% bar([turnover_continents.Oceania.rcp85.increase, turnover_continents.Oceania.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
% title("Oceania")
% 
% nexttile
% bar([turnover_continents.Europe.rcp85.increase, turnover_continents.Europe.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
% title("Europe")
% 
% 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_rcp85_allcontinents.svg');
% 
% %% Create maps
% %% Custom colormap just for grazing turnover rate (central white, two spectrum red and green.)
% % This is actually a very useful color pallet. Can adjst for different colors
% num_colors = 256;  % Total number of colors for smooth transition
% 
% % Define specific color points: dark red, white, dark green
% color_points = [ 
%     0.6, 0.0, 0.0;  % Dark Red
%     1.0, 1.0, 1.0;  % White
%     0.1, 0.8, 0.3;  % Dark Green
% ];
% 
% % Interpolate to create a smooth transition colormap
% reds_to_white = [linspace(color_points(1,1), color_points(2,1), num_colors/2)', ...
%                  linspace(color_points(1,2), color_points(2,2), num_colors/2)', ...
%                  linspace(color_points(1,3), color_points(2,3), num_colors/2)'];
% 
% white_to_greens = [linspace(color_points(2,1), color_points(3,1), num_colors/2)', ...
%                    linspace(color_points(2,2), color_points(3,2), num_colors/2)', ...
%                    linspace(color_points(2,3), color_points(3,3), num_colors/2)'];
% 
% % Combine the two gradients into a single colormap
% custom_colormap = [reds_to_white; white_to_greens];
% %% Try one plot first
% figure2 = figure;
% sgtitle('Turnover', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% %custom_colormap = [1 1 1; addcolorplus(289)];
% custom_colormap = [reds_to_white; white_to_greens];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% caxis([-14000 14000]);
% c = colorbar;  
% R = georefcells([-90,90],[0,360],size(compare_turnover_rcp26));
% geoshow(flipud(compare_turnover_rcp85), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% set(gcf, 'Position',  [584,449,684,537])
% 
% % set(gcf,'renderer','painters');
% % it seems matlab is unable to produce real svg with this map
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover85.svg');
% 
% %% Four rcps one plot
% % Define the RCP scenarios
% rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
% titles = {'Turnover for RCP 2.6', 'Turnover for RCP 4.5', 'Turnover for RCP 6.0', 'Turnover for RCP 8.5'};
% 
% % Create a figure for the subplots
% figure;
% sgtitle('Turnover Across Different RCP Scenarios', 'FontSize', 16);
% 
% % Loop over each RCP scenario to create subplots
% for i = 1:length(rcps)
%     rcp = rcps{i};
%     title_text = titles{i};
%     
%     % Create a subplot for each RCP scenario
%     subplot(2, 2, i);
%     
%     % Set up the world map for the current subplot
%     ax = worldmap('world');
%     setm(ax, 'FFaceColor', [1 1 1]);  % Set the map's background color
%     
%     % Apply the custom colormap
%     set(gcf, 'Colormap', custom_colormap);
%     
%     % Set color axis limits to ensure zero aligns with white
%     caxis([-14000 14000]);  % Adjust according to your data's range
%     
%     % Load and display the data for the current RCP
%     compare_turnover = eval(['compare_turnover_' rcp]); % Dynamically load the variable
%     R = georefcells([-90, 90], [0, 360], size(compare_turnover));
%     geoshow(flipud(compare_turnover), R, 'DisplayType', 'texturemap');
%     
%     % Load and plot coastlines
%     load coastlines;
%     plotm(coastlat, coastlon, 'Color', 'black');
%     
%     % Add a title to each subplot
%     title(title_text, 'FontSize', 12);
%     
%     % Adjust subplot for tight and compact style
% end
% 
% % Adjust the figure's size for a tight, compact display
% set(gcf, 'Position', [100, 100, 1000, 600]);
% 
% % Display a colorbar that applies to all subplots
% colorbar('Position', [0.93, 0.1, 0.02, 0.8]);  % Custom position for a single colorbar
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_allrcp.svg');
% 
% save('./grazingniche/matdata/turnover.mat','compare_turnover_rcp26','compare_turnover_rcp45','compare_turnover_rcp60','compare_turnover_rcp85');     
