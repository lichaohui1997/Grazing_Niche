load('livestockdensity.mat')

load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
aggregate_sfcWind(aggregate_sfcWind>10)=NaN
livestock_all=resizecattle+resizesheep+resizegoats
%% determining modern livestock distribution niche
nexttile
imagesc(aggregate_tas)
nexttile
imagesc(resizecattle)
nexttile
imagesc(livestock)
nexttile
imagesc(aggregate_hurs)
colorbar
%%
sgtitle('present temperature', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(316)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
%c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[-180,180],size(aggregate_tas));
geoshow(flipud(aggregate_hurs), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])
%caxis([0,6000])
%%
vector_tas=reshape(aggregate_tas,[],1)
vector_pr=reshape(aggregate_pr,[],1)
vector_hurs=reshape(aggregate_hurs,[],1)
vector_sfcWind=reshape(aggregate_sfcWind,[],1)

vector_cattle=reshape(resizecattle,[],1)
vector_goats=reshape(resizegoats,[],1)
vector_sheep=reshape(resizesheep,[],1)
vector_livestock=reshape(livestock,[],1)

%% livestock
middletas = vector_tas(vector_livestock >1);
middlepr = vector_pr(vector_livestock >1);
middlesfcWind = vector_sfcWind(vector_livestock >1);
middlehurs = vector_hurs(vector_livestock >1);

%
lower_bound = prctile(middletas, 2.75);
upper_bound = prctile(middletas, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middletas, 0.5);
upper_bound = prctile(middletas, 99.5);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middletas, 5);
upper_bound = prctile(middletas, 95);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);

%
lower_bound = prctile(middlepr, 2.75);
upper_bound = prctile(middlepr, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middlepr, 0.5);
upper_bound = prctile(middlepr, 99.5);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middlepr, 5);
upper_bound = prctile(middlepr, 95);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);

%
lower_bound = prctile(middlesfcWind, 2.75);
upper_bound = prctile(middlesfcWind, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middlesfcWind, 0.5);
upper_bound = prctile(middlesfcWind, 99.5);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middlesfcWind, 5);
upper_bound = prctile(middlesfcWind, 95);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);

%
lower_bound = prctile(middlehurs, 2.75);
upper_bound = prctile(middlehurs, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middlehurs, 0.5);
upper_bound = prctile(middlehurs, 99.5);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
lower_bound = prctile(middlehurs, 5);
upper_bound = prctile(middlehurs, 95);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);

%% sheep
middletas = vector_tas(vector_sheep >1);
middlepr = vector_pr(vector_sheep >1);
middlesfcWind = vector_sfcWind(vector_sheep >1);
middlehurs = vector_hurs(vector_sheep >1);

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
%% goats
middletas = vector_tas(vector_goats > 1);
middlepr = vector_pr(vector_goats > 1);
middlesfcWind = vector_sfcWind(vector_goats > 1);
middlehurs = vector_hurs(vector_goats > 1);

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

%% cattle
middletas = vector_tas(vector_cattle >1);
middlepr = vector_pr(vector_cattle >1);
middlesfcWind = vector_sfcWind(vector_cattle >1);
middlehurs = vector_hurs(vector_cattle >1);

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
%% sheep
middletas = vector_tas(vector_sheep >2);
middlepr = vector_pr(vector_sheep >2);
middlesfcWind = vector_sfcWind(vector_sheep >2);
middlehurs = vector_hurs(vector_sheep >2);

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
%% goats
middletas = vector_tas(vector_goats > 2);
middlepr = vector_pr(vector_goats > 2);
middlesfcWind = vector_sfcWind(vector_goats > 2);
middlehurs = vector_hurs(vector_goats > 2);

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
% Flatten the matrices into vectors
pr_values = aggregate_pr(:); % Temperature per grid cell
cattle_counts = resizecattle(:); % Cattle counts per grid cell

% Filter out cells with zero cattle to focus on cells with cattle presence
cattle_present_indices = cattle_counts > 5;
pr_with_cattle = pr_values(cattle_present_indices);
cattle_with_cattle = cattle_counts(cattle_present_indices);

% Sort temperature and corresponding cattle counts
[sorted_temps, sort_index] = sort(pr_with_cattle);
sorted_cattle_counts = cattle_with_cattle(sort_index);

% Calculate cumulative distribution based on cattle weights
cumulative_cattle = cumsum(sorted_cattle_counts) / sum(sorted_cattle_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_cattle >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_cattle >= 0.975, 1));

% Display the results
fprintf('The lower pr threshold for the central 95%% of cattle distribution is %.2f\n', lower_threshold);
fprintf('The upper pr threshold for the central 95%% of cattle distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
tas_values = aggregate_tas(:); % Temperature per grid cell
cattle_counts = resizecattle(:); % Cattle counts per grid cell

% Filter out cells with zero cattle to focus on cells with cattle presence
cattle_present_indices = cattle_counts > 5;
tas_with_cattle = tas_values(cattle_present_indices);
cattle_with_cattle = cattle_counts(cattle_present_indices);

% Sort temperature and corresponding cattle counts
[sorted_temps, sort_index] = sort(tas_with_cattle);
sorted_cattle_counts = cattle_with_cattle(sort_index);

% Calculate cumulative distribution based on cattle weights
cumulative_cattle = cumsum(sorted_cattle_counts) / sum(sorted_cattle_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_cattle >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_cattle >= 0.975, 1));

% Display the results
fprintf('The lower tas threshold for the central 95%% of cattle distribution is %.2f\n', lower_threshold);
fprintf('The upper tas threshold for the central 95%% of cattle distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
hurs_values = aggregate_hurs(:); % Temperature per grid cell
cattle_counts = resizecattle(:); % Cattle counts per grid cell

% Filter out cells with zero cattle to focus on cells with cattle presence
cattle_present_indices = cattle_counts > 5;
hurs_with_cattle = hurs_values(cattle_present_indices);
cattle_with_cattle = cattle_counts(cattle_present_indices);

% Sort temperature and corresponding cattle counts
[sorted_temps, sort_index] = sort(hurs_with_cattle);
sorted_cattle_counts = cattle_with_cattle(sort_index);

% Calculate cumulative distribution based on cattle weights
cumulative_cattle = cumsum(sorted_cattle_counts) / sum(sorted_cattle_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_cattle >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_cattle >= 0.975, 1));

% Display the results
fprintf('The lower hurs threshold for the central 95%% of cattle distribution is %.2f\n', lower_threshold);
fprintf('The upper hurs threshold for the central 95%% of cattle distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
sfcWind_values = aggregate_sfcWind(:); % Temperature per grid cell
cattle_counts = resizecattle(:); % Cattle counts per grid cell

% Filter out cells with zero cattle to focus on cells with cattle presence
cattle_present_indices = cattle_counts > 5;
sfcWind_with_cattle = sfcWind_values(cattle_present_indices);
cattle_with_cattle = cattle_counts(cattle_present_indices);

% Sort temperature and corresponding cattle counts
[sorted_temps, sort_index] = sort(sfcWind_with_cattle);
sorted_cattle_counts = cattle_with_cattle(sort_index);

% Calculate cumulative distribution based on cattle weights
cumulative_cattle = cumsum(sorted_cattle_counts) / sum(sorted_cattle_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_cattle >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_cattle >= 0.975, 1));

% Display the results
fprintf('The lower wfcWind threshold for the central 95%% of cattle distribution is %.2f\n', lower_threshold);
fprintf('The upper sfcWind threshold for the central 95%% of cattle distribution is %.2f\n', upper_threshold);
%%
%%
% Flatten the matrices into vectors
tas_values = aggregate_tas(:); % Temperature per grid cell
goats_counts = resizegoats(:); % Cattle counts per grid cell

% Filter out cells with zero goats to focus on cells with goats presence
goats_present_indices = goats_counts > 5;
tas_with_goats = tas_values(goats_present_indices);
goats_with_goats = goats_counts(goats_present_indices);

% Sort temperature and corresponding goats counts
[sorted_temps, sort_index] = sort(tas_with_goats);
sorted_goats_counts = goats_with_goats(sort_index);

% Calculate cumulative distribution based on goats weights
cumulative_goats = cumsum(sorted_goats_counts) / sum(sorted_goats_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_goats >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_goats >= 0.975, 1));

% Display the results
fprintf('The lower tas threshold for the central 95%% of goats distribution is %.2f\n', lower_threshold);
fprintf('The upper tas threshold for the central 95%% of goats distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
hurs_values = aggregate_hurs(:); % Temperature per grid cell
goats_counts = resizegoats(:); % Cattle counts per grid cell

% Filter out cells with zero goats to focus on cells with goats presence
goats_present_indices = goats_counts > 5;
hurs_with_goats = hurs_values(goats_present_indices);
goats_with_goats = goats_counts(goats_present_indices);

% Sort temperature and corresponding goats counts
[sorted_temps, sort_index] = sort(hurs_with_goats);
sorted_goats_counts = goats_with_goats(sort_index);

% Calculate cumulative distribution based on goats weights
cumulative_goats = cumsum(sorted_goats_counts) / sum(sorted_goats_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_goats >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_goats >= 0.975, 1));

% Display the results
fprintf('The lower hurs threshold for the central 95%% of goats distribution is %.2f\n', lower_threshold);
fprintf('The upper hurs threshold for the central 95%% of goats distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
sfcWind_values = aggregate_sfcWind(:); % Temperature per grid cell
goats_counts = resizegoats(:); % Cattle counts per grid cell

% Filter out cells with zero goats to focus on cells with goats presence
goats_present_indices = goats_counts > 5;
sfcWind_with_goats = sfcWind_values(goats_present_indices);
goats_with_goats = goats_counts(goats_present_indices);

% Sort temperature and corresponding goats counts
[sorted_temps, sort_index] = sort(sfcWind_with_goats);
sorted_goats_counts = goats_with_goats(sort_index);

% Calculate cumulative distribution based on goats weights
cumulative_goats = cumsum(sorted_goats_counts) / sum(sorted_goats_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_goats >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_goats >= 0.975, 1));

% Display the results
fprintf('The lower sfcWind threshold for the central 95%% of goats distribution is %.2f\n', lower_threshold);
fprintf('The upper sfcWind threshold for the central 95%% of goats distribution is %.2f\n', upper_threshold);
%%
%%
% Flatten the matrices into vectors
pr_values = aggregate_pr(:); % Temperature per grid cell
goats_counts = resizegoats(:); % Cattle counts per grid cell

% Filter out cells with zero goats to focus on cells with goats presence
goats_present_indices = goats_counts > 5;
pr_with_goats = pr_values(goats_present_indices);
goats_with_goats = goats_counts(goats_present_indices);

% Sort temperature and corresponding goats counts
[sorted_temps, sort_index] = sort(pr_with_goats);
sorted_goats_counts = goats_with_goats(sort_index);

% Calculate cumulative distribution based on goats weights
cumulative_goats = cumsum(sorted_goats_counts) / sum(sorted_goats_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_goats >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_goats >= 0.975, 1));

% Display the results
fprintf('The lower pr threshold for the central 95%% of goats distribution is %.2f\n', lower_threshold);
fprintf('The upper pr threshold for the central 95%% of goats distribution is %.2f\n', upper_threshold);
%% sheep
%%
% Flatten the matrices into vectors
tas_values = aggregate_tas(:); % Temperature per grid cell
sheep_counts = resizesheep(:); % Cattle counts per grid cell

% Filter out cells with zero sheep to focus on cells with sheep presence
sheep_present_indices = sheep_counts > 5;
tas_with_sheep = tas_values(sheep_present_indices);
sheep_with_sheep = sheep_counts(sheep_present_indices);

% Sort temperature and corresponding sheep counts
[sorted_temps, sort_index] = sort(tas_with_sheep);
sorted_sheep_counts = sheep_with_sheep(sort_index);

% Calculate cumulative distribution based on sheep weights
cumulative_sheep = cumsum(sorted_sheep_counts) / sum(sorted_sheep_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_sheep >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_sheep >= 0.975, 1));

% Display the results
fprintf('The lower tas threshold for the central 95%% of sheep distribution is %.2f\n', lower_threshold);
fprintf('The upper tas threshold for the central 95%% of sheep distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
hurs_values = aggregate_hurs(:); % Temperature per grid cell
sheep_counts = resizesheep(:); % Cattle counts per grid cell

% Filter out cells with zero sheep to focus on cells with sheep presence
sheep_present_indices = sheep_counts > 5;
hurs_with_sheep = hurs_values(sheep_present_indices);
sheep_with_sheep = sheep_counts(sheep_present_indices);

% Sort temperature and corresponding sheep counts
[sorted_temps, sort_index] = sort(hurs_with_sheep);
sorted_sheep_counts = sheep_with_sheep(sort_index);

% Calculate cumulative distribution based on sheep weights
cumulative_sheep = cumsum(sorted_sheep_counts) / sum(sorted_sheep_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_sheep >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_sheep >= 0.975, 1));

% Display the results
fprintf('The lower hurs threshold for the central 95%% of sheep distribution is %.2f\n', lower_threshold);
fprintf('The upper hurs threshold for the central 95%% of sheep distribution is %.2f\n', upper_threshold);
%%
% Flatten the matrices into vectors
sfcWind_values = aggregate_sfcWind(:); % Temperature per grid cell
sheep_counts = resizesheep(:); % Cattle counts per grid cell

% Filter out cells with zero sheep to focus on cells with sheep presence
sheep_present_indices = sheep_counts > 5;
sfcWind_with_sheep = sfcWind_values(sheep_present_indices);
sheep_with_sheep = sheep_counts(sheep_present_indices);

% Sort temperature and corresponding sheep counts
[sorted_temps, sort_index] = sort(sfcWind_with_sheep);
sorted_sheep_counts = sheep_with_sheep(sort_index);

% Calculate cumulative distribution based on sheep weights
cumulative_sheep = cumsum(sorted_sheep_counts) / sum(sorted_sheep_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_sheep >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_sheep >= 0.975, 1));

% Display the results
fprintf('The lower sfcWind threshold for the central 95%% of sheep distribution is %.2f\n', lower_threshold);
fprintf('The upper sfcWind threshold for the central 95%% of sheep distribution is %.2f\n', upper_threshold);
%%
%%
% Flatten the matrices into vectors
pr_values = aggregate_pr(:); % Temperature per grid cell
sheep_counts = resizesheep(:); % Cattle counts per grid cell

% Filter out cells with zero sheep to focus on cells with sheep presence
sheep_present_indices = sheep_counts > 5;
pr_with_sheep = pr_values(sheep_present_indices);
sheep_with_sheep = sheep_counts(sheep_present_indices);

% Sort temperature and corresponding sheep counts
[sorted_temps, sort_index] = sort(pr_with_sheep);
sorted_sheep_counts = sheep_with_sheep(sort_index);

% Calculate cumulative distribution based on sheep weights
cumulative_sheep = cumsum(sorted_sheep_counts) / sum(sorted_sheep_counts);

% Find the 2.5th and 97.5th percentiles in the weighted cumulative distribution
lower_threshold = sorted_temps(find(cumulative_sheep >= 0.025, 1));
upper_threshold = sorted_temps(find(cumulative_sheep >= 0.975, 1));

% Display the results
fprintf('The lower pr threshold for the central 95%% of sheep distribution is %.2f\n', lower_threshold);
fprintf('The upper pr threshold for the central 95%% of sheep distribution is %.2f\n', upper_threshold);
