%% This script calculates modern niche decided by livestock data and MRI data
%% determining modern livestock distribution niche
% make sure that the data are aligned
load('livestock.mat')
imagesc(data_tas_rcp26_present)
imagesc(resizecattle)
resizecattle = resizecattle(:, [ceil(end/2+1):end, 1:floor(end/2)]);
resizegoats = resizegoats(:, [ceil(end/2+1):end, 1:floor(end/2)]);
resizesheep = resizesheep(:, [ceil(end/2+1):end, 1:floor(end/2)]);

resizecattle=imresize(resizecattle,[160,320])
resizegoats=imresize(resizegoats,[160,320])
resizesheep=imresize(resizesheep,[160,320])

imagesc(resizecattle)
imagesc(resizesheep)
imagesc(aggregate_hurs)
colorbar

%%
vector_tas=reshape(data_tas_rcp26_present,[],1)
vector_pr=reshape(data_pr_rcp26_present,[],1)
vector_hurs=reshape(data_hurs_rcp26_present,[],1)
vector_sfcWind=reshape(data_sfcWind_rcp26_present,[],1)

vector_cattle=reshape(resizecattle,[],1)
vector_goats=reshape(resizegoats,[],1)
vector_sheep=reshape(resizesheep,[],1)

%% cattle
middletas = vector_tas(vector_cattle >2);
middlepr = vector_pr(vector_cattle >2);
middlesfcWind = vector_sfcWind(vector_cattle >2);
middlehurs = vector_hurs(vector_cattle >2);

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
