%% Original code
%% carbon loss due to grazing
load('AGBCBGBC.mat')
load('SOC.mat')
load("worldarea.mat")
load('overgraze_niche.mat')
SOC0_200=double(SOC0_200);
AGBC=double(AGBC);
BGBC=double(BGBC);

carbonloss_SOC_map=(overgraze_insideniche_logic.*SOC0_200.*0.1.*grassland_env./100.*worldarea.*100)+((overgraze_outsideniche_logic.*SOC0_200.*0.2.*grassland_env./100.*worldarea.*100));
carbonloss_AGBC_map=overgraze_logic.*AGBC.*grassland_env./100.*worldarea.*100;
carbonloss_BGBC_map=overgraze_logic.*BGBC.*grassland_env./100.*worldarea.*100;

carbonloss_map=carbonloss_SOC_map+carbonloss_AGBC_map+carbonloss_BGBC_map

carbonloss_outniche=sum(carbonloss_map(landuse_coupled1>1))

carbonloss_SOC=sum(sum(carbonloss_SOC_map));
carbonloss_AGBC=sum(sum(carbonloss_AGBC_map));
carbonloss_BGBC=sum(sum(carbonloss_BGBC_map));
carbonloss_all=carbonloss_SOC+carbonloss_AGBC+carbonloss_BGBC
%% Montecarlo analysis
% Define the number of simulations
num_simulations = 10000;

% Initialize vectors to store the total carbon loss for each simulation
total_carbonloss_SOC = zeros(1, num_simulations);
total_carbonloss_AGBC = zeros(1, num_simulations);
total_carbonloss_BGBC = zeros(1, num_simulations);
total_carbonloss_all = zeros(1, num_simulations);

% Parameters for the normal distribution
mean_inside_niche = 0.1;
mean_outside_niche = 0.2;
std_dev = 0.05; % Adjust this based on your assessment

for i = 1:num_simulations
    % Sample the percentages from the normal distributions
    perc_inside_niche = normrnd(mean_inside_niche, std_dev);
    perc_outside_niche = normrnd(mean_outside_niche, std_dev);
    
    % Calculate carbon loss maps for this simulation
    carbonloss_SOC_map = (overgraze_insideniche_logic .* SOC0_200 .* perc_inside_niche .* grassland_env ./ 100 .* worldarea .* 100) + ...
                         (overgraze_outsideniche_logic .* SOC0_200 .* perc_outside_niche .* grassland_env ./ 100 .* worldarea .* 100);
    carbonloss_AGBC_map = overgraze_logic .* AGBC .* grassland_env ./ 100 .* worldarea .* 100;
    carbonloss_BGBC_map = overgraze_logic .* BGBC .* grassland_env ./ 100 .* worldarea .* 100;
    
    % Calculate total carbon loss for this simulation
    total_carbonloss_SOC(i) = sum(sum(carbonloss_SOC_map));
    total_carbonloss_AGBC(i) = sum(sum(carbonloss_AGBC_map));
    total_carbonloss_BGBC(i) = sum(sum(carbonloss_BGBC_map));
    total_carbonloss_all(i) = total_carbonloss_SOC(i) + total_carbonloss_AGBC(i) + total_carbonloss_BGBC(i);
end

% After all simulations are done, analyze the results
% Calculate mean, standard deviation, and other statistics if needed
mean_total_carbonloss_all = mean(total_carbonloss_all);
std_total_carbonloss_all = std(total_carbonloss_all);

% For visualization, you can plot histograms of the results
% histogram(total_carbonloss_all);
% title('Distribution of Total Carbon Loss Across Simulations');
% xlabel('Total Carbon Loss');
% ylabel('Frequency');
%
% Calculate the means
mean_SOC = mean(total_carbonloss_SOC);
mean_AGBC = mean(total_carbonloss_AGBC);
mean_BGBC = mean(total_carbonloss_BGBC);
mean_all = mean(total_carbonloss_all);

% Calculate the standard deviations
std_SOC = std(total_carbonloss_SOC);
std_AGBC = std(total_carbonloss_AGBC);
std_BGBC = std(total_carbonloss_BGBC);
std_all = std(total_carbonloss_all);

% Number of simulations
n = length(total_carbonloss_SOC); % Assuming all arrays are the same length

% Calculate the margin of error (MOE) for 95% confidence interval
% Z-score for 95% confidence is approximately 1.96
Z = 1.96;
MOE_SOC = Z * (std_SOC / sqrt(n));
MOE_AGBC = Z * (std_AGBC / sqrt(n));
MOE_BGBC = Z * (std_BGBC / sqrt(n));
MOE_all = Z * (std_all / sqrt(n));

% Calculate the 95% confidence intervals
CI_SOC = [mean_SOC - MOE_SOC, mean_SOC + MOE_SOC];
CI_AGBC = [mean_AGBC - MOE_AGBC, mean_AGBC + MOE_AGBC];
CI_BGBC = [mean_BGBC - MOE_BGBC, mean_BGBC + MOE_BGBC];
CI_all = [mean_all - MOE_all, mean_all + MOE_all];

% Display the results
fprintf('95%% CI for Total Carbon Loss from SOC: [%f, %f]\n', CI_SOC(1), CI_SOC(2));
fprintf('95%% CI for Total Carbon Loss from AGBC: [%f, %f]\n', CI_AGBC(1), CI_AGBC(2));
fprintf('95%% CI for Total Carbon Loss from BGBC: [%f, %f]\n', CI_BGBC(1), CI_BGBC(2));
fprintf('95%% CI for Total Combined Carbon Loss: [%f, %f]\n', CI_all(1), CI_all(2));

%% Original code

load('livestockdensity.mat') %variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep')
load('grassland_env.mat')
load('AGB.mat')
load("worldarea.mat")

nicheland_world=landuse_coupled1.*worldarea./100;
nichelandarea_world=sum(sum(nicheland_world));

%
cattlegrass=resizecattle*455*0.03*365*0.001;
goatsgrass=resizegoats*91*0.03*365*0.001;
sheepgrass=resizesheep*91*0.03*365*0.001;

livestockgrass=cattlegrass+goatsgrass+sheepgrass;

livestockgrass(grassland_env==0)=0

livestockpercentage=(livestockgrass)./AGB

livestockpercentage(grassland_env<3)=0

livestockpercentage(livestockpercentage>1) = 1;
livestockpercentage(livestockpercentage<0) = 0;

%
overgrazed_logic=double(livestockpercentage >= 0.65 & livestockpercentage <= 1)
mediumgrazed_logic=double(livestockpercentage >= 0.20 & livestockpercentage <= 0.65)
lightgrazed_logic=double(livestockpercentage >= 0 & livestockpercentage <= 0.20)

% how much of world's grassland is overgrazed
overgraze_world = double(livestockpercentage >= 0.65 & livestockpercentage <= 1);
overgraze_world_area=sum(sum(overgraze_world.*grassland_env.*worldarea/100));
percentage_overgraze=overgraze_world_area./grasslandarea_world;
percentage_overgraze_niche=overgraze_world_area./nichelandarea_world;
%% Montecarlo analysis
% Define parameters for Monte Carlo simulation
num_simulations = 10000; % Number of Monte Carlo simulations
mean_consumption_rate = 0.03; % Mean consumption rate
std_dev_consumption_rate = 0.005; % Standard deviation of consumption rate, adjust based on your assessment

% Initialize an array to store the results of each simulation
overgraze_world_areas = zeros(num_simulations, 1);

for i = 1:num_simulations
    % Sample the consumption rate from a normal distribution
    sampled_consumption_rate = normrnd(mean_consumption_rate, std_dev_consumption_rate);
    
    % Calculate livestockgrass with the new sampled rate
    cattlegrass = resizecattle * 455 * sampled_consumption_rate * 365 * 0.001;
    goatsgrass = resizegoats * 91 * sampled_consumption_rate * 365 * 0.001;
    sheepgrass = resizesheep * 91 * sampled_consumption_rate * 365 * 0.001;
    livestockgrass = cattlegrass + goatsgrass + sheepgrass;
    
    % Apply your existing logic with the new livestockgrass
    livestockgrass(grassland_env == 0) = 0;
    livestockpercentage = livestockgrass ./ AGB;
    livestockpercentage(grassland_env < 3) = 0;
    livestockpercentage(livestockpercentage > 1) = 1;
    livestockpercentage(livestockpercentage < 0) = 0;
    
    % Calculate overgrazed area for this iteration
    overgraze_world = double(livestockpercentage >= 0.65 & livestockpercentage <= 1);
    overgraze_world_areas(i) = sum(sum(overgraze_world .* grassland_env .* worldarea / 100));
end

% Calculate the 95% confidence interval for overgraze_world_area
mean_overgraze_area = mean(overgraze_world_areas);
std_overgraze_area = std(overgraze_world_areas);
ci_lower_bound = mean_overgraze_area - 1.96 * std_overgraze_area / sqrt(num_simulations);
ci_upper_bound = mean_overgraze_area + 1.96 * std_overgraze_area / sqrt(num_simulations);

fprintf('95%% Confidence Interval for Overgrazed World Area: [%f, %f]\n', ci_lower_bound, ci_upper_bound);


