close all;
%% In this code I am downloading and processing 4 new climate models, but from the same
% institution from the model in run1_historical_simulation, so there is
% similarity 
%% The downloading code is emitted.
%% I have downloaded r1i1p121 to 124 into my directory. I will now process them 
% this section of code takes around 30s to run. You have to start from this
% section of code in this .mat file
base_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/climatenew';
fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/climatenew';

% Define list of variables and models
variables = {'pr', 'tas', 'sfcWind', 'hurs'};
versions = {'r1i1p121', 'r1i1p122', 'r1i1p123', 'r1i1p124'};

% Loop over each variable
for var = variables
    % Loop over each model version
    for ver = versions
        % Get a list of all .nc files for the current variable and version
        files = dir(fullfile(base_path, [var{:} '_*_' ver{:} '_*.nc']));
        
        % Loop over each file
        for idx = 1:length(files)
            % Full path to the current file
            filename = fullfile(files(idx).folder, files(idx).name);
            
            % Load the data
            data = ncread(filename, var{:});
            
           % here I am making the above 100 values of relative humidity 100
        if strcmp(var, 'hurs')
            data(data > 100) = 100;
        end


        % Compute the mean for the first time point and convert to Celsius if temperature data
        if strcmp(var, 'tas') % If the variable is 'tas', convert from Kelvin to Celsius
            data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
            mean_data = rot90(data(:,:,1)) - 273.15;
            mean_data=imresize(mean_data,[90,180],'nearest')
            mean_data=reshape(mean_data,[],1)
        elseif strcmp(var, 'pr') % If the variable is 'pr', sum up monthly precipitations
            data = squeeze(sum(reshape(data, 144, 90, 12, []), 3)); % Turn monthly aggregate to yearly aggregate
            mean_data = rot90(data(:,:,1));
            mean_data=86400*30*mean_data
            mean_data=imresize(mean_data,[90,180],'nearest')
            mean_data=reshape(mean_data,[],1)
        else
            data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
            mean_data = rot90(data(:,:,1));
            mean_data=imresize(mean_data,[90,180],'nearest')
            mean_data=reshape(mean_data,[],1)

        end
           
            % Define the name for the new mean_data variable and for the saved files
            varName = ['mean_' var{:} '_p' ver{1}(6:end) '_' num2str(idx)];
            
            % Store the mean_data to a variable
            eval([varName ' = mean_data;']);
                      
            % Save the variable to a .mat file
            save_path = fullfile(base_path, [varName '.mat']);
         
            % Save the variable
            save(save_path, varName);      
        end
    end
end
%% Now I am putting my data into variables by only taking the odd years so I am only taking into years
% 1000, 1100, 1200, .. 1800 into account and not years 1050,etc.
% remember the data structure in run1_historical_simulation?
% For sfcWind
sfcWind_p_124 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['sfcWind_p_124 = [sfcWind_p_124; mean_sfcWind_p124_' num2str(i) '];']);
end

% For pr
pr_p_124 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['pr_p_124 = [pr_p_124; mean_pr_p124_' num2str(i) '];']);
end

% For tas
tas_p_124 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['tas_p_124 = [tas_p_124; mean_tas_p124_' num2str(i) '];']);
end

% For hurs
hurs_p_124 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['hurs_p_124 = [hurs_p_124; mean_hurs_p124_' num2str(i) '];']);
end
%% repeat the process for p123
% For sfcWind
sfcWind_p_123 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['sfcWind_p_123 = [sfcWind_p_123; mean_sfcWind_p123_' num2str(i) '];']);
end

% For pr
pr_p_123 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['pr_p_123 = [pr_p_123; mean_pr_p123_' num2str(i) '];']);
end

% For tas
tas_p_123 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['tas_p_123 = [tas_p_123; mean_tas_p123_' num2str(i) '];']);
end

% For hurs
hurs_p_123 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['hurs_p_123 = [hurs_p_123; mean_hurs_p123_' num2str(i) '];']);
end
%% repeat the process for 122
% For sfcWind
sfcWind_p_122 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['sfcWind_p_122 = [sfcWind_p_122; mean_sfcWind_p122_' num2str(i) '];']);
end

% For pr
pr_p_122 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['pr_p_122 = [pr_p_122; mean_pr_p122_' num2str(i) '];']);
end

% For tas
tas_p_122 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['tas_p_122 = [tas_p_122; mean_tas_p122_' num2str(i) '];']);
end

% For hurs
hurs_p_122 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['hurs_p_122 = [hurs_p_122; mean_hurs_p122_' num2str(i) '];']);
end
%% repeat the process for p121
% For sfcWind
sfcWind_p_121 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['sfcWind_p_121 = [sfcWind_p_121; mean_sfcWind_p121_' num2str(i) '];']);
end

% For pr
pr_p_121 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['pr_p_121 = [pr_p_121; mean_pr_p121_' num2str(i) '];']);
end

% For tas
tas_p_121 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['tas_p_121 = [tas_p_121; mean_tas_p121_' num2str(i) '];']);
end

% For hurs
hurs_p_121 = [];
for i = [1, 3, 5, 7, 9, 11, 13, 15]
    eval(['hurs_p_121 = [hurs_p_121; mean_hurs_p121_' num2str(i) '];']);
end
%% close to the final product
pr_5=[pr_p_121,pr_p_122,pr_p_123,pr_p_124];
sfcWind_5=[sfcWind_p_121,sfcWind_p_122,sfcWind_p_123,sfcWind_p_124];
tas_5=[tas_p_121,tas_p_122,tas_p_123,tas_p_124];
hurs_5=[hurs_p_121,hurs_p_122,hurs_p_123,hurs_p_124];
%% max analysis
pr_max1=[pr_max,pr_5];
tas_max1=[tas_max,tas_5];
sfcWind_max1=[sfcWind_max,sfcWind_5];
hurs_max1=[hurs_max,hurs_5];

pr_max_final=sum(pr_max1,2)/7
tas_max_final=sum(tas_max1,2)/7
sfcWind_max_final=sum(sfcWind_max1,2)/7
hurs_max_final=sum(hurs_max1,2)/7
%% final step
% my pr_5 is 12960*4 variable. resizedhyde_vector is 12960*1 variable. I want the rows values in the resizedhyde_vector is smaller than 0, correspondingly, the rows inpr_5 will be deleted
pr_5(resizedhyde_vector<=0,:)=[]
sfcWind_5(resizedhyde_vector<=0,:)=[]
tas_5(resizedhyde_vector<=0,:)=[]
hurs_5(resizedhyde_vector<=0,:)=[]
%% my variable
pr_final=[pr,pr_5];
tas_final=[tas,tas_5];
sfcWind_final=[sfcWind,sfcWind_5];
hurs_final=[hurs,hurs_5];

%%
load('resizedhyde_vector_s.mat')
%%

pr_final1=sum(pr_5,2)/7
tas_final1=sum(tas_5,2)/7
sfcWind_final1=sum(sfcWind_5,2)/7
hurs_final1=sum(hurs_5,2)/7

save('./grazingniche/matdata/final_data.mat', 'pr_final','pr_max_final', 'tas_final', 'tas_max_final','sfcWind_final', 'sfcWind_max_final','hurs_final','hurs_max_final','resizedhyde_vector_s');

%% finding the niche!
%% finding the niche!
% Filter out data points where resizedhyde_vector < 1
% The number at the end of these four lines specify which model you want to
% test. In the next steps where you want to find out the 95 interval and
% produce nice presentable graphs from it. You only need to load the
% final_data.mat
% load the final.mat
% this mat file contains pr_final,tas_final,sfcWind_final, and hurs_final,
% and resizedhyde_vector_s
load('final_data.mat') %variables: 'pr_final', 'tas_final', 'sfcWind_final', 'hurs_final','resizedhyde_vector_s');
%pr_final is a 10602*7 variable, each column is a different data source but
%each are the same precipitation
% producing the niche values in the command window


middletas = tas_final1(resizedhyde_vector_s >= 1,1);
middlepr = pr_final1(resizedhyde_vector_s >= 1,1);
middlesfcWind = sfcWind_final1(resizedhyde_vector_s >= 1,1);
middlehurs = hurs_final1(resizedhyde_vector_s >= 1,1);
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

%% Max suggestion: grassland/distribution (Lenton paper Fig1G)
%% final step
% my pr_5 is 12960*4 variable. resizedhyde_vector is 12960*1 variable. I want the rows values in the resizedhyde_vector is smaller than 0, correspondingly, the rows inpr_5 will be deleted
load('grassland_env.mat')
grassland_env_resized = imresize(grassland_env, [90 144]);
imagesc(grassland_env_resized)
colorbar

grassland_env_vector=reshape(grassland_env_resized,[],1)

pr_max_final(grassland_env_vector<=0,:)=[]
tas_max_final(grassland_env_vector<=0,:)=[]
sfcWind_max_final(grassland_env_vector<=0,:)=[]
hurs_max_final(grassland_env_vector<=0,:)=[]

%% scatter图 max
subplot(2,3,1)
X=pr_max_final;
Y=tas_max_final;
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
title({'Grassland Distribution'; 'Precipitation and Temperature'});
ylim([-40,40])
defualtAxes()


subplot(2,3,2)
X=hurs_max_final;
Y=tas_max_final;
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
title({'Grassland Distribution'; 'Humidity and Temperature'});
defualtAxes()
ylim([-40,40])


subplot(2,3,3)
X=sfcWind_max_final;
Y=tas_max_final;
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
title({'Grassland Distribution'; 'Windspeed and Temperature'});
defualtAxes()
ylim([-40,40])


% Save the figure as a SVG and specify the size
set(gcf, 'Position',  [751,163,1092,753])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter_max.svg');
