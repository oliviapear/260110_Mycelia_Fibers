% Olivia Pear
% Tensile test code v4

% Last updated: 08-11-2025

%% INSTRUCTIONS FOR USE

% 1. Remove initial indexing column when copying Exponent Raw Data to excel
% 2. Import data as output type: cell array 
% 3. Rename imported raw data table as 'Tensile1' 
% 4. Import fiber diameters as output type: table and rename 'diameters'. Format this as a single column corresponding to sample and
% replicate order
% 5. Young's modulus is currently calculated by linear fit between 1-5%
% strain for hydrated, 0.5-1% strain for dry


% get the entire first row of the cellarray, that cellarray row is
% all the samples and their replicates eg (Name_#)
row_one = Tensile1(1,:);
sample_name = '';
sample_names = [];
sample_replicates = [];

array_index = 0;
replicate_count = 0;
cols_per_replicate = 4;

for i = row_one    
    name_expression = '(.*)_';
    tmp_sample_name = regexp(string(i),name_expression,'match');
    
    % tmp_sample_name = extractBetween(string(i),1,strlength(string(i))-1);

    if sample_name ~= tmp_sample_name
        replicate_count = 0;
        array_index = array_index + 1;
        sample_name = tmp_sample_name;
        sample_names = [sample_names, tmp_sample_name];
    end
    
    replicate_count = replicate_count + 1;
    sample_replicates(array_index) = replicate_count/cols_per_replicate;
end

corrected_sample_names = [];
for i = 1:length(sample_names)
    for j = 1:sample_replicates(i)
        corrected_sample_names = [corrected_sample_names, strcat(sample_names(i),'_',num2str(j))];
    end
end


%%Geometries
% diameter[mm]
% cross section area = pi*(diameter/2)^2;
total_replicates = sum(sample_replicates);
raw_data_index = 1;
area = (pi/(10^6))*(diameters{:, 1}./2).^2; % [m^2]

gage_length = 5; %mm
force_index = 1;
distance_index = 2;
increment = 4;
sz = 1;


% create cell array for corrected data
corrected_stress_strain_cell_array = cell(length(Tensile1(:,1)), total_replicates*2);

%create cell array for all replicates
all_replicates_summary_cell_array = cell(total_replicates, 4);

% store summary data per sample in order to average
% total_replicates
for i = 1:total_replicates
    strain_linear = [];
    stress_linear = [];
    raw_distance = Tensile1(:,distance_index);
    raw_force = Tensile1(:,force_index);
    distance = [];
    force = [];

    % remove strings from data
    for j = 1:length(raw_distance)
        % only use numbers for distance and force arrays
        if isnumeric(raw_distance{j})
            distance = [distance;raw_distance{j}];
            force = [force;raw_force{j}];
        end
    end

    %convert to strain[mm/mm] and stress[Pa]
    strain = distance/gage_length;  
    stress = force*9.8/(1000*area(i)); %i is replicate index

    

    corrected_strain = [];
    corrected_stress = [];
    % correct data
    
    % ~ismember(stress(j), corrected_stress) &
    for j = 1:length(stress)
        % max_index = find(stress==max(stress));
        % if stress(j) >=0 & j <= max_index
        corrected_strain = [corrected_strain;strain(j)];
        corrected_stress = [corrected_stress;stress(j)];
        % end
    end
    
    odds = 2*i-1;
    evens = 2*i;

    ultimate_stress = max(corrected_stress);
   
    toughness = trapz(strain, stress);

   %% Young's modulus - linear fit
    % Identify the linear portion (based on data)
 
    for j = 1:length(corrected_strain)
        if corrected_strain(j) >= 0.1 && corrected_strain(j) <= 0.15
            strain_linear = [strain_linear;corrected_strain(j)];
            stress_linear = [stress_linear;corrected_stress(j)];
        end
    end
    % Extract the x and y values for the linear portion   
    % strain_linear = corrected_strain(linear_start:linear_end);
    % stress_linear = corrected_stress(linear_start:linear_end);
    % Perform linear regression
    coefficients = polyfit(strain_linear, stress_linear, 1); % Fit a first-degree polynomial (line)
    slope = coefficients(1); % Extract the slope
    
    %plot(strain_linear,stress_linear, 'linestyle' , '--')
    %hold on
    %xlabel('Linear Strain (%)')
    %ylabel('Linear Stress (Pa)')
    %legend


    all_replicates_summary_cell_array{i,1} = corrected_sample_names(i);
    all_replicates_summary_cell_array{i,2} = ultimate_stress;
    all_replicates_summary_cell_array{i,3} = toughness;
    all_replicates_summary_cell_array{i,4} = slope;  
    
    %set all data in cell array
    corrected_stress_strain_cell_array(3:length(corrected_stress)+2,evens) = num2cell(corrected_stress);
    corrected_stress_strain_cell_array(1,evens) = cellstr(corrected_sample_names(i));
    corrected_stress_strain_cell_array(2,evens) = cellstr('Stress (Pa)');

    corrected_stress_strain_cell_array(3:length(corrected_strain)+2,odds) = num2cell(corrected_strain*100);
    corrected_stress_strain_cell_array(1,odds) = cellstr(corrected_sample_names(i));
    corrected_stress_strain_cell_array(2,odds) = cellstr('Strain (%)');
    
 
   
    scatter(strain*100, stress)
    hold on
    xlabel('Strain (%)')
    ylabel('Stress (Pa)')
    legend(corrected_sample_names)

    distance_index = distance_index + increment;
    force_index = force_index + increment;
end

%convert cellarray to table
all_replicates_summary_table = cell2table(all_replicates_summary_cell_array, "VariableNames",["Sample Name","Ultimate Stress (Pa)","Toughness(Pa)","Youngs modulus (Pa)"]);

% find averages
avg_std_cell_array = cell(length(sample_names), 7);

for i = 1:length(sample_names)
    ult_array = [];
    tough_array = [];
    youngs_array = [];
    avg_std_cell_array{i,1} = sample_names(i);
    for j = 1:length(all_replicates_summary_cell_array(:,1))
        if contains(all_replicates_summary_cell_array{j, 1}, sample_names(i))
            ult_array = [ult_array; all_replicates_summary_cell_array{j, 2}];
            tough_array = [tough_array; all_replicates_summary_cell_array{j, 3}];
            youngs_array = [youngs_array; all_replicates_summary_cell_array{j, 4}];
        end
    end
    avg_std_cell_array{i,2} = mean(ult_array);
    avg_std_cell_array{i,3} = mean(tough_array);
    avg_std_cell_array{i,4} = mean(youngs_array);
    avg_std_cell_array{i,5} = std(ult_array);
    avg_std_cell_array{i,6} = std(tough_array);
    avg_std_cell_array{i,7} = std(youngs_array);
end

avg_std_summary_table = cell2table( ...
    avg_std_cell_array, ...
    "VariableNames", ...
    [ ...
        "Sample Name", ...
        "Average Ultimate Stress (Pa)", ...
        "Average Toughness(Pa)", ...
        "Average Youngs modulus (Pa)", ...
        "StDev Ultimate Stress (Pa)", ...
        "StDev Toughness(Pa)", ...
        "StDev Youngs modulus (Pa)" ...
    ] ...
);

% Specify the file name and export as csv
% writetable(all_replicates_summary_table,'Summary_Table.csv')
filename = "Summary_table_"+string(datetime('now','Format','yyyy-MM-dd_HH_mm_ss'))+".xlsx";
writecell(corrected_stress_strain_cell_array,filename,'Sheet','Raw Data');
writetable(all_replicates_summary_table,filename,'Sheet','Analyzed Data','WriteVariableNames',true);
writetable(avg_std_summary_table,filename,'Sheet','Statistics','WriteVariableNames',true);