%% ============================================================
%  UFEMISM CONFIGURATION CREATOR
%
%  Copies a template configuration file and changes the
%  experiment-specific settings.
%
%  The template is never modified.
% =============================================================

clear;
clc;


%% ============================================================
%  USER SETTINGS
% =============================================================

% -------------------------------------------------------------
% Template
% -------------------------------------------------------------

template_file = 'config-files_tetralith/config_R-LIS_init_template.cfg';


% -------------------------------------------------------------
% Experiment settings
%
% Change these when creating a new experiment
% -------------------------------------------------------------

% Ice thickness dataset
ice_thickness = 'H1'; % H1 = Bedmachine, H2 = Bedmap3

switch ice_thickness
    case 'H1'
        ice_thickness_file = ...
    '/home/frare/data/Bedmachine_Antarctica/BedMachineAntarctica_v3_2km.nc';
    case 'H2'
        ice_thickness_file = ...
    '/home/frare/data/Bedmap3_Antarctica/bedmap3_2km.nc';
end

% Ocean dataset
ocean = 'O1';

% SMB dataset
SMB = 'SMB1';

% Geothermal heat flux dataset
GHF = 'GHF1';

% Initialisation method
%
% Examples:
%   'BCAP'
%   'BCAPdH' matches the dH/dt observed 
%   'STD'
%   'STDdH' matches the dH/dt observed
%
initialisation = 'STDdH';


% -------------------------------------------------------------
% Model resolutions of ROI grounding line
% -------------------------------------------------------------

resolutions = [20 10]; %10 5 2];


% -------------------------------------------------------------
% Output directory
% -------------------------------------------------------------

path_to_result = '/nobackup/proj/disk/bolinc/personal/frare/outputs/';

output_directory = 'generated_configs';

%% ============================================================
%  EXPERIMENT NAME
% =============================================================

experiment_name = sprintf( ...
    'R-LIS_init_%s_%s_%s_%s_%s', ...
    initialisation, ice_thickness, ocean, SMB, GHF);


%% ============================================================
%  CHECK TEMPLATE
% =============================================================

if ~isfile(template_file)
    error('Template file not found: %s', template_file);
end


%% ============================================================
%  CREATE OUTPUT DIRECTORY
% =============================================================

if ~isfolder(output_directory)
    mkdir(output_directory);
end


%% ============================================================
%  READ TEMPLATE
% =============================================================

fid = fopen(template_file, 'r');

if fid == -1
    error('Could not open template file: %s', template_file);
end

template = fread(fid, '*char')';

fclose(fid);


%% ============================================================
%  CREATE CONFIGURATION FILES
% =============================================================

fprintf('\nCreating UFEMISM configuration files...\n\n');

for i = 1:length(resolutions)

    resolution = resolutions(i);

    % ---------------------------------------------------------
    % Create filename
    % ---------------------------------------------------------

    output_filename = sprintf( ...
        'config_%s_%dkm.cfg', ...
        experiment_name, resolution);

    result_name = sprintf( ...
        'results_%s_%dkm', ...
        experiment_name, resolution);

    output_file = fullfile( ...
        output_directory, output_filename);


    % ---------------------------------------------------------
    % Start from the original template
    % ---------------------------------------------------------

    config = template;


    % =========================================================
    % CHANGE SETTINGS IN THE CONFIGURATION
    %
    % These are deliberately kept together so you can easily
    % adapt them to your actual UFEMISM configuration.
    % =========================================================

    % ---------------------------------------------------------
    % Write output name in the config file
    % ---------------------------------------------------------
    output_folder = fullfile( ...
        path_to_result,result_name);
    
    config = regexprep(config, ...
    'fixed_output_dir_config\s*=\s*''[^'']*''', ...
    ['fixed_output_dir_config = ''' output_folder '''']);

    % ---------------------------------------------------------
    % Resolution
    % ---------------------------------------------------------
    %
    % Replace the value associated with your resolution
    % parameter in the template.
    %
    % Example:
    if resolution == 20

        config = regexprep(config, ...
    'ROI_maximum_resolution_uniform_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_uniform_config = ''' 50e3 '''']);

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounded_ice_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_grounded_ice_config = ''' 20e3 '''']);

        config = regexprep(config, ...
    'ROI_maximum_resolution_floating_ice_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_floating_ice_config = ''' 20e3 '''']);

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounding_line_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_grounding_line_config = ''' 20e3 '''']);

        config = regexprep(config, ...
    'ROI_maximum_resolution_calving_front_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_calving_front_config = ''' 20e3 '''']);

        config = regexprep(config, ...
    'ROI_maximum_resolution_ice_front_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_ice_front_config = ''' 20e3 '''']);

        config = regexprep(config, ...
    'ROI_maximum_resolution_coastline_config\s*=\s*''[^'']*''', ...
    ['ROI_maximum_resolution_coastline_config = ''' 20e3 '''']);

        config = regexprep(config, ...
    'ROI_grounding_line_width_config\s*=\s*''[^'']*''', ...
    ['ROI_grounding_line_width_config = ''' 10e3 '''']);

        config = regexprep(config, ...
    'ROI_calving_front_width_config\s*=\s*''[^'']*''', ...
    ['ROI_calving_front_width_config = ''' 10e3 '''']);

        config = regexprep(config, ...
    'ROI_ice_front_width_config\s*=\s*''[^'']*''', ...
    ['ROI_ice_front_width_config = ''' 10e3 '''']);

        config = regexprep(config, ...
    'ROI_coastline_width_config\s*=\s*''[^'']*''', ...
    ['ROI_coastline_width_config = ''' 10e3 '''']);

    % ============= 10 km RESOLUTION =========================
    % ========================================================
    elseif resolution == 10
        config = regexprep(config, ...
    'ROI_maximum_resolution_uniform_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_uniform_config = 40e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounded_ice_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_grounded_ice_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_floating_ice_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_floating_ice_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounding_line_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_grounding_line_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_calving_front_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_calving_front_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_ice_front_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_ice_front_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_coastline_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_coastline_config = 10e3 ');

        config = regexprep(config, ...
    'ROI_grounding_line_width_config\s*=\s*[^\r\n]*', ...
    'ROI_grounding_line_width_config = 10e3 ');

        config = regexprep(config, ...
    'ROI_calving_front_width_config\s*=\s*[^\r\n]*', ...
    'ROI_calving_front_width_config = 10e3 ');

        config = regexprep(config, ...
    'ROI_ice_front_width_config\s*=\s*[^\r\n]*', ...
    'ROI_ice_front_width_config = 10e3 ');

        config = regexprep(config, ...
    'ROI_coastline_width_config\s*=\s*[^\r\n]*', ...
    'ROI_coastline_width_config = 10e3');

    % ============= 5 GL km RESOLUTION =========================
    % ==========================================================
    elseif resolution == 5
        config = regexprep(config, ...
    'ROI_maximum_resolution_uniform_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_uniform_config = 30e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounded_ice_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_grounded_ice_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_floating_ice_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_floating_ice_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounding_line_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_grounding_line_config = 5e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_calving_front_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_calving_front_config = 5e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_ice_front_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_ice_front_config = 10e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_coastline_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_coastline_config = 10e3 ');

        % ===== Widths for each resolution band
        config = regexprep(config, ...
    'ROI_grounding_line_width_config\s*=\s*[^\r\n]*', ...
    'ROI_grounding_line_width_config = 8e3 ');

        config = regexprep(config, ...
    'ROI_calving_front_width_config\s*=\s*[^\r\n]*', ...
    'ROI_calving_front_width_config = 8e3 ');

        config = regexprep(config, ...
    'ROI_ice_front_width_config\s*=\s*[^\r\n]*', ...
    'ROI_ice_front_width_config = 8e3 ');

        config = regexprep(config, ...
    'ROI_coastline_width_config\s*=\s*[^\r\n]*', ...
    'ROI_coastline_width_config = 8e3');

    % ============= 2 GL km RESOLUTION =========================
    % ==========================================================
    elseif resolution == 2
        config = regexprep(config, ...
    'ROI_maximum_resolution_uniform_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_uniform_config = 20e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounded_ice_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_grounded_ice_config = 5e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_floating_ice_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_floating_ice_config = 5e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_grounding_line_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_grounding_line_config = 2e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_calving_front_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_calving_front_config = 2e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_ice_front_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_ice_front_config = 5e3');

        config = regexprep(config, ...
    'ROI_maximum_resolution_coastline_config\s*=\s*[^\r\n]*', ...
    'ROI_maximum_resolution_coastline_config = 5e3 ');

        % ===== Widths for each resolution band
        config = regexprep(config, ...
    'ROI_grounding_line_width_config\s*=\s*[^\r\n]*', ...
    'ROI_grounding_line_width_config = 8e3 ');

        config = regexprep(config, ...
    'ROI_calving_front_width_config\s*=\s*[^\r\n]*', ...
    'ROI_calving_front_width_config = 8e3 ');

        config = regexprep(config, ...
    'ROI_ice_front_width_config\s*=\s*[^\r\n]*', ...
    'ROI_ice_front_width_config = 8e3 ');

        config = regexprep(config, ...
    'ROI_coastline_width_config\s*=\s*[^\r\n]*', ...
    'ROI_coastline_width_config = 8e3');

    else % ===== ELSE ERROR
        sprintf('No config for this resolution in Matlab file')
        break
    end

    % ---------------------------------------------------------
    % Ice thickness
    % ---------------------------------------------------------
    %
    config = regexprep(config, ...
    'filename_refgeo_init_ANT_config\s*=\s*''[^'']*''', ...
    ['filename_refgeo_init_ANT_config = ''' ice_thickness_file '''']);

    config = regexprep(config, ...
    'filename_refgeo_PD_ANT_config\s*=\s*''[^'']*''', ...
    ['filename_refgeo_PD_ANT_config = ''' ice_thickness_file '''']);

    config = regexprep(config, ...
    'filename_refgeo_GIAeq_ANT_config\s*=\s*''[^'']*''', ...
    ['filename_refgeo_GIAeq_ANT_config = ''' ice_thickness_file '''']);
    %
    % ---------------------------------------------------------


    % ---------------------------------------------------------
    % SMB
    % ---------------------------------------------------------
    %
    % Example:
    %
    % config = strrep(config, ...
    %     'SMB1', SMB);
    %
    % ---------------------------------------------------------


    % ---------------------------------------------------------
    % Ocean
    % ---------------------------------------------------------
    %
    % Example:
    %
    % config = strrep(config, ...
    %     'O1', ocean);
    %
    % ---------------------------------------------------------


    % ---------------------------------------------------------
    % GHF
    % ---------------------------------------------------------
    %
    % Example:
    %
    % config = strrep(config, ...
    %     'GHF1', GHF);
    %
    % ---------------------------------------------------------

    % ---------------------------------------------------------
    % Bed roughness for higher resolutions that use previous run
    % ---------------------------------------------------------
    if resolution < 20

        % change ending time of simulation end_time_of_run_config
        config = regexprep(config, ...
        'end_time_of_run_config\s*=\s*[^\r\n]*', ...
        'end_time_of_run_config = 10000.0');

        config = regexprep(config, ...
        'protect_grounded_mask_t_end_config\s*=\s*[^\r\n]*', ...
        'protect_grounded_mask_t_end_config = 5000.0');

        config = regexprep(config, ...
        'choice_bed_roughness_config\s*=\s*''[^'']*''', ...
        'choice_bed_roughness_config = ''read_from_file''');

        % Previous resolution
        previous_resolution = resolutions(i-1);

        % Name of the previous simulation
        previous_result_name = sprintf( ...
        'results_%s_%dkm', experiment_name, previous_resolution);

        % Full path to the previous simulation
        previous_bed_roughness_file = fullfile( ...
            path_to_result, ...
            previous_result_name, ...
            'main_output_ANT_LAST.nc');

        % Replace filename in configuration
        config = regexprep(config, ...
        'filename_bed_roughness_ANT_config\s*=\s*''[^'']*''', ...
        ['filename_bed_roughness_ANT_config = ''' ...
        previous_bed_roughness_file '''']);
    end

    % ---------------------------------------------------------
    % Initialisation
    % ---------------------------------------------------------
    %
    switch initialisation
        case 'STD' % standard
            config = regexprep(config, ...
        'choice_bed_roughness_nudging_limit_method_config\s*=\s*''[^'']*''', ...
        'choice_bed_roughness_nudging_limit_method_config = ''uniform''');
        case 'STDdH'
            config = regexprep(config, ...
        'choice_bed_roughness_nudging_limit_method_config\s*=\s*''[^'']*''', ...
        'choice_bed_roughness_nudging_limit_method_config = ''uniform''');

            config = regexprep(config, ...
        'do_target_dHi_dt_config\s*=\s*[^\r\n]*', ...
        'do_target_dHi_dt_config = .TRUE.');
        case 'BCAP'
            config = regexprep(config, ...
        'choice_bed_roughness_nudging_limit_method_config\s*=\s*''[^'']*''', ...
        'choice_bed_roughness_nudging_limit_method_config = ''Martin2011''');
        case 'BCAPdH'
            config = regexprep(config, ...
        'choice_bed_roughness_nudging_limit_method_config\s*=\s*''[^'']*''', ...
        'choice_bed_roughness_nudging_limit_method_config = ''Martin2011''');

            config = regexprep(config, ...
        'do_target_dHi_dt_config\s*=\s*[^\r\n]*', ...
        'do_target_dHi_dt_config = .TRUE.');
    end % END SWITCH INITIALISATION

    %% --------------------------------------------------------
    % WRITE CONFIGURATION FILE
    % ---------------------------------------------------------

    fid = fopen(output_file, 'w');

    if fid == -1
        error('Could not create file: %s', output_file);
    end

    fprintf(fid, '%s', config);

    fclose(fid);


    fprintf('Created: %s\n', output_file);

end


%% ============================================================
%  FINISHED
% =============================================================

fprintf('\n---------------------------------------------\n');
fprintf('Finished creating configuration files.\n');
fprintf('Output directory: %s\n', output_directory);
fprintf('Experiment: %s\n', experiment_name);
fprintf('---------------------------------------------\n\n');