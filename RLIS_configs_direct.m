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

template_file = 'config-files_tetralith/config_R-LIS_init_template_direct.cfg';


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
% SMB1 for RACMO
% SMB2 for HIRHAM5
SMB = 'SMB2';

switch SMB
    case 'SMB1'
        SMB_file1 = ...
    '/home/frare/data/RACMO/Antarctica/RACMO_SMB_1979_2016.nc';
        SMB_file2 = ...
    '/home/frare/data/climate/RACMO_clim_1979_2021_winds_noNaN.nc';
    case 'SMB2'
        SMB_file1 = ...
    '/home/frare/data/SMB/HIRHAM5_SMB_1979_2021.nc';
        SMB_file2 = ...
    '/home/frare/data/climate/HIRHAM5_clim_1979_2021.nc';
end

% Geothermal heat flux dataset
% GHF1 for Shapiro 2004
% GHF2 for LoesingEbbing 2021
% GHF3 for Staal 2020 - not yet implemented
GHF = 'GHF1';

switch GHF
    case 'GHF1'
        GHF_file = ...
    '/home/frare/data/GHF/ShapiroRitzwoller2004_global.nc';
    case 'GHF2'
        GHF_file = ...
    '/home/frare/data/GHF/LoesingEbbing_2021_topocorr_4km.nc';
end

% Initialisation method
%
% Examples:
%   'BCAP'
%   'BCAPdH' matches the dH/dt observed 
%   'STD'
%   'STDdH' matches the dH/dt observed
%
initialisation = 'BCAP';

sliding = 'SL1';
% SL1: 'ZI' Zoet-Iverson used in the template
% SL2: 'Budd' 
% SL3: 'Weertman'

% -------------------------------------------------------------
% Model resolutions of ROI grounding line
% -------------------------------------------------------------

resolutions = 2; % 2 km

% -------------------------------------------------------------
% Output directory
% -------------------------------------------------------------

path_to_result = '/nobackup/proj/disk/bolinc/personal/frare/outputs/';

output_directory = 'generated_configs';

%% ============================================================
%  EXPERIMENT NAME
% =============================================================

experiment_name = sprintf( ...
    'direct_R-LIS_init_%s_%s_%s_%s_%s_%s', ...
    initialisation, ice_thickness, ocean, SMB, GHF, sliding);


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
    % ============= 2 GL km RESOLUTION =========================
    % ==========================================================
    if resolution == 2
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
    config = regexprep(config, ...
    'filename_SMB_prescribed_ANT_config\s*=\s*''[^'']*''', ...
    ['filename_SMB_prescribed_ANT_config = ''' SMB_file1 '''']);

    config = regexprep(config, ...
    'filename_climate_snapshot_ANT_config\s*=\s*''[^'']*''', ...
    ['filename_climate_snapshot_ANT_config = ''' SMB_file2 '''']);
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
   
    config = regexprep(config, ...
    'filename_geothermal_heat_flux_config\s*=\s*''[^'']*''', ...
    ['filename_geothermal_heat_flux_config = ''' GHF_file '''']);
    %
    % ---------------------------------------------------------

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

            config = regexprep(config, ...
        'target_dHi_dt_t_end_config\s*=\s*[^\r\n]*', ...
        'target_dHi_dt_t_end_config = 30000.0');
            
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

            config = regexprep(config, ...
        'target_dHi_dt_t_end_config\s*=\s*[^\r\n]*', ...
        'target_dHi_dt_t_end_config = 30000.0');
    end % END SWITCH INITIALISATION

    % ---------------------------------------------------------
    % sliding law
    % ---------------------------------------------------------
    switch sliding
        case 'SL1' % ZI
            % do nothing, same as template
        case 'SL2' % Budd
            config = regexprep(config, ...
        'choice_sliding_law_config\s*=\s*''[^'']*''', ...
        'choice_sliding_law_config = ''Budd''');
        case 'SL3' % Wertman
            config = regexprep(config, ...
        'choice_sliding_law_config\s*=\s*''[^'']*''', ...
        'choice_sliding_law_config = ''Weertman''');
    end


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