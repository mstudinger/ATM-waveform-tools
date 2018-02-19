function [poly_lon, poly_lat, t_start, t_end, data_fmt] = atm_wvfm_info(f_name_inp, verbose)
%
% SYNTAX: [poly_lon, poly_lat, t_start, t_end, data_fmt] = ATM_WVFM_INFO(f_name_inp, verbose)
%
%                   ATM_WVFM_INFO reads and optionally displays some basic information about the contents
%                   of an ATM HDF5 file with laser waveforms. 
%                   The output can be used as input for spatial and temporal searches that import all waveforms of a file.
%
% INPUT PARAMETERS:
%   f_name_inp      File name of HDF5 file to be read including path.
%   verbose         1 displays basic information on the console. 0 is silent.
%
% OUTPUT PARAMETERS:
%   poly_lon        longitude vertices of min and max longitude polygon
%   poly_lat        latitude  vertices of min and max latitude  polygon
%   t_start         time of first laser shot in datestr format 31
%   t_end           time of last  laser shot in datestr format 31
%   data_fmt        type of data format:
%                   0: cannot determine data format type from HDF5 file
%                   1: legacy  ATM waveform format stored in subgroup /waveforms/vld
%                   2: current ATM waveform format stored in subgroup /waveforms/twv
%
% Michael Studinger, NASA Goddard Space Flight Center, Greenbelt, MD.
% Last changed:     02/14/2018
% Version:          3.01
% See also:         atm_wvfm_reader
%

%% basic error checking

if (nargin ~= 2)
    error('atm_wvfm_info:nargin', ['\n\tERROR: Number of input arguments must be 2:\n', ...
        '\tUSAGE: [poly_lon, poly_lat, t_start, t_end, data_fmt_type] = atm_wvfm_info_v1(f_name_inp, verbose);'])
end

if (nargout ~= 5)
    error('atm_wvfm_info:nargout', ['\n\tERROR: Number of out arguments must be 4:\n', ...
        '\tUSAGE: [poly_lon, poly_lat, t_start, t_end, data_fmt_type] = atm_wvfm_info_v1(f_name_inp, verbose);'])
end

if (exist(f_name_inp) == 0)
    warndlg({'Input file not found!';'Script aborted.'},'!! Warning !!');
    error('atm_wvfm_info:file_chk', '\n\tInput file not found. Script aborted.')
end

if (verbose ~= 0 && verbose ~= 1)
    error('atm_wvfm_info:verbose', '\n\tERROR: verbose must be either 0 or 1.')
end

%% read file

% time of first and last shot
first_shot_str = h5read(f_name_inp,'/ancillary_data/time/first_shot'); % 2017-05-10T13:26:28-0000 in ISO-8601 standard
last_shot_str = h5read(f_name_inp,'/ancillary_data/time/last_shot');

first_shot_numdate = datenum(first_shot_str,'yyyy-mm-ddTHH:MM:SS-0000'); % just in case fields are duplicated. 
last_shot_numdate = datenum(last_shot_str,'yyyy-mm-ddTHH:MM:SS-0000');

t_start = datestr(first_shot_numdate,31);
t_end = datestr(last_shot_numdate + 1/(60*60*24),31); % need to add 1 seconds to make sure all laser shots get selected

% get latitude and longitude boundaries
min_lon = wrapTo180(h5read(f_name_inp,'/ancillary_data/meta_data/min_longitude'));
max_lon = wrapTo180(h5read(f_name_inp,'/ancillary_data/meta_data/max_longitude'));

min_lat = h5read(f_name_inp,'/ancillary_data/meta_data/min_latitude');
max_lat = h5read(f_name_inp,'/ancillary_data/meta_data/max_latitude');

% make closed polygon for spatial search
% arrange the vertices in clockwise order so it doesn't trigger a warning in atm_wvfm_reader.m
poly_lon = [min_lon max_lon max_lon min_lon min_lon];
poly_lat = [max_lat max_lat min_lat min_lat max_lat];

% determine the version of ATM waveform data format
info = h5info(f_name_inp);

if (size(info.Groups,1) == 7)
    if (isfield(info,'Groups') && strcmp(info.Groups(7).Groups(1).Name,'/waveforms/twv'))
        data_fmt = 2; % current ATM waveform format stored in subgroup /waveforms/twv
        data_fmt_str = 'ATM Waveform Format Vers. 2.0';
    elseif (isfield(info,'Groups') && strcmp(info.Groups(8).Groups.Name,'/waveforms/vld'))
        data_fmt = 1; % legacy ATM waveform format stored in subgroup /waveforms/vld
        data_fmt_str = 'ATM Waveform Format Vers. 1.0';
        error('atm_wvfm_reader:wvfm_fmt_vers', '\n\tERROR: ATM waveform format version 1 not yet implemented.')
    else
        error('atm_wvfm_reader:wvfm_fmt_vers', '\n\tERROR: Cannot determine ATM waveform format version from HDF5 file.')
    end
else
    error('atm_wvfm_info:input_struct','\n\tERROR: ATM waveform struct must have 7 groups.\n')
end

if (verbose == 1)
    
    [~,name,ext] = fileparts(f_name_inp); 
    
    ele = h5read(f_name_inp,'/footprint/elevation');

    % get sampling interval
    sample_int_ns = h5read(f_name_inp,'/waveforms/twv/ancillary_data/sample_interval'); % sampling interval in nano seconds
    % sample_int = sample_int_ns*1E-9; % sampling interval in seconds


    fprintf('\n-----------------------------------------------------------------------\n');
    fprintf('Contents of file %s\n',[name ext]);
    fprintf('-----------------------------------------------------------------------\n');
    fprintf('UTC start: %s  |  UTC end:\t%s\n',datestr(first_shot_numdate,'yyyy/mm/dd HH:MM:SS.FFF'),datestr(last_shot_numdate,'yyyy/mm/dd HH:MM:SS.FFF'));
    fprintf('Length: %s                    |  N shots: %d\n',datestr(last_shot_numdate - first_shot_numdate,'HH:MM:SS'),size(ele,1));
    fprintf('------------------------------------|----------------------------------\n');
    fprintf('Sampling interval:    %4.2f ns       |  Data format type: %d \n',sample_int_ns,data_fmt);
    fprintf('------------------------------------|----------------------------------\n');
    fprintf('Longitude minimum: %9.4f°\n',min_lon);
    fprintf('Longitude maximum: %9.4f°\n',max_lon);
    fprintf('Latitude  minimum: %9.4f°\n',min_lat);
    fprintf('Latitude  maximum: %9.4f°\n',max_lat);
    fprintf('-----------------------------------------------------------------------\n');
    
end

end

