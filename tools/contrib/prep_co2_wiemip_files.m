% prep_co2_wiemip_files.m
%
% Reminder: MODIFY file_raw and file_wiemip below for each individual case
%
% Workflow:
% - Before running this matlab script, I renamed the TRENDY2025 co2 file (fco2_datm_global_simyr_1700-2024_TRENDY_c250625.nc) to the name of the wiemip co2 file (used below)
% - The script, extends CO2 and time vars in the file from 2024 (TRENDY2025) to 2300 (wiemip)
% - and brings in wiemip co2 to the file
% - After running the matlab script, I used nco to append to each file's history, for example:
%   ncatted -h -a history,global,o,c,"06/25/2025 21:30: converted by TRENDY2024_Data_Prep.ipynb; 06/04/2026: slevis used matlab script tools/contrib/prep_co2_wiemip_files.m to extend the co2 and time variables from 2024 to 2300 and to append the co2 from WIEMIP_hl_co2_ann_2024_2300.txt" WIEMIP_hl_co2_ann_2024_2300_copied_to_fco2_datm_global_simyr_1700-2024_TRENDY_c250625.nc
%
% More information appears in issues
% github.com/ESCOMP/CTSM/issues/4072
% github.com/ESCOMP/CTSM/issues/3936
%
% slevis 2026/06/03

clear

% get trendy CO2, time, time_bnds
file_trendy = 'fco2_datm_global_simyr_1700-2024_TRENDY_c250625.nc';
co2_trendy = ncread(file_trendy, 'CO2');
time_trendy = ncread(file_trendy, 'time');
time_bnds_trendy = ncread(file_trendy, 'time_bnds');

% extend CO2, time, time_bnds
co2_wiemip = co2_trendy;  % orig. to 2024
time_wiemip = time_trendy;  % orig. to 2024
time_bnds_wiemip = time_bnds_trendy;  % orig. to 2024
for yr = 1:276  % out to 2300
  co2_wiemip(:,:,end+1) = co2_wiemip(:,:,end);  % dims (lon, lat, time)
  time_wiemip(end+1) = time_wiemip(end) + 365;  % dims (time)
  time_bnds_wiemip(:,end+1) = time_bnds_wiemip(:,end);  % dims (bnds, time)
end
% Fix time_bnds preexisting glitch in 2021
time_bnds_wiemip(:,322) = time_bnds_wiemip(:,321);  % dims (bnds, time)

% get wiemip co2 for the years 2024-2300
file_raw = '/glade/derecho/scratch/swensosc/WIEMIP/co2/WIEMIP_m_co2_ann_2024_2300.txt';
co2 = readmatrix(file_raw);
co2_wiemip(1,1,325:end) = squeeze(co2(:,2));

% write modified time, time_bnds, and CO2 to the renamed trendy file
file_wiemip = 'WIEMIP_m_co2_ann_2024_2300_copied_to_fco2_datm_global_simyr_1700-2024_TRENDY_c250625.nc';
ncwrite(file_wiemip, 'time', time_wiemip);
ncwrite(file_wiemip, 'time_bnds', time_bnds_wiemip);
ncwrite(file_wiemip, 'CO2', co2_wiemip);

