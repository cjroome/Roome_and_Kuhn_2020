tic
clear all
close all
load('cell_list.mat')
for cell_list_number = 8:size(cell_list,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load imaging raw data (dendrites) cell by cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell_number = cell_list(cell_list_number,1);
cell_number
trace_length = 100000;
start_trace = 1;
control_only = 1;
if control_only == 1
traces = cell_list(cell_list_number,2);
else
traces = cell_list(cell_list_number,5);
end
start_trace_control = 1;
traces_control = cell_list(cell_list_number,2);  
lines =  cell_list(cell_list_number,3);
start_trace_drug = cell_list(cell_list_number,4);
traces_drug = cell_list(cell_list_number,5);
start_line = cell_list(cell_list_number,6);
end_line = cell_list(cell_list_number,7);
s1 = 'cell';
s2 = num2str(cell_number);
s3 = '.xlsx';
s4 = '_';
filename = strcat(s1,s2,s3);
cell = strcat(s1,s2,s4);
interpol = 5;
segments = 4;
edges = 25;
seg = floor((lines - 2*edges)/segments);
spectral_mixing_ratio_a = 0.3;
spectral_mixing_ratio_b = 1.5;
start_time = 4991;
time_window = 100000;
end_time = time_window + 4990;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral mixing correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1
V_data = zeros(lines,time_window,traces);
C_data = zeros(lines,time_window,traces);
for i = start_trace:traces
tracenum = num2str(i);
V_data_raw(:,:) = im2double(imread(strcat(cell,tracenum,'.tif'),1));
C_data_raw(:,:) = im2double(imread(strcat(cell,tracenum,'.tif'),2));
V_data_inter(:,:) = kron(V_data_raw(:,:),ones(1,interpol));
C_data_inter(:,:) = kron(C_data_raw(:,:),ones(1,interpol));
V_data_raw_in(:,:) = V_data_inter(:,start_time:end_time);
C_data_raw_in(:,:) = C_data_inter(:,start_time:end_time);
V_data(:,:,i) = V_data_raw_in;
C_data(:,:,i) = C_data_raw_in;
end
end
clear V_data_raw C_data_raw V_data_inter C_data_inter V_data_raw_in C_data_raw_in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movement correction raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1
    mov_xcorr_shift = zeros(1,size(V_data,2));
for i = start_trace:traces
   ref_profile = mean(V_data(:,trace_length,i),2);
   for tt = 1:size(V_data,2)
   mov_xcorr = xcorr(ref_profile,V_data(:,tt,i));
   [max_val, pos] = max(mov_xcorr);
   mov_xcorr_shift(1,tt) = pos - (lines);
   V_data(:,tt,i) = circshift(V_data(:,tt,i),mov_xcorr_shift(1,tt),1);
   C_data(:,tt,i) = circshift(C_data(:,tt,i),mov_xcorr_shift(1,tt),1);
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CS detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1
DCSbinary = zeros(size(V_data,2),traces);
DCSbinary_corrected = zeros(size(V_data,2),traces);
locs_for_excel = ones(50,traces)*105000;
for i = start_trace_control:traces
    filtered_V_data = movmean(mean(V_data(:,:,i),1),50) - movmean(mean(V_data(:,:,i),1),500);
    filtered_C_data = (movmean(mean(C_data(:,:,i),1),50) - mean(mean(C_data(:,:,i),1),2))./mean(mean(C_data(:,:,i),1),2);
    filtered_V_data_fine = movmean(mean(V_data(:,:,i),1),5) - movmean(mean(V_data(:,:,i),1),500);
    V_CS_detection_threshold = 3*std(filtered_V_data);
    V_CS_detection_threshold_fine = 3*std(filtered_V_data_fine);
    [~, locs] = findpeaks(filtered_V_data,'MinPeakHeight',V_CS_detection_threshold,'MinPeakDistance',200);
    locs(locs<51) = [];
    DCSbinary(locs-50,i) = 1;
    for tt = 1+200:trace_length-400
    if DCSbinary(tt,i) == 1 && sum(DCSbinary(tt-199:tt+200,i),1) == 1
        C_CS_section = filtered_C_data(1,tt-199:tt+400) - mean(filtered_C_data(1,tt-199:tt),2);
        if max(C_CS_section) >= 0.02
    V_CS_section = filtered_V_data_fine(tt-199:tt+200);
    binary_CS_section = DCSbinary(tt-199:tt+200 ,i);
    [pks1, locs1] = findpeaks(V_CS_section,'MinPeakHeight',V_CS_detection_threshold_fine,'MinPeakDistance',10);
    if size(locs1,2) == 0
    locs1(1,1) = 0; 
    end
    binary_shift = 200-locs1(1,1) + 20;
    DCSbinary_corrected(tt-binary_shift,i) = 1;
        end
    end
    end
    DCSbinary(:,i) = DCSbinary_corrected(:,i);
    [pks, locs] = findpeaks(DCSbinary(:,i),'MinPeakHeight',0.9,'MinPeakDistance',10);
    locs_for_excel(1:size(locs,1),i) = locs;
  
end
%figure(1)
%plot(movmean(mean(V_data(:,1:50000,i),1),1) - movmean(mean(V_data(:,1:50000,i),1),500));
%hold on
%plot(DCSbinary(1:50000,i)./100);
%plot(DCSbinary(1:size(V_data,2),i));
%ylim([-0.01 0.01]);
%figure(2)
%plot(mean(C_data(:,1:50000,i),1))
%hold on  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ephys data and binary traces (soma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cell_number < 8
E_data = xlsread(filename,1);
else
go = 0; % for auto DCS detection use go = 0   
E_data = xlsread(filename,1);
SSbinary = xlsread(filename,2);
DCSpositions = xlsread(filename,3);
DCSpositions(isnan(DCSpositions)) = 105000;
DCSpositions = int64(DCSpositions);
DSpositions = xlsread(filename,4);
DSpositions(isnan(DSpositions)) = 105000;
DSpositions = int64(DSpositions);
if go == 1
DCSbinary = zeros(105000,traces);
DSbinary = zeros(105000,traces);
end
windowsize = 500;
SSsum = 0;
DCSevents = 0;
DCSevents_traces = zeros(1,traces);
DSevents = 0;
SSbinary_cut = SSbinary;
SSbinary_shifted = SSbinary;
DCSbinary_cut = ones(lines,100000,traces);
DCSbinary_full = zeros(lines,100000,traces);
for i = start_trace:traces
if go == 1
DCSbinary(DCSpositions(:,i),i) = 1;
DSbinary(DSpositions(:,i),i) = 1;
else
DSbinary  = DCSbinary;
end
for K = windowsize:100000-(4*windowsize)
  if DCSbinary(K,i) == 1
  DCSevents = DCSevents + 1;
  DCSevents_traces(1,i) = DCSevents_traces(1,i) + 1;
  SSbinary_cut(K-99:K+400,i) = 0;
  DCSbinary_cut(:,K-99:K+400,i) = 0;
  DCSbinary_full(:,K-99:K+400,i) = 10000;
  end
  if DSbinary(K,i) == 1
  DSevents = DSevents + 1;
  SSbinary_cut(K-99:K+400,i) = 0;
  DCSbinary_cut(:,K-99:K+400,i) = 0;
  DCSbinary_full(:,K-99:K+400,i) = 10000;
  end
end
end
end
DCSbinary = DCSbinary(1:100000,:);
DSbinary = DSbinary(1:100000,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% behaviour data input if it exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cell_number >= 8  
interpol_behaviour = 100;
Eye_puff = xlsread(filename,5);
Eye_no_puff = xlsread(filename,6);
Eye_puff = Eye_puff(50:end,:);
Eye_no_puff = Eye_no_puff(50:end,:);
Eye_puff_diff = diff(-(Eye_puff));
Eye_no_puff_diff = diff(-(Eye_no_puff));
Eye_puff = Eye_puff';
Eye_no_puff = Eye_no_puff';
Eye_puff_diff = Eye_puff_diff';
Eye_no_puff_diff = Eye_no_puff_diff';
Eye_puff_int = kron(Eye_puff,ones(1,interpol_behaviour));
Eye_puff_diff_int = kron(Eye_puff_diff,ones(1,interpol_behaviour));
Eye_no_puff_int = kron(Eye_no_puff,ones(1,interpol_behaviour));
Eye_no_puff_diff_int = kron(Eye_no_puff_diff,ones(1,interpol_behaviour));
Eye_puff_int = Eye_puff_int(start_trace:traces_control,1:trace_length);
Eye_puff_diff_int = Eye_puff_diff_int(start_trace:traces_control,1:trace_length);
Eye_no_puff_int = Eye_no_puff_int(start_trace:traces_control,1:trace_length);
Eye_no_puff_diff_int = Eye_no_puff_diff_int(start_trace:traces_control,1:trace_length);
Eye_closed_min = min(Eye_puff_int(Eye_puff_int~= 0));
if isempty(Eye_closed_min)
    Eye_closed_min = 0;
end
Eye_open_max = max(Eye_puff_int(:,50000));
no_Eye_closed_min = min(Eye_no_puff_int(Eye_puff_int~= 0));
if isempty(no_Eye_closed_min)
    no_Eye_closed_min = 0;
end
no_Eye_open_max = max(Eye_no_puff_int(:,50000));
for i = start_trace:traces_control 
    Eye_puff_int(i,:) = (Eye_puff_int(i,:) - Eye_closed_min)./(Eye_open_max - Eye_closed_min);
    Eye_no_puff_int(i,:) = (Eye_no_puff_int(i,:) - no_Eye_closed_min)./(no_Eye_open_max - no_Eye_closed_min);
end
Eye_puff = xlsread(filename,7);
Eye_no_puff = xlsread(filename,8);
Eye_puff = Eye_puff(50:end,:);
Eye_no_puff = Eye_no_puff(50:end,:);
Eye_puff_diff = diff(-(Eye_puff));
Eye_no_puff_diff = diff(-(Eye_no_puff));
Eye_puff = Eye_puff';
Eye_no_puff = Eye_no_puff';
Eye_puff_diff = Eye_puff_diff';
Eye_no_puff_diff = Eye_no_puff_diff';
Eye_puff_int_drug = kron(Eye_puff,ones(1,interpol_behaviour));
Eye_puff_diff_int_drug = kron(Eye_puff_diff,ones(1,interpol_behaviour));
Eye_no_puff_int_drug = kron(Eye_no_puff,ones(1,interpol_behaviour));
Eye_no_puff_diff_int_drug = kron(Eye_no_puff_diff,ones(1,interpol_behaviour));
Eye_puff_int_drug = Eye_puff_int_drug(1:(traces_drug-traces_control),1:trace_length);
Eye_puff_diff_int_drug = Eye_puff_diff_int_drug(1:(traces_drug-traces_control),1:trace_length);
Eye_no_puff_int_drug = Eye_no_puff_int_drug(1:(traces_drug-traces_control),1:trace_length);
Eye_no_puff_diff_int_drug = Eye_no_puff_diff_int_drug(1:(traces_drug-traces_control),1:trace_length);
for i = start_trace:(traces_drug-traces_control)
    Eye_puff_int_drug(i,:) = (Eye_puff_int_drug(i,:) - Eye_closed_min)./(Eye_open_max - Eye_closed_min);
    Eye_no_puff_int_drug(i,:) = (Eye_no_puff_int_drug(i,:) - no_Eye_closed_min)./(no_Eye_open_max - no_Eye_closed_min);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatio-temporal averaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_data_space_av = mean(V_data(:,1:trace_length,:),1);
C_data_space_av = mean(C_data(:,1:trace_length,:),1);
V_data_space_temp_av = movmean(V_data_space_av(:,1:trace_length,:),100,2);
C_data_space_temp_av = movmean(C_data_space_av(:,1:trace_length,:),100,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hot pixel removal part1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1
V_data_linescan = V_data;
else
V_data_linescan = movmean((movmean(V_data,10,1)),50,2);
end
C_data_linescan = movmean((movmean(C_data,10,1)),50,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use DCS to calulate V C baselines and calcium signal crosstalk scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1
DCSeventcounter = 1;
beta_out_traces = zeros(1,traces_control);
DCSbinary_prestim = zeros(trace_length,traces_control);
DCSbinary_prestim(:,:) = DCSbinary(:,1:traces_control);
for i = start_trace_control:traces_control
  V_DFF_test_in = mean(V_data(:,1:(trace_length/2),i),1);
  V_DFF_test_in = (V_DFF_test_in - movmean(V_DFF_test_in,10000))./mean(V_DFF_test_in);
  V_DFF_test_in_std = std(V_DFF_test_in);
    if V_DFF_test_in_std < 0.03
     for K = windowsize:(trace_length/2)-(4*windowsize)
        if DCSbinary_prestim(K,i) == 1 && sum(DCSbinary_prestim(K:K+4*windowsize,i),1) == 1
         DCS_V = V_data(:,K-(windowsize-1):(K+4*windowsize),i);
         DCS_C = C_data(:,K-(windowsize-1):(K+4*windowsize),i);
         DCS_V_single = mean(DCS_V,1);
         DCS_V_single(1,windowsize:windowsize+200) = mean(DCS_V_single(1,1:windowsize),2);
         DCS_V_single = movmean(DCS_V_single,200);
         DCS_C_single = mean(DCS_C,1);
         DCS_C_single = movmean(DCS_C_single,200);
         DCS_C_single_DFF = (DCS_C_single - mean(DCS_C_single(1,1:windowsize),2))./mean(DCS_C_single(1,1:windowsize),2);
         DCS_V_single_DFF = (DCS_V_single - mean(DCS_V_single(1,1:windowsize),2))./mean(DCS_V_single(1,1:windowsize),2);
            if max(DCS_C_single_DFF) > 0.05 && mean(DCS_V_single_DFF(1,2*windowsize:3*windowsize),2) > 0.01
             y1 = DCS_V_single(1:windowsize)';
             y2 = DCS_V_single(windowsize+100:2*windowsize)';
             y3 = [y1;y2];
             y = y3(:);
             x1 = DCS_C_single(1:windowsize)';
             x2 = DCS_C_single(windowsize+100:2*windowsize)';
             x3 = [x1;x2];
             x = x3(:);
             X = [ones(size(x)) x];
             beta = regress(y,X);
             beta_out(DCSeventcounter) = beta(2);
             if beta_out(DCSeventcounter) > 1
                beta_out(DCSeventcounter) = 0;
             end
             DCSeventcounter = DCSeventcounter+1;
             else
             beta_out(DCSeventcounter) = 0;
             DCSeventcounter = DCSeventcounter+1;
            end     
        end
     end
   end
if exist('beta_out') == 1
beta_out_traces(1,i) = mean(beta_out(beta_out~=0));
if isnan(beta_out_traces(1,i))
beta_out_traces(1,i) = 0;
end
clear beta_out
end
end
Beta_all_mean = mean(beta_out_traces(beta_out_traces~=0));
if Beta_all_mean < 0 
Beta_all_mean = 0;
end
else
Beta_all_mean = 0;
end
%Beta_all_mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% or use e-stim to calulate V C baselines and calcium signal crosstalk scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 0;
if go == 1
ESeventcounter = 1;
beta_out_traces = zeros(1,5);
for i = 1:5
for K = windowsize:trace_length-(4*windowsize)
  if K == 50000
     DCS_V = V_data(:,K-(windowsize-1):(K+4*windowsize),i);
     DCS_C = C_data(:,K-(windowsize-1):(K+4*windowsize),i);
     DCS_V_single = mean(DCS_V,1);
     DCS_C_single = mean(DCS_C,1);
     y1 = DCS_V_single(1:windowsize)';
     y2 = DCS_V_single(windowsize+100:2*windowsize)';
     y3 = [y1;y2];
     y = y3(:);
     x1 = DCS_C_single(1:windowsize)';
     x2 = DCS_C_single(windowsize+100:2*windowsize)';
     x3 = [x1;x2];
     x = x3(:);
     X = [ones(size(x)) x];
     beta = regress(y,X);
     beta_out(ESeventcounter) = beta(2);
     DCSeventcounter = ESeventcounter+1;
  end
end
beta_out_traces(1,i) = mean(beta_out);
end
Beta_all_mean = mean(beta_out_traces);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All DF/F calculations with movement and calcium crosstalk removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear V_data C_data 
V_DFF_space_av = zeros(1,trace_length,traces);
C_DFF_space_av = zeros(1,trace_length,traces);
for i = start_trace:traces
C_data_crosstalk = (C_data_space_av(1,:,i)*Beta_all_mean);
V_data_space_av(:,:,i) = (V_data_space_av(:,:,i)-C_data_crosstalk(:,:));
V_DFF_space_av(:,:,i) = (V_data_space_av(:,:,i)-mean2(V_data_space_av(:,:,i)))./mean2(V_data_space_av(:,:,i));
V_DFF_movements = movmean(V_DFF_space_av(:,:,i),25000,2) - mean(V_DFF_space_av(:,:,i));
V_DFF_space_av(:,:,i) = (V_DFF_space_av(:,:,i)-V_DFF_movements(:,:));
C_DFF_space_av(:,:,i) = (C_data_space_av(:,:,i)-mean2(C_data_space_av(:,:,i)))./mean2(C_data_space_av(:,:,i));
C_DFF_space_av(:,:,i) = (C_DFF_space_av(:,:,i)-V_DFF_movements(:,:));
end
V_DFF_linescan = zeros(lines,trace_length,traces);
V_data_linescan_corr = zeros(lines,trace_length,traces);
C_DFF_linescan = zeros(lines,trace_length,traces);
for i = start_trace:traces
for ll = 1:lines
C_data_crosstalk = (C_data_linescan(ll,1:trace_length,i)*Beta_all_mean);
V_data_linescan_corr(ll,:,i) = (V_data_linescan(ll,1:trace_length,i)-C_data_crosstalk(:,1:trace_length));
V_DFF_linescan(ll,:,i) = (V_data_linescan_corr(ll,:,i)-mean(V_data_linescan_corr(ll,:,i)))./mean(V_data_linescan_corr(ll,:,i));
V_DFF_movements = movmean(V_DFF_linescan(ll,:,i),25000,2) - mean(V_DFF_linescan(ll,:,i));
V_DFF_linescan(ll,:,i) = (V_DFF_linescan(ll,:,i)-V_DFF_movements(:,:));
C_DFF_linescan(ll,:,i) = (C_data_linescan(ll,1:trace_length,i)-mean(C_data_linescan(ll,1:trace_length,i)))./mean(C_data_linescan(ll,1:trace_length,i));
C_DFF_linescan(ll,:,i) = (C_DFF_linescan(ll,:,i))-V_DFF_movements(:,:);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear C_data_inter C_data_raw C_data_raw_in C_data_space_av C_data_space_temp_av...
      V_data_linescan_corr V_data_inter V_data_raw V_data_raw_in V_data_space_av V_data_space_temp_av
      %C_data_linescan V_data_linescan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hot pixel removal part2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1
for i = start_trace:traces    
V_DFF_linescan_mean = mean2(V_DFF_linescan(:,:,i));
V_DFF_linescan(V_DFF_linescan > 0.5) = V_DFF_linescan_mean;
V_DFF_linescan(V_DFF_linescan < -0.5) = V_DFF_linescan_mean;
end
V_DFF_linescan_unfiltered = V_DFF_linescan;
V_DFF_linescan = movmean((movmean(V_DFF_linescan,10,1)),50,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal analysis DCS DS detection and sub and supra threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go = 1;
if go == 1  
spikelet_count_trace = zeros(1,trace_length,traces);
DCS_spikelet_count_trace = zeros(1,trace_length,traces);
DCS_first_spikelet_count_trace = zeros(1,trace_length,traces);
DS_spikelet_count_trace = zeros(1,trace_length,traces);
DCS_cal_peak_trace = zeros(1,trace_length,traces);
DCS_cal_locs_trace = zeros(1,trace_length,traces);
DS_cal_peak_trace = zeros(1,trace_length,traces);
DS_cal_locs_trace = zeros(1,trace_length,traces);
V_DFF_space_av_filt_in = zeros(1,trace_length,traces);
C_DFF_space_av_filt_in = zeros(1,trace_length,traces);
C_DFF_space_av_filt_in_baseline_mins = zeros(1,10); 
C_DFF_space_av_filt_in_baseline = zeros(1,trace_length,traces);
V_DFF_space_av_sub = zeros(1,trace_length,traces);
V_DFF_hotspot_map = zeros(end_line-start_line+1,10000,traces);
DSeventcounter = 1;
DSeventcounter_stim = 1;
DCSeventcounter = 1;
DCSeventcounter_stim = 1;
min_hotspot_time = (30);
min_hotspot_space = (3);
max_hotspot_V = 0.07;
hotspot_threshold = 0.07;
for i = start_trace:traces   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % all spikelets/cal analysis 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    V_DFF_space_av_filt_in(:,:,i) = movmean(V_DFF_space_av(1,:,i),10);
    C_DFF_space_av_filt_in(:,:,i) = movmean(C_DFF_space_av(1,:,i),100);
    for secs = 1:10
    C_DFF_space_av_filt_in_baseline_mins(secs) = min(C_DFF_space_av_filt_in(1,(secs-1)*10000+1:secs*10000,i));
    end
    pline = polyfit(1:1:10,C_DFF_space_av_filt_in_baseline_mins,1);
    pline(1) = pline(1)/10000;
    C_DFF_space_av_filt_in_baseline(1,:,i) = ((1:1:trace_length)'*pline(1) + pline(2));
    C_DFF_space_av_filt_in(:,:,i) = C_DFF_space_av_filt_in(:,:,i) - C_DFF_space_av_filt_in_baseline(1,:,i);
    V_DFF_space_av_filt_in(:,:,i) = V_DFF_space_av_filt_in(:,:,i) - movmean(V_DFF_space_av_filt_in(:,:,i),500,2);
    [vpks, vlocs] = findpeaks(V_DFF_space_av_filt_in(:,:,i),'MinPeakHeight',std(V_DFF_space_av_filt_in(:,:,i))*3,'MinPeakDistance',5,'MinPeakProminence',std(V_DFF_space_av_filt_in(:,:,i)));
    spikelet_count_trace(1,vlocs,i) = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subthreshold analysis (hotspots)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_DFF_space_av_sub(:,:,i) = V_DFF_space_av(1,:,i);
    for tt = 101:trace_length-50
       if spikelet_count_trace(1,tt,i) == 1
       V_DFF_space_av_sub(:,tt-50:tt+50,i) = mean(V_DFF_space_av_sub(:,tt-100:tt-51,i));
       end
    end
    V_no_DCS = V_DFF_linescan(start_line:end_line,50200-4999:50200+5000,i);
    V_thres = V_no_DCS > hotspot_threshold;
    stats1 = regionprops('table',V_thres,V_no_DCS,'Centroid','MajorAxisLength','MinoraxisLength','Area','MaxIntensity','MeanIntensity','PixelList');
    stats1_out = stats1(stats1.MajorAxisLength > min_hotspot_time,:);
    stats1_out = stats1_out(stats1_out.MinorAxisLength > min_hotspot_space, :);
    stats1_out = stats1_out(stats1_out.MaxIntensity > max_hotspot_V, :);
    for hotspot = 1:size(stats1_out,1)
    V_DFF_hotspot_map(stats1_out.PixelList{hotspot}(:,2),stats1_out.PixelList{hotspot}(:,1),i) = 1;
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DSC DS spikelet counting only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for tt = 1+windowsize:trace_length-2*windowsize
        if spikelet_count_trace(1,tt,i) == 1 && sum(spikelet_count_trace(1,tt-100:tt+100,i),2) == 1 && mean(C_DFF_space_av_filt_in(1,tt:tt+200,i),2) > ((mean(C_DFF_space_av_filt_in(1,tt-200:tt,i),2)+0.01)) && mean(V_DFF_space_av_filt_in(:,tt:tt+100,i),2) < (std(V_DFF_space_av_filt_in(:,:,i),0,2))
        DS_spikelet_count_trace(1,tt,i) = 1;
        end
        if spikelet_count_trace(1,tt,i) == 1 && sum(spikelet_count_trace(1,tt-100:tt+100,i),2) > 1 && mean(C_DFF_space_av_filt_in(1,tt:tt+200,i),2) > ((mean(C_DFF_space_av_filt_in(1,tt-200:tt,i),2)+0.01))
        DCS_spikelet_count_trace(1,tt,i) = 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DSC DS splitting and sorting 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for tt = 1+windowsize:trace_length-2*windowsize  
    if spikelet_count_trace(1,tt,i) == 1 && sum(spikelet_count_trace(1,tt-100:tt+100,i),2) == 1 && mean(C_DFF_space_av_filt_in(1,tt:tt+200,i),2) > (mean(C_DFF_space_av_filt_in(1,tt-200:tt,i),2)+0.01) && mean(V_DFF_space_av_filt_in(:,tt:tt+100,i),2) < (std(V_DFF_space_av_filt_in(:,:,i),0,2))
    DS_spikelet_count_trace(1,tt,i) = 1;
    C_DFF_space_av_filt_in_section = (C_DFF_space_av_filt_in(:,tt:tt+windowsize,i) - mean(C_DFF_space_av_filt_in(:,tt-200:tt,i),2));
    [cpks, clocs] = max(C_DFF_space_av_filt_in_section);
    DS_cal_peak_trace(1,tt,i) = cpks(1,1);
    DS_cal_locs_trace(1,tt,i) = clocs(1,1);
    end   
    if spikelet_count_trace(1,tt,i) == 1 && sum(spikelet_count_trace(1,tt-100:tt+100,i),2) > 1 && mean(C_DFF_space_av_filt_in(1,tt:tt+200,i),2) > ((mean(C_DFF_space_av_filt_in(1,tt-200:tt,i),2)+0.01))
    DCS_spikelet_count_trace(1,tt,i) = 1;
    if spikelet_count_trace(1,tt,i) == 1 && sum(spikelet_count_trace(1,tt-100:tt,i),2) == 1
    DCS_first_spikelet_count_trace(1,tt,i) = 1;
    C_DFF_space_av_filt_in_section = (C_DFF_space_av_filt_in(:,tt:tt+windowsize,i) - mean(C_DFF_space_av_filt_in(:,tt-200:tt,i),2));
    [cpks, clocs] = max(C_DFF_space_av_filt_in_section);
    DCS_cal_peak_trace(1,tt,i) = cpks(1,1);
    DCS_cal_locs_trace(1,tt,i) = clocs(1,1);
    end
    end
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save traces for each cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
s1 = num2str(cell_number);
s2 = '.mat';
s3 = 'V_data_linescan_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'V_data_linescan','-v7.3');
s3 = 'C_data_linescan_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'C_data_linescan','-v7.3');
s3 = 'V_DFF_linescan_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'V_DFF_linescan_unfiltered','-v7.3');
s3 = 'C_DFF_linescan_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'C_DFF_linescan','-v7.3');
s3 = 'V_DFF_hotspot_map_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'V_DFF_hotspot_map');
s3 = 'V_DFF_space_av_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'V_DFF_space_av');
s3 = 'C_DFF_space_av_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'C_DFF_space_av');
s3 = 'spikelet_count_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'spikelet_count_trace');
s3 = 'DCS_spikelet_count_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'DCS_spikelet_count_trace');
s3 = 'DS_spikelet_count_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'DS_spikelet_count_trace');
s3 = 'DCS_cal_peak_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'DCS_cal_peak_trace');
s3 = 'DCS_cal_locs_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'DCS_cal_locs_trace');
s3 = 'DS_cal_peak_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'DS_cal_peak_trace');
s3 = 'DS_cal_locs_trace_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'DS_cal_locs_trace');
s3 = 'V_DFF_space_av_sub_cell_';
s_all = strcat(s3,s1,s2);
save(s_all,'V_DFF_space_av_sub');
if cell_number >= 8  
s3 = 'Eye_puff_int_cell_';
s_all = strcat(s3,s1,s2);  
save(s_all,'Eye_puff_int');
s3 = 'Eye_no_puff_int_cell_';
s_all = strcat(s3,s1,s2);  
save(s_all,'Eye_no_puff_int');
s3 = 'Eye_puff_drug_int_cell_';
s_all = strcat(s3,s1,s2);  
save(s_all,'Eye_puff_int_drug');
s3 = 'Eye_no_puff_drug_int_cell_';
s_all = strcat(s3,s1,s2);  
save(s_all,'Eye_no_puff_int_drug');
end
if cell_number < 8  
s3 = 'E_data_cell_';
s_all = strcat(s3,s1,s2);  
save(s_all,'E_data');
end
end
end
beep()
toc