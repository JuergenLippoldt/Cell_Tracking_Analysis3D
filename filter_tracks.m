
%% filter out the first 50 frames, because tissue is too dim there 

filtered_Tracks=Tracks;

for tr=1:length(filtered_Tracks)
     filtered_Tracks{tr}(filtered_Tracks{tr}(:,1)<50,:)=[];
    if size(filtered_Tracks{tr},2)<5 | size(filtered_Tracks{tr},1)<5
        filtered_Tracks{tr}=[];
    end
end


%% filter out whole tracks that are too bad

%set the parameters to what you think is appropriate
%Track_cell_vol_var_norm > 300
%Track_cell_vol_mad > 800
%Track_median_diff_vol_norm > 1.6


for tr=1:length(filtered_Tracks)
    if ~isempty(filtered_Tracks{tr})
        if Track_median_diff_vol_norm(tr) > 1.6 | Track_nuc_vol_var_norm(tr) > 300 | Track_nuc_vol_mad(tr) > 800 | min(median(filtered_Tracks{tr}(:,2:3)))==0 %| median(filtered_Tracks{tr}(:,2)) > 478
            filtered_Tracks{tr}=[];
        end
    end
end

sum(arrayfun(@(tr) ~isempty(filtered_Tracks{tr}),1:length(filtered_Tracks)))


time_filtered_Tracks=filtered_Tracks;

%% filter out individual times where good tracks are bad 

for tr=1:length(filtered_Tracks) % get rid of tracks with a certain difference to mean volume raito
    if ~isempty(time_filtered_Tracks{tr})
    time_filtered_Tracks{tr}(max([(time_filtered_Tracks{tr}(:,8)/median(time_filtered_Tracks{tr}(:,8))),(median(time_filtered_Tracks{tr}(:,8))./(time_filtered_Tracks{tr}(:,8)))]') >1.5,:)=[];
    end
end

for tr=1:length(filtered_Tracks) %get rid of tracks lower than a certian convex hull ratio
    if ~isempty(time_filtered_Tracks{tr})
    time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(:,19) < 0.8,:)=[];
    end
end

for tr=1:length(filtered_Tracks) %get rid of tracks lower than a certian convex hull ratio and with a certain difference to mean volume raito
    if ~isempty(time_filtered_Tracks{tr})
    time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(:,19) < 0.85 &  max([(time_filtered_Tracks{tr}(:,8)/median(time_filtered_Tracks{tr}(:,8))),(median(time_filtered_Tracks{tr}(:,8))./(time_filtered_Tracks{tr}(:,8)))]')' >1.4 ,:)=[] ; 
    end
end

for tr=1:length(filtered_Tracks) %get rid of tracks above and lower than certain volumes
    if ~isempty(time_filtered_Tracks{tr})
    time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(:,8) < 300 | time_filtered_Tracks{tr}(:,8) > 2000 ,:)=[];
    end
end

%% filter out individual times where good tracks are bad (volume only)
%filtering tracks based on volume changes, if the volume change too
%drastically the function will observe the frame before and after the
%change and delete the one that is too different from the mean.
%mild success


for tr=1:length(filtered_Tracks)
    if ~isempty(time_filtered_Tracks{tr})
       if size(time_filtered_Tracks{tr},2)>48
       vol_change = 200;
       tracks_with_big_vol_jumps = find(time_filtered_Tracks{tr}(:,49)> vol_change | time_filtered_Tracks{tr}(:,49)< -vol_change ) ; 
        for t = tracks_with_big_vol_jumps(end:-1:1)'
          if time_filtered_Tracks{tr}(t-1,8) - Track_median_cell_vol(tr)  > vol_change || time_filtered_Tracks{tr}(t-1,8) - Track_median_cell_vol(tr) < -vol_change 
             time_filtered_Tracks{tr}(t-1,:) = [];
          elseif time_filtered_Tracks{tr}(t,8) - Track_median_cell_vol(tr)  > vol_change || time_filtered_Tracks{tr}(t,8) - Track_median_cell_vol(tr) < -vol_change 
             time_filtered_Tracks{tr}(t,:) = [];
          end
        end
       end
    end
end


