% get some idea of overall statistics and which tracks are interesting 
%working with NaN is better here than zeros so that we don't inflate the velocity with zeros 
TracksAreFiltered = true;
AddVelocities = true; %we do have to filter to add the velocities %sometimes I need track stats without additional velocities, might be a good idea to incorporate the add veloctiies code in here too
Tracks_for_stats = time_filtered_Tracks;
if TracksAreFiltered
        TrackCell_length=NaN(1,length(Tracks_for_stats));
        TrackCell_median__cell_vol=NaN(1,length(Tracks_for_stats));
        TrackCell_median_dif_to_median_cell_vol=NaN(1,length(Tracks_for_stats));
        TrackCell_vol_var=NaN(1,length(Tracks_for_stats));
        TrackCell_vol_var_norm=NaN(1,length(Tracks_for_stats));
        TrackCell_vol_mad=NaN(1,length(Tracks_for_stats)); %median absloute deviation
        TrackCell_most_diff_vol_norm=NaN(1,length(Tracks_for_stats));
        TrackCell_median_diff_vol_norm=NaN(1,length(Tracks_for_stats));
        TrackCell_most_diff_time=NaN(1,length(Tracks_for_stats));
        TrackCell_median_shape=NaN(1,length(Tracks_for_stats));
        if AddVelocities 
            

            TrackCell_Pearson_shape_vel=NaN(1,length(Tracks_for_stats)); %Pearson correlation between Vel and Shape

            
        end
    %find(Track_cell_vol_var_norm>3000);

    for tr=1:length(Tracks_for_stats)
        if size(Tracks_for_stats{tr},2)>7 && size(Tracks_for_stats{tr},1)>10
            TrackCell_length(tr)=length(Tracks_for_stats{tr});
            TrackCell_median__cell_vol(tr)=median(Tracks_for_stats{tr}(:,23));
            TrackCell_median_dif_to_median_cell_vol(tr)=median(abs(Tracks_for_stats{tr}(:,23)-median(Tracks_for_stats{tr}(:,23))))/median(Tracks_for_stats{tr}(:,23));
            TrackCell_vol_var(tr)=var(Tracks_for_stats{tr}(:,23));
            TrackCell_vol_var_norm(tr)=var(Tracks_for_stats{tr}(:,23))/mean(Tracks_for_stats{tr}(:,23));
            TrackCell_vol_mad(tr)=mad(Tracks_for_stats{tr}(:,23),1);    
            TrackCell_most_diff_vol_norm(tr)=max(max(Tracks_for_stats{tr}(:,23)/median(Tracks_for_stats{tr}(:,23))),max(median(Tracks_for_stats{tr}(:,23))./(Tracks_for_stats{tr}(:,23))));
            TrackCell_median_diff_vol_norm(tr)=median(max([(Tracks_for_stats{tr}(:,23)/median(Tracks_for_stats{tr}(:,23))),(median(Tracks_for_stats{tr}(:,23))./(Tracks_for_stats{tr}(:,23)))]'));
            [~,TrackCell_most_diff_time(tr)]=max(abs(Tracks_for_stats{tr}(:,23)-median(Tracks_for_stats{tr}(:,23))));    
            TrackCell_median_shape(tr)=median(Tracks_for_stats{tr}(:,22));  
            if AddVelocities
                intermed=corrcoef(Tracks_for_stats{tr}(1:end-1,22),sum((Tracks_for_stats{tr}(2:end,2:4)-Tracks_for_stats{tr}(1:end-1,2:4)).^2,2).^.5);
                TrackCell_Pearson_shape_vel(tr)=intermed(1,2); %A=shape,B=velocity                
                
            end
        end
    end 
else
        TrackCell_length=NaN(1,length(Tracks_for_stats));
        TrackCell_median__cell_vol=NaN(1,length(Tracks_for_stats));
        TrackCell_median_dif_to_median_cell_vol=NaN(1,length(Tracks_for_stats));
        TrackCell_vol_var=NaN(1,length(Tracks_for_stats));
        TrackCell_vol_var_norm=NaN(1,length(Tracks_for_stats));
        TrackCell_vol_mad=NaN(1,length(Tracks_for_stats)); %median absloute deviation
        TrackCell_most_diff_vol_norm=NaN(1,length(Tracks_for_stats));
        TrackCell_median_diff_vol_norm=NaN(1,length(Tracks_for_stats));
        TrackCell_most_diff_time=NaN(1,length(Tracks_for_stats));
        TrackCell_median_shape=NaN(1,length(Tracks_for_stats));
        TrackCell_Pearson_shape_vel=NaN(1,length(Tracks_for_stats)); %Pearson correlation between Vel and Shape
        
    for tr=1:length(Tracks_for_stats)
        if size(Tracks_for_stats{tr},2)>7 && size(Tracks_for_stats{tr},1)>10
            TrackCell_length(tr)=length(Tracks_for_stats{tr});
            TrackCell_median__cell_vol(tr)=median(Tracks_for_stats{tr}(:,23));
            TrackCell_median_dif_to_median_cell_vol(tr)=median(abs(Tracks_for_stats{tr}(:,23)-median(Tracks_for_stats{tr}(:,23))))/median(Tracks_for_stats{tr}(:,23));
            TrackCell_vol_var(tr)=var(Tracks_for_stats{tr}(:,23));
            TrackCell_vol_var_norm(tr)=var(Tracks_for_stats{tr}(:,23))/mean(Tracks_for_stats{tr}(:,23));
            TrackCell_vol_mad(tr)=mad(Tracks_for_stats{tr}(:,23),1);    
            TrackCell_most_diff_vol_norm(tr)=max(max(Tracks_for_stats{tr}(:,23)/median(Tracks_for_stats{tr}(:,23))),max(median(Tracks_for_stats{tr}(:,23))./(Tracks_for_stats{tr}(:,23))));
            TrackCell_median_diff_vol_norm(tr)=median(max([(Tracks_for_stats{tr}(:,23)/median(Tracks_for_stats{tr}(:,23))),(median(Tracks_for_stats{tr}(:,23))./(Tracks_for_stats{tr}(:,23)))]'));
            [~,TrackCell_most_diff_time(tr)]=max(abs(Tracks_for_stats{tr}(:,23)-median(Tracks_for_stats{tr}(:,23))));    
            TrackCell_median_shape(tr)=median(Tracks_for_stats{tr}(:,22));
            intermed=corrcoef(Tracks_for_stats{tr}(1:end-1,22),sum((Tracks_for_stats{tr}(2:end,2:4)-Tracks_for_stats{tr}(1:end-1,2:4)).^2,2).^.5);
            TrackCell_Pearson_shape_vel(tr)=intermed(1,2); %A=shape,B=velocity
        end
    end
end 
