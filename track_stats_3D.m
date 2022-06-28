% get some idea of overall statistics and which tracks are interesting 
%working with NaN is better here than zeros so that we don't inflate the velocity with zeros 
TracksAreFiltered = false;
AddVelocities = false; %we do have to filter to add the velocities %sometimes I need track stats without additional velocities, might be a good idea to incorporate the add veloctiies code in here too
if TracksAreFiltered 
    OG_Tracks= Tracks; 
    Tracks = time_filtered_Tracks;
        Track_length=NaN(1,length(Tracks));
        Track_median_nuc_vol=NaN(1,length(Tracks));
        Track_mean_nuc_vol=NaN(1,length(Tracks));
        Track_median_dif_to_median_nuc_vol=NaN(1,length(Tracks));
        Track_nuc_vol_var=NaN(1,length(Tracks));
        Track_nuc_vol_var_norm=NaN(1,length(Tracks));
        Track_nuc_vol_mad=NaN(1,length(Tracks)); %median absloute deviation
        Track_most_diff_vol_norm=NaN(1,length(Tracks));
        Track_median_diff_vol_norm=NaN(1,length(Tracks));
        Track_most_diff_time=NaN(1,length(Tracks));
        Track_median_nuc_shape=NaN(1,length(Tracks));
        Track_mean_nuc_shape=NaN(1,length(Tracks));
        Track_mean_Track_velocity=NaN(1,length(Tracks));
        Track_median_Track_velocity=NaN(1,length(Tracks));
        Track_max_velocity=NaN(1,length(Tracks));
        Track_mean_rel_Track_velocity=NaN(1,length(Tracks));

        if AddVelocities 
            
            Track_mean_centroid_rel_velocity=NaN(1,length(Tracks));
            Track_median_centroid_rel_velocity=NaN(1,length(Tracks));
            Track_max_rel_centroid_velocity=NaN(1,length(Tracks));
            Track_Pearson_shape_vel=NaN(1,length(Tracks)); %Pearson correlation between Vel and Shape
            Track_Pearson_shape_vel_1h=NaN(1,length(Tracks));
            Track_Pearson_shape_vel_2h=NaN(1,length(Tracks));
            Track_Pearson_shape_rel_vel=NaN(1,length(Tracks)); 
            Track_Pearson_ar_rel_vel=NaN(1,length(Tracks)); 
            Track_Pearson_shape_rel_vel_1h=NaN(1,length(Tracks)); 
            Track_Pearson_ar_rel_vel_1h=NaN(1,length(Tracks)); 
            Track_Pearson_shape_rel_vel_2h=NaN(1,length(Tracks)); 
            Track_Pearson_ar_rel_vel_2h=NaN(1,length(Tracks)); 
            Track_Pearson_nucvol_rel_vel=NaN(1,length(Tracks)); 
            Track_Pearson_avg_shape_avg_rel_vel2h=NaN(1,length(Tracks));
            Track_max_shape_change=NaN(1,length(Tracks));
            all_filtered_Tracks=[];
          	Track_mean_centroid_rel_velocity1h=NaN(1,length(Tracks));
            Track_median_centroid_rel_velocity1h=NaN(1,length(Tracks));
            Track_mean_centroid_rel_velocity2h=NaN(1,length(Tracks));
            Track_median_centroid_rel_velocity2h=NaN(1,length(Tracks));
          	Track_var_centroid_rel_velocity1h=NaN(1,length(Tracks));
            Track_var_centroid_rel_velocity2h=NaN(1,length(Tracks));
          	Track_std_centroid_rel_velocity1h=NaN(1,length(Tracks));
            Track_std_centroid_rel_velocity2h=NaN(1,length(Tracks));
        end
    %find(Track_cell_vol_var_norm>3000);

    for tr=1:length(Tracks)
        if size(Tracks{tr},2)>7 && size(Tracks{tr},1)>10
            Track_length(tr)=size(Tracks{tr},1);
            Track_median_nuc_vol(tr)=median(Tracks{tr}(:,8));
            Track_mean_nuc_vol(tr)=nanmean(Tracks{tr}(:,8));
            Track_median_dif_to_median_nuc_vol(tr)=median(abs(Tracks{tr}(:,8)-median(Tracks{tr}(:,8))))/median(Tracks{tr}(:,8));
            Track_nuc_vol_var(tr)=var(Tracks{tr}(:,8));
            Track_nuc_vol_var_norm(tr)=nanvar(Tracks{tr}(:,8))/nanmean(Tracks{tr}(:,8));
            Track_nuc_vol_mad(tr)=mad(Tracks{tr}(:,8),1);    
            Track_most_diff_vol_norm(tr)=max(max(Tracks{tr}(:,8)/median(Tracks{tr}(:,8))),max(median(Tracks{tr}(:,8))./(Tracks{tr}(:,8))));
            Track_median_diff_vol_norm(tr)=median(max([(Tracks{tr}(:,8)/median(Tracks{tr}(:,8))),(median(Tracks{tr}(:,8))./(Tracks{tr}(:,8)))]'));
            [~,Track_most_diff_time(tr)]=max(abs(Tracks{tr}(:,8)-median(Tracks{tr}(:,8))));    
            Track_median_nuc_shape(tr)=median(Tracks{tr}(:,7));
            Track_mean_nuc_shape(tr)=nanmean(Tracks{tr}(:,7));
            Track_mean_Track_velocity(tr)=nanmean(sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);
            Track_median_Track_velocity(tr)=median(sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);
            Track_max_velocity(tr)=max(sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5); 
            Track_mean_rel_Track_velocity(tr)=nanmean(Tracks{tr}(:,50));
            if AddVelocities
                
                Track_mean_centroid_rel_velocity(tr)=nanmean(Tracks{tr}(:,36));
                Track_median_centroid_rel_velocity(tr)=nanmedian(Tracks{tr}(:,36));
                Track_max_rel_centroid_velocity(tr)=max(Tracks{tr}(:,36));
                
                intermed=corrcoef(Tracks{tr}(1:end-1,7),sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);
                Track_Pearson_shape_vel(tr)=intermed(1,2); %A=shape,B=velocity                

                Track_max_shape_change(tr)=max(abs(time_filtered_Tracks{tr}(:,48)));  
                
                intermed=corrcoef(Tracks{tr}(1:end-1,8),Tracks{tr}(1:end-1,36),'Rows','complete');
                Track_Pearson_nucvol_rel_vel(tr)=intermed(1,2); %A=shape,B=velocity
                
                intermed=corrcoef(Tracks{tr}(1:end-1,7),Tracks{tr}(1:end-1,36),'Rows','complete');
                Track_Pearson_shape_rel_vel(tr)=intermed(1,2); %A=shape,B=velocity

                intermed=corrcoef(Tracks{tr}(1:end-1,6),Tracks{tr}(1:end-1,36),'Rows','complete');
                Track_Pearson_ar_rel_vel(tr)=intermed(1,2); %A=shape,B=velocity


                if size(Tracks{tr},2)>7 && size(Tracks{tr},1)>18
                intermed=corrcoef(Tracks{tr}(1:end-1,7),Tracks{tr}(1:end-1,37),'Rows','complete');
                Track_Pearson_shape_vel_1h(tr)=intermed(1,2); %A=shape,B=velocity
                
                intermed=corrcoef(Tracks{tr}(1:end-1,7),Tracks{tr}(1:end-1,39),'Rows','complete');
                Track_Pearson_shape_rel_vel_1h(tr)=intermed(1,2); %A=shape,B=velocity

                intermed=corrcoef(Tracks{tr}(1:end-1,6),Tracks{tr}(1:end-1,39),'Rows','complete');     
                Track_Pearson_ar_rel_vel_1h(tr)=intermed(1,2); %A=shape,B=velocity

                end
                
                if size(Tracks{tr},2)>7 && size(Tracks{tr},1)>30

                intermed=corrcoef(Tracks{tr}(1:end-1,7),Tracks{tr}(1:end-1,38),'Rows','complete');
                Track_Pearson_shape_vel_2h(tr)=intermed(1,2); %A=shape,B=velocity


                intermed=corrcoef(Tracks{tr}(1:end-1,7),Tracks{tr}(1:end-1,40),'Rows','complete');
                Track_Pearson_shape_rel_vel_2h(tr)=intermed(1,2); %A=shape,B=velocity

                intermed=corrcoef(Tracks{tr}(1:end-1,6),Tracks{tr}(1:end-1,40),'Rows','complete');
                Track_Pearson_ar_rel_vel_2h(tr)=intermed(1,2); %A=shape,B=velocity


                intermed=corrcoef(Tracks{tr}(1:end-1,43),Tracks{tr}(1:end-1,46),'Rows','complete');
                Track_Pearson_avg_shape_avg_rel_vel2h(tr)=intermed(1,2); %A=shape,B=velocity
                
                end
                
             	Track_mean_centroid_rel_velocity1h(tr)=nanmean(Tracks{tr}(:,39));
                Track_median_centroid_rel_velocity1h(tr)=nanmedian(Tracks{tr}(:,39));
                Track_mean_centroid_rel_velocity2h(tr)=nanmean(Tracks{tr}(:,40));
                Track_median_centroid_rel_velocity2h(tr)=nanmedian(Tracks{tr}(:,40));
                Track_var_centroid_rel_velocity1h(tr)=nanvar(Tracks{tr}(:,39));
                Track_var_centroid_rel_velocity2h(tr)=nanvar(Tracks{tr}(:,40));
          
                Track_std_centroid_rel_velocity1h(tr)=nanstd(Tracks{tr}(:,39));
                Track_std_centroid_rel_velocity2h(tr)=nanstd(Tracks{tr}(:,40));  
                
                all_filtered_Tracks=[all_filtered_Tracks; Tracks{tr}];
            end
        end
    end
Tracks = OG_Tracks;    

motility_measure2= (Track_mean_centroid_rel_velocity2h(:)-nanmean(Track_mean_centroid_rel_velocity2h(:)))/nanstd(Track_mean_centroid_rel_velocity2h) + (Track_var_centroid_rel_velocity2h(:)-nanmean(Track_var_centroid_rel_velocity2h(:)))/nanstd(Track_var_centroid_rel_velocity2h(:)) ;

else
OG_Tracks = Tracks;
        Track_length=NaN(1,length(Tracks));
        Track_median_nuc_vol=NaN(1,length(Tracks));
        Track_median_dif_to_median_nuc_vol=NaN(1,length(Tracks));
        Track_nuc_vol_var=NaN(1,length(Tracks));
        Track_nuc_vol_var_norm=NaN(1,length(Tracks));
        Track_nuc_vol_mad=NaN(1,length(Tracks)); %median absloute deviation
        Track_most_diff_vol_norm=NaN(1,length(Tracks));
        Track_median_diff_vol_norm=NaN(1,length(Tracks));
        Track_most_diff_time=NaN(1,length(Tracks));
        Track_median_nuc_shape=NaN(1,length(Tracks));
        Track_mean_Track_velocity=NaN(1,length(Tracks));
        Track_median_Track_velocity=NaN(1,length(Tracks));
        Track_max_velocity=NaN(1,length(Tracks));
        Track_Pearson_shape_vel=NaN(1,length(Tracks)); %Pearson correlation between Vel and Shape
       
        
    for tr=1:length(Tracks)
        if size(Tracks{tr},2)>7 && size(Tracks{tr},1)>10
            Track_length(tr)=size(Tracks{tr},1);
            Track_median_nuc_vol(tr)=median(Tracks{tr}(:,8));
            Track_median_dif_to_median_nuc_vol(tr)=median(abs(Tracks{tr}(:,8)-median(Tracks{tr}(:,8))))/median(Tracks{tr}(:,8));
            Track_nuc_vol_var(tr)=var(Tracks{tr}(:,8));
            Track_nuc_vol_var_norm(tr)=nanvar(Tracks{tr}(:,8))/nanmean(Tracks{tr}(:,8));
            Track_nuc_vol_mad(tr)=mad(Tracks{tr}(:,8),1);    
            Track_most_diff_vol_norm(tr)=max(max(Tracks{tr}(:,8)/median(Tracks{tr}(:,8))),max(median(Tracks{tr}(:,8))./(Tracks{tr}(:,8))));
            Track_median_diff_vol_norm(tr)=median(max([(Tracks{tr}(:,8)/median(Tracks{tr}(:,8))),(median(Tracks{tr}(:,8))./(Tracks{tr}(:,8)))]'));
            [~,Track_most_diff_time(tr)]=max(abs(Tracks{tr}(:,8)-median(Tracks{tr}(:,8))));    
            Track_median_nuc_shape(tr)=median(Tracks{tr}(:,7));
            Track_mean_Track_velocity(tr)=nanmean(sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);
            Track_median_Track_velocity(tr)=median(sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);
            Track_max_velocity(tr)=max(sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);   
            intermed=corrcoef(Tracks{tr}(1:end-1,7),sum((Tracks{tr}(2:end,2:4)-Tracks{tr}(1:end-1,2:4)).^2,2).^.5);
            Track_Pearson_shape_vel(tr)=intermed(1,2); %A=shape,B=velocity
        end
    end
end 
Tracks = OG_Tracks;

if false 
  	long_tracks=[];
    for tr=1:length(time_filtered_Tracks)
        if length(time_filtered_Tracks{tr})==max(Track_length)
            long_tracks{end+1}=time_filtered_Tracks{tr};
        end
    end
end

