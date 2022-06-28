shape_file='C:\Users\sm69w\Documents\MATLAB\14th tumor s1\nuc_shapes2_resliced f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff';
save_folder='C:\Users\sm69w\Documents\MATLAB\14th tumor s1\Individual filtered nucs movie\'%'C:\Users\sm69w\Documents\MATLAB\2nd scene worked on in 14th tumor\Individual nucs movie\';
save_name='14th_MCA_s1'
z_dim=84;
%this function will produce matlab files that can be accessed with imshow4
%to insepct indiviual tracked nuclei

%% find which tracks to plot depending on track_stats 
if false 
track_list_to_plot=[find(Track_max_shape_change>.5 ,5)] ;
end

%% 

track_list_to_plot = [3433] %put the track numbers of the tracks you want to observe here
for tr_num=track_list_to_plot(1:end)  %fixed a bug and started with the track that gave the error JL
    track=Tracks{tr_num};%Tracks{tr_num};
   % disp(tr)
    plot_tracked_nuc(track,shape_file,save_folder,save_name,z_dim,tr_num)
end
