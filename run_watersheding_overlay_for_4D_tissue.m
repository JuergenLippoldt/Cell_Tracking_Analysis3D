
% >> tracks should be loaded before and delete im_org if aquired previously
% since it might cause ram restrictions issue
% >> play around with diskradius in parameter section
% >> switch between fgm only for bw(black and white)(when you only want the
%nuclei for further analysis and between overlay when you need it for
%tracking 
% >> play around with fgm in '=0.X*graythresh' in 'finding nuclei' section for better fgm and with
% diskradius here at the start 
% >> change zslices in im_org_re depending on what you want and on the
% tracks in z 

%% parameter
im_org=[];
want_overlay_for_track=true; want_black_and_white_nuclei_segments=false; %choose between just a highlighting overlay or black and white nucsegments (flat image)
want_overlay_of_watersheding=false; Provide_Outlines = false; %cell outlines are only provided by doing the watershedding that is originally used for the visualization
NumberOfFrames=10;

x_umperpixel=0.867; %when tracks are in um then 0.867 
z_umperpixel=4; %depending on number of microns between each slice
ratio=[x_umperpixel,x_umperpixel,z_umperpixel];
disk_radius=round(12*.58/ratio(1));
upper_bound_nuc_area=round(20000*(.58/ratio(1))^2 );

tracks_in_um=true; %adjust this acordingly 

% Only needed if you need watershedding for want_overlay_of_watersheding or Provide_Outlines
original_non_sliced_filename='f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff' %insert orginal ometiff used for tracking here
original_non_sliced_filename_directory='C:\Users\sm69w\Documents\MATLAB\14th tumor s1\f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff'

% needed always
resliced_filename = 'resliced f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff' %insert resliced version here
resliced_filename_directory = 'C:\Users\sm69w\Documents\MATLAB\14th tumor s1\resliced f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff'

%% load original image
if want_overlay_of_watersheding  %isnt needed for tracking
   
filename = original_non_sliced_filename; 
fname = original_non_sliced_filename_directory;
info = imfinfo(fname);
num_images = numel(info);
for n = 1:num_images-1
%for n=1+50*21:51*21
    t=1+floor((n-1)/21);
    z=n-(t-1)*21;
    l=double(imread(fname, n, 'Info', info));
    im_org_watershed(:,:,z,t) = l/max(l(:));
  %  imorg(:,:,z) = l/max(l(:));
end
end
%some parameters


%%Key step for all processes. For loop over time - saving every timestep because of memory constrains  
for t=1:NumberOfFrames
        im_org_re=[];
        filename=resliced_filename; %insert resliced version here
        fname = resliced_filename_directory;
        info = imfinfo(fname);
        num_images = numel(info);
        for k_0 = 1:84
            k=k_0+(t-1)*84;    
            im_org_re(:,:,k_0) = imresize(imread(fname, k, 'Info', info), ratio(1));
        end

        % noise adjusting nuclei images
        nuc_org=im_org_re;
        for z=1:size(nuc_org,3)
            nuc_t_sl(:,:,z)=nuc_org(:,:,z)-imopen(nuc_org(:,:,z),strel('disk',2*disk_radius));
            nuc_t_sl(:,:,z)=adpmedian(nuc_t_sl(:,:,z), 3); 
        end
        nuc_t_sl = reshape(imadjust(nuc_t_sl(:)/max(nuc_t_sl(:))),size(nuc_t_sl));
        nuc = nuc_t_sl; %adjusted time slice of nuclei

        nuc=double(nuc)/double(max(nuc(:)));

%% Image processing/overlaying
    if want_overlay_for_track
        
        
        % intensity mask
        int_mask=[];
        for z=1:1:size(nuc,3)
            int_mask(:,:,z)=imclose(nuc(:,:,z),strel('disk',3*disk_radius));
        end    
        int_mask=int_mask/max(int_mask(:));
        %int_mask = reshape(double(int_mask(:))/max(double(int_mask(:))),size(int_mask));
        int_mask = imdilate(imdilate(int_mask >.4*graythresh(int_mask(int_mask>0)),strel('sphere',1)),strel('disk',3));    

        ws=round(disk_radius/2);
        sig_norm=imgaussfilt3(nuc,ws);
        outside2=imerode(sig_norm<.2*graythresh(sig_norm),strel('sphere',1));
        Gaussmask=imerode(imdilate(~outside2,strel('disk',3)),strel('disk',1));
        for z=1:size(Gaussmask,3)
            Gaussmask(:,:,z)=imfill_up2_size(Gaussmask(:,:,z),disk_radius^2);
        end
        Mask=logical(Gaussmask.*int_mask); %% tissue_mask
        
        
        % finding nuclei 
        fgm= nuc>0.5*graythresh(nuc(Mask)) & Mask; %fgm => fore-ground mask  one can play with number before greythresh
        for ws=1:10
            for z=1:size(fgm,3) %doing adaptive thresholding over multiple neighbourhoods 
                fgm(:,:,z) = fgm(:,:,z).*adaptivethreshold(nuc(:,:,z),3*ws^2*disk_radius);
            end
        end
        % imprioving mask - firstly fill hole
        for z=1:size(fgm,3)
        fgm(:,:,z)=imfill_up2_size(fgm(:,:,z),disk_radius^2);
        %fgm(:,:,z) = imfill(fgm(:,:,z),'holes');
        end
        % second step - morphological cleaning 
        fgm = imerode(fgm, strel('sphere',1));    
        fgm = imdilate(fgm, strel('sphere',1));
        for z=1:size(fgm,3)
            fgm(:,:,z)=imfill_up2_size(fgm(:,:,z),disk_radius^2);
        %fgm(:,:,z) = imfill(fgm(:,:,z),'holes');
        end    
        fgm = bwareaopen(fgm,10,6); % deleting small points    
        fgm = imerode(fgm, strel('sphere',1));
        fgm = imdilate(fgm, strel('sphere',1));

        %bw_nucs(:,:,:,t)=fgm;    %% building 4D image - might be to large -
        %thats why its commented out 

        %building overlay for tracking 
        overlay=0.4*nuc+0.5*fgm;
        %overlay=imresize(overlay,1/ratio(1)); %if you have issues with
        %size
        for z=1:size(overlay,3) %saves time-slice of overlay
           if want_black_and_white_nuclei_segments
                imwrite(double(fgm(:,:,z)),['bw_nucsegment_',filename],'WriteMode','append')
           else
                imwrite(overlay(:,:,z),['overlay_nucsegment_',filename],'WriteMode','append')
           end
        end
    end
    
%% Watershedding for visualization and cell outlines
    if want_overlay_of_watersheding
        %% acuring binary images from FrameTrackCoordinates
        %init volumes
        
        %nuc= imresize(nuc, ratio(1));
        
        
        % intensity mask
        int_mask=[];
        for z=1:1:size(nuc,3)
            int_mask(:,:,z)=imclose(nuc(:,:,z),strel('disk',3*disk_radius));
        end    
        int_mask=int_mask/max(int_mask(:));
        %int_mask = reshape(double(int_mask(:))/max(double(int_mask(:))),size(int_mask));
        int_mask = imdilate(imdilate(int_mask >.4*graythresh(int_mask(int_mask>0)),strel('sphere',1)),strel('disk',3));    

        ws=round(disk_radius/2);
        sig_norm=imgaussfilt3(nuc,ws);
        outside2=imerode(sig_norm<.2*graythresh(sig_norm),strel('sphere',1));
        Gaussmask=imerode(imdilate(~outside2,strel('disk',3)),strel('disk',1));
        for z=1:size(Gaussmask,3)
            Gaussmask(:,:,z)=imfill_up2_size(Gaussmask(:,:,z),disk_radius^2);
        end
        Mask=logical(Gaussmask.*int_mask); %% tissue_mask
        
        

            track_points = zeros(size(nuc));

        if tracks_in_um
            FrameTrackCoordinates_discretized = cellfun(@round, FrameTrackCoordinates, 'UniformOutput', false);
        else
            FrameTrackCoordinates_discretized=cellfun(@(x) x.*[x_umperpixel,x_umperpixel,1],FrameTrackCoordinates,'un',0);  %Z-skale for pixel is always corrected
            FrameTrackCoordinates_discretized = cellfun(@round, FrameTrackCoordinates_discretized, 'UniformOutput', false);
      
        end
        %index the first track
        currentTrack = FrameTrackCoordinates_discretized{t}+1;%add one because fiji starts counting with 0

        %linearIndex = sub2ind(size(FrameTrackCoordinates_discretized), currentTrack(:, 1), currentTrack(:, 2), currentTrack(:, 3));
        for i = 1:length(currentTrack)
            track_points( currentTrack(i, 2), currentTrack(i, 1), currentTrack(i, 3)) = 1;
        end

        track_points=imdilate(track_points,strel('sphere',2)); %so that outlines smooth has it that the duplicated tracks over one point are merged to avoid issues
        track_points=imerode(track_points,strel('sphere',1)); %erode them so they arent that big

        %actual watershedding
        L = watershed(~track_points);
        L(Mask==0)=0;
        %get outlines as shape estimates
        I3 = zeros(size(track_points));
        I3(imdilate(L == 0, ones(1, 1))) = 1;

        outlines=zeros(size(I3));
        outlines_smooth=zeros(size(I3));
        %some parameters
        smoothy=1;
        lower_bound_cell_area=40; 
        %% smooth and adjust
        for k=1:size(I3,3)
        outlines(:,:,k)=bwareaopen(I3(:,:,k), 80);% 
        outlines_smooth(:,:,k)=outlines(:,:,k);
        for i=1:smoothy
        outlines_smooth(:,:,k)=bwmorph(outlines_smooth(:,:,k),'dilate');
        end
        for i=1:smoothy
        outlines_smooth(:,:,k)=bwmorph(outlines_smooth(:,:,k),'erode');
        end
        holes=ones(size(outlines_smooth(:,:,k)))-outlines_smooth(:,:,k);
        bigholes = bwareaopen(holes, lower_bound_cell_area);
        smallholes = holes & ~bigholes;
        outlines_smooth(:,:,k)=outlines_smooth(:,:,k) | smallholes;
        
        outlines_smooth(1,:,k)=1;
        outlines_smooth(end,:,k)=1;
        outlines_smooth(:,1,k)=1;
        outlines_smooth(:,end,k)=1;
        end


        overlay_watershed=.5*nuc+.5*outlines_smooth;
        for z=1:size(overlay_watershed,3) %saves time-slice of overlay
            if Provide_Outlines == false
            imwrite(overlay_watershed(:,:,z),['overlay_watershed_',filename],'WriteMode','append')
            elseif Provide_Outlines == true
            imwrite(outlines_smooth(:,:,z),['outlines_smooth2_',filename],'WriteMode','append')
            end
        end


    end
   
end

if false
   test=(cellfun(@(tr) max(Tracks{tr}(:,2),1:length(Tracks),'UniformOutput',false)));
    max(cell2mat(cellfun(@max ,Tracks ,'UniformOutput' ,false)));
end
