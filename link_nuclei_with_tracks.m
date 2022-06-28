        % Tracks 1 time in frames starting with 0 
        % Tracks 2 x postition 
        % Tracks 3 y postition 
        % Tracks 4 z postition 
        % Tracks 5 internal control variable 
        % Tracks 6 aspect ratio nuclei
        % Tracks 7 ellipsoid shape parameter
        % Tracks 8 volume in voxel
        % Tracks 9-11 xyz of centroid of nucleus

%% changes to keep in mind, delete if not resolved
%added . to L/max(L) in line 62+3

%% define properties 

%% code is writen for tracks-file in pixel unit !!!
%% look at ratios - something is strange 
ratio=[0.867, 0.867,1]; % micrometer per pixel 
z_dim=84;
load_cell_outlines=false %in case you have cell outlines data and want to anaylse cells and not just nuclei, this option was mostly kept false in previous analysis since the focus was on nuclei

trackname='C:\Users\sm69w\Documents\MATLAB\14th tumor s1\Full sequence\resliced 14th tumor MCA106 5h-post-OP 5min_s1_Tracks 9 40 15 mean 1 std 1 quality 1 min 1 median 1.1 - 20 2 mean 2.5 std 2.5 quality 2 min 2.5 median 2.6 [um].xml'; %track file obtained from fiji
d4_bw_image_name='C:\Users\sm69w\Documents\MATLAB\14th tumor s1\Full sequence\bw_nucsegment_resliced 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff'; %insert bw nucsegment file name here
resliced_nuc_name_f='resliced 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff'; %insert resliced file name here
outlines_smooth_name_f=''; %optional, leave empty if load_cell_outlines=false
direct='C:\Users\sm69w\Documents\MATLAB\14th tumor s1\Full sequence\'; %file directory where all the files ares
resliced_nuc_name=[direct,resliced_nuc_name_f];
outlines_smooth_name=[direct,outlines_smooth_name_f];

delete_shorter_tracks=true; shorter_than=11; %insert the length of tracks you want omitted from the data. Keep in mind its better to keep shorter tracks for the segmentation and to not include them later in the anaylsis




%% load tracks and bw_nuclei

[Tracks, ~] = importTrackMateTracks(trackname);
%     if trackname =='C:\Users\sm69w\Documents\MATLAB\2nd scene worked on in 14th tumor\overlay_14th tumor MCA 106 5h-post-OP 5min_Tracks.xml' % special case only 
%         
%         Tracks=cellfun(@(x) x./[1,1.07,1.07,1],Tracks,'un',0);  %Z-skale for pixel is always corrected
%     end 
info = imfinfo(d4_bw_image_name);
num_images = numel(info);
for n = 1:num_images-1
%for n=1+50*21:51*21
    t=1+floor((n-1)/z_dim);
    z=n-(t-1)*z_dim;
    im_org(:,:,z,t)=logical(imread(d4_bw_image_name, n, 'Info', info)); %no need to resize since bw is already resized to be in micron in run_watershedding
    
  %  imorg(:,:,z) = l/max(l(:));
end


%% load outlines_smooth
if load_cell_outlines==true
info = imfinfo(outlines_smooth_name);
num_images = numel(info);
for n = 1:num_images-1
%for n=1+50*21:51*21
    t=1+floor((n-1)/z_dim);
    z=n-(t-1)*z_dim;
    outlines_smooth(:,:,z,t)=logical(imread(outlines_smooth_name, n, 'Info', info)); %no need to resize since it is already resized to be in micron in run_watershedding
end
end

%% clean up tracks 

Tracks=fill_gaps3D(Tracks);

% take out tracks shorter than 10 frames
if delete_shorter_tracks
Tracks(cellfun('length',Tracks) <shorter_than)=[];
end

T = Tracks;
% This tells for each frame the tracks of which it is made up
[FrameTracks, FrameTrackCoordinates] = FindAllTracksInFrames(T);
clear T
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );


for f = 1:length(FrameTracks)
    %disp(f);
    
    % find the vertex coordinates and the cells of the tesselation
    [V,C]	= voronoin( FrameTrackCoordinates{f} );
    
    %C(cellfun(@isempty,C))=[];
    
    % find the neighbors of each vertex and each cell of the tesselation
    [VN,CN] = FindNeighbors( C );
    
    % collect data
    Vertices{f}			= V;
    VertexNeighbors{f}	= VN;
    Cells{f}			= C;
    CellNeighbors{f}	= CN;
    
    clearvars C V VN CN
    
    %	waitbar(f/length(FrameTracks));
end
%close(waitbar(0));

% now we have to associate this to the actual Tracks!
TrackNeighborTracks = cell( size(Tracks) );

TrackNeighborCells = cell( size(Tracks) );
for f=1:length(FrameTracks)
    for t=1:length(FrameTracks{f})
        TrackNeighborCells{FrameTracks{f}(t)}{end+1}	= CellNeighbors{f}{t};
        TrackNeighborTracks{FrameTracks{f}(t)}{end+1}	= FrameTracks{f}(CellNeighbors{f}{t});
    end
end


%% segment nuclei using tracks 
FrameTrackCoordinates_discretized = cellfun(@round, FrameTrackCoordinates, 'UniformOutput', false);
%index the first track


for t=1:3%size(im_org,4) % 1:
    
    track_points=zeros(size(im_org,1),size(im_org,2),size(im_org,3));

currentTrack = FrameTrackCoordinates_discretized{t}+1;%add one because fiji starts counting with 0

	%linearIndex = sub2ind(size(FrameTrackCoordinates_discretized), currentTrack(:, 1), currentTrack(:, 2), currentTrack(:, 3));
	for i = 1:length(currentTrack)
        track_points( currentTrack(i, 2), currentTrack(i, 1), currentTrack(i, 3)) = 1;
    end
    
    
    
    %% load nuc_frame
if false   %we tried out a segmentation based on the original images and failed so far
    info = imfinfo(resliced_nuc_name);
    num_images = numel(info);
    for n = 1:num_images-1
    %for n=1+50*z_dim:51*z_dim
    t_nuc=1+floor((n-1)/z_dim);
    if t_nuc==t
        z=n-(t_nuc-1)*z_dim;
        nuc_org(:,:,z)=double(imresize(imread(resliced_nuc_name, n, 'Info', info),ratio(1))); %make the ratio in micron since tracks are in micron
    end
  %  imorg(:,:,z) = l/max(l(:));
    end
    nuc=nuc_org/max(nuc_org(:));

    [Gx_nuc,Gy_nuc,Gz_nuc] = imgradientxyz(nuc);

    gradmag_nuc = Gx_nuc.^2 + Gy_nuc.^2 + Gz_nuc.^2;% + Gx_nuc.^2 + Gy_nuc.^2 + Gz_nuc.^2; %using this to help segment the nuclei more accurately

    IM=nuc-imgaussfilt3(nuc,10);

    lehisto=imhist(reshape(IM,[size(nuc,1)*size(nuc,3),size(nuc,2)]),100);
         level=otsuthresh(lehisto);
    outside = IM<level;
    outside = ~bwareaopen(~outside, 25);
    outside = bwareaopen(outside, 25);

    outside3 = ~imerode(imdilate(imdilate(~outside+track_points,strel('sphere',6)),strel('disk',9)),strel('disk',1));

    gradmag2 = imimposemin(gradmag_nuc,outside3 | track_points); % use this if you estimate the background
    %
    L = watershed(gradmag2);
    L = double(L).*double(squeeze(im_org(:,:,:,t)));
    I = zeros(size(gradmag2));
    I(L == 0) = 1;

    outlines= I > 0.5;
end

        pseudo_gradmag= (bwdist(~imdilate(squeeze(im_org(:,:,:,t)),strel('sphere',2)))).^0.5; %watersheding works as if you are filling a landscape, it looks like we get better results when we make the landscape less steeper
        pseudo_gradmag(pseudo_gradmag>0)=max(pseudo_gradmag(:))./pseudo_gradmag(pseudo_gradmag>0);
        pseudo_gradmag=pseudo_gradmag/max(pseudo_gradmag(:));
        


        towatershed = imimposemin(pseudo_gradmag, imdilate(track_points,strel('sphere',2)));
        L = watershed(towatershed);
        L = double(L).*double(squeeze(im_org(:,:,:,t)));
        L=L/max(L(:));
        boundaries=L==0;
        
        nuc_shapes=L>0;
        
      	for z=1:size(nuc_shapes,3) %saves time-slice of overlay
            imwrite(double(nuc_shapes(:,:,z)),['nuc_shapes2_',resliced_nuc_name_f],'WriteMode','append')

        end

        
    if false
    overlay(:,:,:,1)=imdilate(track_points,strel('sphere',2));
    overlay(:,:,:,2)=L;
    overlay(:,:,:,3)=track_points*0;
    imshow3D(overlay)
    end


%% cleanup 

%% get nuclei properties 
CC_nuc=bwconncomp(nuc_shapes,6); 
nuc_shape_props=regionprops(CC_nuc);
nuc_shape_props3=regionprops3(CC_nuc,'ConvexVolume');

if load_cell_outlines==true
CC_cell=bwconncomp(~squeeze(outlines_smooth(:,:,:,t)),6); 
cell_shape_props=regionprops(CC_cell);
cell_shape_props3=regionprops3(CC_cell,'ConvexVolume');
end


 seeds_com=FrameTrackCoordinates_discretized{t}'+1;
 seeds_com=[seeds_com(2,:);seeds_com(1,:);seeds_com(3,:)];
 seeds_com=idx_vec2num(seeds_com,size(track_points))'; %writes it as an index that spans the whole image where each pixel is uniquely numbered

cell_list=NaN(length(seeds_com),34);


for c=1:length(CC_nuc.PixelIdxList)
    if sum(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}))>0 %picks nuc with one track only
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),1)=c;
        blub=zeros(size(track_points));
        blub(CC_nuc.PixelIdxList{1,c})=1;
        [boxes, labels]=imBoundingBox(blub);
%imshow3D(blub(boxes(3)+0.5:boxes(4)-0.5,boxes(1)-0.5:boxes(2)+0.5,boxes(5)-0.5:boxes(6)+0.5))
[ar,eG,el_shape,dG,axis,eV]=get_axis_of_shape_in_box(blub(boxes(3)+0.5:boxes(4)-0.5,boxes(1)+0.5:boxes(2)-0.5,boxes(5)+0.5:boxes(6)-0.5),[1,1,1]);  %its already in microns now
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),2)=ar;
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),3)=el_shape;
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),4)=length(CC_nuc.PixelIdxList{1,c});%*prod(ratio);  %its already in microns now
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),5)=nuc_shape_props(c).Centroid(1);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),6)=nuc_shape_props(c).Centroid(2);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),7)=nuc_shape_props(c).Centroid(3);        
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),8)=boxes(1);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),9)=boxes(2);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),10)=boxes(3);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),11)=boxes(4);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),12)=boxes(5);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),13)=boxes(6);
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),14)= table2array(nuc_shape_props3(c,1));%*prod(ratio); %its already in microns now
        cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),15)= length(CC_nuc.PixelIdxList{1,c})/table2array(nuc_shape_props3(c,1));
      
       %% save the nuclei properties in track structure 
        % find the correct track and time 
        Track_inq=FrameTracks{t}(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}));        
        cell_list_inq=cell_list(ismember(seeds_com, CC_nuc.PixelIdxList{1,c}),1:15); %cell_list_in question, so that for repeated track points issues it assigns the trackpoint info to each track in question since they are supposed to be the same 
        for tr=1:length(Track_inq) %as stated above
        %Tracks{Track_inq}(find(Tracks{Track_inq}(:,1)==t-1),2:4)   %finds correct time in tracks 
        % put the data in 
        Tracks{Track_inq(tr)}(find(Tracks{Track_inq(tr)}(:,1)==t-1),5:19)=cell_list_inq(tr,:);
        end
        % Tracks 1 time in frames starting with 0 
        % Tracks 2 x postition 
        % Tracks 3 y postition 
        % Tracks 4 z postition 
        % Tracks 5 internal control variable 
        % Tracks 6 aspect ratio nuclei
        % Tracks 7 ellipsoid shape parameter
        % Tracks 8 volume in voxel
        % Tracks 9-11 xyz of centroid of nucleus
        % Tracks 12-17 bounding box indices
        % Tracks 18-19 convex volume and convex hull
    end
end
    
if load_cell_outlines==true  
for c=1:length(CC_cell.PixelIdxList)
    if sum(ismember(seeds_com, CC_cell.PixelIdxList{1,c}))>0 %picks cells with one track only
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),16)=c;
        blub2=zeros(size(track_points));
        blub2(CC_cell.PixelIdxList{1,c})=1;
        [boxes, labels]=imBoundingBox(blub2);
%imshow3D(blub(boxes(3)+0.5:boxes(4)-0.5,boxes(1)-0.5:boxes(2)+0.5,boxes(5)-0.5:boxes(6)+0.5))
        [ar,eG,el_shape,dG,axis,eV]=get_axis_of_shape_in_box(blub2(boxes(3)+0.5:boxes(4)-0.5,boxes(1)+0.5:boxes(2)-0.5,boxes(5)+0.5:boxes(6)-0.5),[1,1,1]); %its already in microns now
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),17)=ar;
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),18)=el_shape;
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),19)=length(CC_cell.PixelIdxList{1,c});%*prod(ratio);  %its already in microns now
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),20)=cell_shape_props(c).Centroid(1);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),21)=cell_shape_props(c).Centroid(2);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),22)=cell_shape_props(c).Centroid(3);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),23)=boxes(1);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),24)=boxes(2);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),25)=boxes(3);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),26)=boxes(4);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),27)=boxes(5);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),28)=boxes(6);
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),29)= table2array(cell_shape_props3(c,1));%*prod(ratio); %its already in microns now
        cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),30)= length(CC_cell.PixelIdxList{1,c})/table2array(cell_shape_props3(c,1));

               % find the correct track and time 
               
        Track_inq=FrameTracks{t}(ismember(seeds_com, CC_cell.PixelIdxList{1,c}));    
        cell_list_inq=cell_list(ismember(seeds_com, CC_cell.PixelIdxList{1,c}),16:30);
        %Tracks{Track_inq}(find(Tracks{Track_inq}(:,1)==t-1),2:4)   %finds correct time in tracks 
        % put the data in 
        for tr=1:length(Track_inq)
            Tracks{Track_inq(tr)}(find(Tracks{Track_inq(tr)}(:,1)==t-1),20:34)=cell_list_inq(tr,:);
        end
    end
end
end
%% some analysis and inspection
if false
for tr=1:length(Tracks)
    if max(isinf(Tracks{tr}(:)))
        print(tr)
    end
    
end 
end


if false %test overlay nuclei %colors all nuclei that dont have tracks in them with different color
    test_overlay=zeros([size(track_points),3]);
  	test_overlay(:,:,:,1)=track_points;      
  	nucs_with_track=0*track_points;
  	nucs_wo_track=0*track_points;
    for c=1:length(CC_nuc.PixelIdxList)
        if ismember(c,cell_list(:,1))==0
            nucs_wo_track(CC_nuc.PixelIdxList{1, c})=1;
        else
            nucs_with_track(CC_nuc.PixelIdxList{1, c})=1;
        end
    end
    test_overlay(:,:,:,2)=nucs_with_track;
    test_overlay(:,:,:,3)=nucs_wo_track;
    imshow3D(test_overlay)
    
end

end


