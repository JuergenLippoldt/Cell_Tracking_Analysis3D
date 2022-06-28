%% Dont forget to add path

addpath(genpath('C:\Users\sm69w\Documents\MATLAB'))
addpath(genpath('C:\Users\sm69w\Documents\MATLAB\Matlab code'))
%% giving the parameters to the function
pathy='C:\Users\sm69w\Documents\MATLAB\14th tumor s1';
trackname='resliced 3 f200-300 14th tumor MCA 106 5h-post-OP 5min_Tracks2 9 4 15 mean 1 std 1 quality 1 min 1 median 1.1 - 20 2 mean 2.5 std 2.5 quality 2 min 2.5 median 2.6 [um]'; %insert trackfile name here
num_Zslices= 21; %number of z slices here
original_filename = 'f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff' %put the one used for tracking, has to have same zslices as mentioned above
delete_shorter_tracks=true; shorter_than=11; %insert the length of tracks you want omitted from the data.

load_org_images=true;
demix=false;
autocorr=false;
% if contains(trackname,'.xml')
%     error('Pls dont include ".xml" in the trackname')
% % end

%% here it starts
cd(pathy)
[Tracks, ~] = importTrackMateTracks([pathy, '\' ,trackname,'.xml']);

%% some explainations
% most results are stored for each track in the Tracks structure
% Tracks{i} is the ith Track and a (n,x) arreay where n is the number of
% spots in the Tracks                                        % i know its excessive
% Tracks{i}(:,1) is the frame of the spot
% Tracks{i}(:,2) is the x position of the spot
% Tracks{i}(:,3) is the y position of the spot
% Tracks{i}(:,4) is the z position of the spot (0 for 2D)
% Tracks{i}(:,5) is the mean area of the cell and it's neighbors
% Tracks{i}(:,6) T1-trans=>(1 true,0 false, NaN neighbors (dis)appearing)
% Tracks{i}(:,7) shape parameter of Voronoi-tesselation of cell
% Tracks{i}(:,8) is cell squeezing through 2 other cells?
% Tracks{i}(:,9) is cell beeing squeezed through by another cell?
% Tracks{i}(:,10) length of starting squeeze event in that frame
% Tracks{i}(:,11) mean shape parameter of Voronoi-tesselation of cell neighborhood
% Tracks{i}(:,12) # of new neighbors from far away
% Tracks{i}(:,13) what kind of cell population
% Tracks{i}(:,14) *Orderparameter of Normalvectors to centre-vectors
% Tracks{i}(:,15) area of the cell
% Tracks{i}(:,16) at border of large cluster in that time point?
% Tracks{i}(:,17) # of new neighbors from medium far away
% Tracks{i}(:,18) mean shape of immediate neighbourhood
% Tracks{i}(:,19) mean shape of neighbourhood containing 3rd neighbours
% Tracks{i}(:,20) mean shape of cells within 3rd neighbours
% Tracks{i}(:,21) mean shape of cell over the coming time interval
% Tracks{i}(:,22) length of track over the coming time interval
% Tracks{i}(:,23) mean shape of immediate neighbourhood over the coming time interval first time than spacial average
% Tracks{i}(:,24) mean shape of neighbourhood containing 4th neighbours
% Tracks{i}(:,25) mean shape of cells within 4th neighbours
% Tracks{i}(:,26) mean shape of neighbourhood containing 5th neighbours
% Tracks{i}(:,27) mean shape of cells within 5th neighbours
% Tracks{i}(:,28) velocity of cell compared to prior image
% Tracks{i}(:,29) Voronoi cell savely inside frame?
% Tracks{i}(:,30) number of T1's in next hour
% Tracks{i}(:,31) percentage of 'demixed white' inside Voronoi cell
% Tracks{i}(:,32) velocity of cell over the hour around this time point
% Tracks{i}(:,33) angle difference from hexagonal order
% Tracks{i}(:,34) inside? all vertices inside margin?
%%

interval=3*6;  %% important - is system fluid over timeframe of interval frames?


data_in_mueh=false;
%demix=true; % => de_bw is needed to differentiate
edge_length_check = false;
divide_fluid_jammed=false;
ergodic_check=false;

if load_org_images
    fname = 'f200-300 14th tumor MCA 106 5h-post-OP 5min._s1.ome.tiff';
    info = imfinfo(fname);
    num_images = numel(info);
    for k = 1:num_images
        t=floor((k-1)/num_Zslices)+1;
        z=k-(t-1)*num_Zslices;
        im_org(:,:,z,t) = imread(fname, k, 'Info', info);
    end

%     for k = 1:num_images
%         t=floor((k-1)/num_Zslices)+1;
%         z=k-(t-1)*num_Zslices;
%         im_org(:,:,z,t) = rgb2gray(imread(fname, k, 'Info', info));
%     end

end
%% if you want to make it double
%     if class(im_org) == 'uint16'
%         im_org=double(im_org)/max(im_org(:));
%     end
 
%take out tracks shorter than 10 frames
if delete_shorter_tracks
Tracks(cellfun('length',Tracks) <shorter_than)=[];
end

%fill gaps
Tracks=fill_gaps3D(Tracks);


%%

T = Tracks;
% This tells for each frame the tracks of which it is made up
[FrameTracks, FrameTrackCoordinates] = FindAllTracksInFrames(T);
clear T
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );

for f = 1:length(FrameTracks)
    
    % find the vertex coordinates and the cells of the tesselation
    [V,C]	= voronoin( FrameTrackCoordinates{f} );
    
    % find the neighbors of each vertex and each cell of the tesselation
    [VN,CN] = FindNeighbors( C );
    
    % collect data
    Vertices{f}			= V;
    VertexNeighbors{f}	= VN;
    Cells{f}			= C;
    CellNeighbors{f}	= CN;
    
    clearvars C V VN CN
    
    
end

% now we have to associate this to the actual Tracks!
TrackNeighborTracks = cell( size(Tracks) );

TrackNeighborCells = cell( size(Tracks) );
for f=1:length(FrameTracks)
    for t=1:length(FrameTracks{f})
        TrackNeighborCells{FrameTracks{f}(t)}{end+1}	= CellNeighbors{f}{t};
        TrackNeighborTracks{FrameTracks{f}(t)}{end+1}	= FrameTracks{f}(CellNeighbors{f}{t});
    end
end