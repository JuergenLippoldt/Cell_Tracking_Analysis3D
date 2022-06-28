function plot_tracked_nuc(track,shape_file,save_folder,save_name,z_dim,tr_num)
x_umperpixel=0.867; %when tracks are in um then 0.867
z_umperpixel=4; 
ratio=[x_umperpixel,x_umperpixel,z_umperpixel];
track(:,1:4)=track(:,1:4)+1;  % fiji starts with counting at 0
%boxs dimensions §
l_y_edge=track(:,14); %lower y edge
u_y_edge=track(:,15); %upper y edge
l_x_edge=track(:,12); 
u_x_edge=track(:,13);
l_z_edge=track(:,16);
u_z_edge=track(:,17);

box_boundaries=[max(min(l_y_edge(l_y_edge>0))-0.5,1),max(u_y_edge(u_y_edge>0))-0.5,max(min(l_x_edge(l_x_edge>0))-0.5,1),max(u_x_edge(u_x_edge>0))+0.5, max(min(l_z_edge(l_z_edge>0))-0.5,1),max(u_z_edge(u_z_edge>0))-0.5,track(1,1),track(end,1)]; %setting up the boundries of the box relative to each track, that indicates the nuclei in question has to be closest to the boundries!

image_box=zeros(1+box_boundaries(2)-box_boundaries(1),1+box_boundaries(4)-box_boundaries(3),1+box_boundaries(6)-box_boundaries(5),box_boundaries(8)-box_boundaries(7)); %4D segemnent of the zoomed in area where the nuclei in question lies

%imshow3D(boxes(3)+0.5:boxes(4)-0.5,boxes(1)-0.5:boxes(2)+0.5,boxes(5)-0.5:boxes(6)+0.5)

for t=1:size(track,1)%track(1,1)+1:track(end,1)+1 %fix for starting at 0 case
    trt=track(t,1)-1%t-track(1,1); %local track timeframe
    
	info = imfinfo(shape_file);
    num_images = numel(info);
    for n = 1:num_images-1
    t_nuc=1+floor((n-1)/z_dim);
        if t_nuc==trt+1%t
            z=n-(t_nuc-1)*z_dim;
            shape_org(:,:,z)=double(imread(shape_file, n, 'Info', info)); %change ratio of picture here to make sure its in the micron dimensions
        end
    end
    
    shape_box=shape_org(box_boundaries(1):box_boundaries(2),box_boundaries(3):box_boundaries(4),box_boundaries(5):box_boundaries(6)); %zooming in into the area of interest
    CC_nuc=bwconncomp(shape_box,6);
    
    %indexing the nuclei
   	seed_com= round(track(t,2:4));  %rounding xyz for this track at the local track time  %t was trt before
    seed_com=[seed_com(2);seed_com(1);seed_com(3)];
    seed_com=seed_com-box_boundaries(1:2:5)'; %since we are zooming in on one area, the pixel index will change and we have to account for that
    seed_com=idx_vec2num(seed_com,size(shape_box)); %indexing the nuceli so that we know where it is in refrence to the pixel index list
        

    shape_box_tr=shape_box;


    for c=1:length( CC_nuc.PixelIdxList) %go through nucs and find the one that is tracked
        if sum(ismember(seed_com, CC_nuc.PixelIdxList{1,c}))==1
        % make it red
 
            shape_box_tr(CC_nuc.PixelIdxList{1,c})=75; 
        end
    end


    
    
    image_box(:,:,:,trt+1)=shape_box_tr; 
end


    maximum = max(image_box(:));
    minimum = min(image_box(:));
    image_box = (image_box - minimum)./(maximum - minimum);
    range = 0;
%imshow4(image_box,1,1,range,'hot')
save([save_folder, save_name,'-track ', num2str(tr_num),' t ',num2str(track(1,1)),'-',num2str(track(end,1)),' nuc_movie.mat'],'image_box');


end