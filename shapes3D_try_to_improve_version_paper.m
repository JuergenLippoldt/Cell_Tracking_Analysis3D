disk_radius=12;             % radius of a disk used to blur. The number of pixels should roughly equal 10�m. (whole number)
upper_bound_nuc_area=20000; % nuclei are smaller than .. in square pixel - should rather be to high. (whole number)
smoothy=1;                  % how often do you want to delate and erode to smooth the resulting grid (whole number)
lower_bound_cell_area=40;   % cells are larger than .. in square pixel. (whole number)
shrink_bound=true;          % do you want to shrink the boundary? will take 1-2 hours
multi_knuedel=1;            % how many kn�del in one go? (overnight) 
get_shape_isometry=false;  
lify=true;
insideholes=false;
pathy='E:\10A_13 monolayer\10A_13 day3 good 40X';


if lify
cd(pathy)
name_source_file='10A_13 day3 good 40X.lei';
org=bfopen(name_source_file);
end


if ~lify
    cd(pathy)
org=dir('*.mat');
end

% ratio=[0.608, 0.608, 0.503*.87];
%     % the scale
%     zScalingFactor=.87;  %% average scale factor
%     Scale_Barbara = [0.608, 0.608, 0.503*zScalingFactor]; % not EXACTLY known
%     Scale_LSM = [0.586, 0.586, 0.644*zScalingFactor];    
%     Scale_Mierke =[0.445, 0.445, 0.363*zScalingFactor];
%     Scale_Mierke_63x =[0.308 ,0.308 , 0.46*zScalingFactor];
% Voxel-Width	[�m]		0.585437 
% Voxel-Height	[�m]		0.585437 
% Voxel-Depth	[�m]		0.644273 
h = waitbar(0,'Spheroid Segmentation');
if lify
stacks=[];
for kn=1:length(org) 
    if length(org{kn,1})>50
    stacks=[stacks, kn];
    end
    
end
else 
    stacks=1:length(org);
end


for kn=stacks  %% kn=1; 
   clearvars -except pathy insideholes stacks lify maty h name_source_file get_shape_isometry ratio org kn disk_radius upper_bound_nuc_area smoothy lower_bound_cell_area shrink_bound multi_knuedel Scale_Barbara Scale_LSM Scale_Mierke knsleft
%  
%load('original_spheroid_data.mat')
if lify
  name= org{kn,2}.values.toArray;
  keys= org{kn,2}.keySet.toArray;
  keys_logical=zeros(size(1,length(keys)));
  for fu=1:length(keys)
    keys_logical(fu)=strcmp(keys(fu),'Image name') || strcmp(keys(fu),'Series name');
  end
  name=name(keys_logical==1);
  namestart=strfind(name,'/');
  if ~isempty(namestart)
  name=name(namestart+1:end);
  end
% %   
   act_org=zeros(size(org{kn,1}{1,1},1),size(org{kn,1}{1,1},2),length(org{kn,1})/2);
   nuc_org=zeros(size(org{kn,1}{1,1},1),size(org{kn,1}{1,1},2),length(org{kn,1})/2);
% 
for k=1:1:size(nuc_org,3)
    nuc_org(:,:,k)=org{kn,1}{2*k,1};%imadjust(org{1,1}{2*k,1});
    act_org(:,:,k)=org{kn,1}{2*k-1,1};%imadjust(org{1,1}{2*k-1,1});
    nuc(:,:,k)=nuc_org(:,:,k)-imopen(nuc_org(:,:,k),strel('disk',disk_radius));
    act(:,:,k)=act_org(:,:,k)-imopen(act_org(:,:,k),strel('disk',disk_radius));
    act(:,:,k) = adpmedian(act(:,:,k), 3); 
    nuc(:,:,k) = adpmedian(nuc(:,:,k), 3); 
 %   nuc(:,:,k) = wiener2(nuc(:,:,k));
end
nuc = reshape(imadjust(nuc(:)/max(nuc(:))),size(nuc));
act = reshape(imadjust(act(:)/max(act(:))),size(act));
% clearvars org
% 
% pathy='\\pwm039\F_Knoedel2000\KnoedelEdge STED LSM\2016-08-17 Knoedel SiRDNA Ph488 10A 436\Knoedelshape';
% name='Knoedelshape_10A-d1-1';

if lify && contains('*lif',name_source_file)
omeMeta = org{kn,4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

voxelSizeXdefaultValue = double(omeMeta.getPixelsPhysicalSizeX(0).value());    % returns value in default unit
voxelSizeYdefaultValue = double(omeMeta.getPixelsPhysicalSizeY(0).value());    % returns value in default unit
voxelSizeZdefaultValue = double( omeMeta.getPixelsPhysicalSizeZ(0).value());    % returns value in default unit

ratio=[voxelSizeXdefaultValue, voxelSizeYdefaultValue, voxelSizeZdefaultValue];
else
    ratio=[.644,.644,.644];
end

disk_radius=round(12*.58/ratio(1));
upper_bound_nuc_area=round(20000*(.58/ratio(1))^2 );
end
if ~lify
    load(org(kn).name)
    disk_radius=round(12*.58/ratio(1) );
    upper_bound_nuc_area=round(20000*(.58/upper_bound_nuc_area)^2 );
    act=double(act)/max(double(act(:)));
    nuc=double(nuc)/max(double(nuc(:)));    
    name=org(kn).name;
end
%      if ~isempty(strfind(name,'S001'))
%          ratio = Scale_Barbara;
%      elseif ~isempty(strfind(name,'-M-'))
%          ratio = Scale_Mierke;
%      else
%          ratio = Scale_LSM;
%      end 
% oldp=cd(pathy);
%  act_dir=dir([name,'*ch00*']);
%  nuc_dir=dir([name,'*ch01*']);
% %org_overlay_dir=dir('*.tif');
% cd(oldp);
% 
%  for k=1:1:length(nuc_dir)
%      nuc_org(:,:,k)=imread([pathy,'\',nuc_dir(k).name]);
%      act_org(:,:,k)=imread([pathy,'\',act_dir(k).name]);
%      nuc(:,:,k)=nuc_org(:,:,k)-imopen(nuc_org(:,:,k),strel('disk',3*disk_radius));
%      act(:,:,k)=act_org(:,:,k)-imopen(act_org(:,:,k),strel('disk',3*disk_radius));
%      act(:,:,k) = adpmedian(act(:,:,k), 3); 
%      nuc(:,:,k) = adpmedian(nuc(:,:,k), 3); 
% end
% nuc = reshape(imadjust(double(nuc(:))/max(double(nuc(:)))),size(nuc));
% act = reshape(imadjust(double(act(:))/max(double(act(:)))),size(act));
if lify
int_mask=[];
for k=1:1:size(act,3)
    int_mask(:,:,k)=imclose(.5*act(:,:,k)+.5*nuc(:,:,k),strel('disk',3*disk_radius));
end
int_mask=int_mask/max(int_mask(:));
%int_mask = reshape(double(int_mask(:))/max(double(int_mask(:))),size(int_mask));
int_mask = imdilate(imdilate(int_mask >.4*graythresh(int_mask(int_mask>0)),strel('sphere',1)),strel('disk',3));


Eim = entropyfilt(act);
%G = fspecial('gaussian',[10 10],3);
Enmask=zeros(size(Eim));
for z=1:size(Eim,3)
    Eim(:,:,z)=Eim(:,:,z)/max(max(Eim(:,:,z)));
Enmask(:,:,z) = im2bw(Eim(:,:,z),graythresh(Eim(mean(Eim(:,:,z))>0.1,mean(Eim(:,:,z),2)>0.1,z))); 
Enmask(:,:,z)=imfill(Enmask(:,:,z));
holes=imfill(Enmask(:,:,z))-Enmask(:,:,z);
Enmask(:,:,z)=Enmask(:,:,z) + holes;
end
Enmask=Enmask.*int_mask;
% imshow3D(BW1)

%% adjust actin-signal loss in z direction
signalstrength=(act.*Enmask);
signalstrength(signalstrength==0)=NaN;
signalstrength=squeeze(nanmean(squeeze(nanmean(signalstrength,1)),1));
for z=1:size(act,3)
    if ~isnan(signalstrength(z))==1
        act(:,:,z)=act(:,:,z)*nanmax(signalstrength)/max(signalstrength(z),.5*nanmax(signalstrength));
    end
end
act=act/max(act(:))*3/2;  % 3/2 is arbitrary ... just makes the image slighte brighter and 'overexposes' brigth spots 
act(act>1)=1;
end

%% entropy
ws=round(disk_radius/4);

sig_norm=imgaussfilt3(nuc+act,ws);
%sig_norm=sig_norm/max(sig_norm(:));
if insideholes
    outside2=imerode(sig_norm<.8*graythresh(sig_norm),strel('sphere',1));
else
    outside2=imerode(sig_norm<.2*graythresh(sig_norm),strel('sphere',1));
end

Enmask=imerode(imdilate(~outside2,strel('disk',3)),strel('disk',1));
for z=1:size(Enmask,3)
    Enmask(:,:,z)=imfill_up2_size(Enmask(:,:,z),disk_radius^2);
end

%% lets test something 

%     
% IM=imgaussfilt3(nuc,2);
% IM=IM-nuc; 
% En=entropyfilt(IM);
%     lehisto=imhist(reshape(En,[size(nuc,1)*size(nuc,3),size(nuc,2)]),100);
%     level=otsuthresh(lehisto);
% seeds=En < level;
% seeds=seeds .*Enmask;

%% find foreground
% fgm=zeros(size(nuc));
% 
% blub=.1*graythresh(nuc);
% for f=1:5
% for i=1:size(fgm,3)
%     fgm(:,:,i) = adaptivethreshold(nuc(:,:,i),f^2*disk_radius);%im2bw(nuc(:,:,i),blub);%imbinarize(nuc(:,:,i));
%     %fgm(:,:,i) = bradley(nuc(:,:,i),[35 35], 0);
% end
% end
% 
% fgm=fgm.*Enmask.*(nuc>blub);
% 
% fgm = imdilate(fgm, strel('sphere',1));
% fgm = imfill(fgm,'holes');
% fgm = imerode(fgm, strel('sphere',4));
% %fgm = imdilate(fgm, strel('sphere',1));
% CN=bwconncomp(fgm,6);
% fgm = zeros(size(nuc));
% for i=1:length(CN.PixelIdxList)
%     if length(CN.PixelIdxList{i})>500
%         fgm(CN.PixelIdxList{i})=1;
%     end
% end
% %fgm = imdilate(fgm, strel('sphere',2));
% for i=1:3
%     disty=bwdist(~fgm);
%     fgm=disty>1;
% end
% 
% %fgm=fgm.*BW1;
% CN=bwconncomp(fgm,6);
% disty=bwdist(~fgm);
% disty=imgaussfilt3(disty,3);
% 
% %     lehisto=imhist(reshape(nuc,[size(nuc,1)*size(nuc,3),size(nuc,2)]),100);
% %     level=otsuthresh(lehisto);
% %     bright_mask=nuc > level;
% %     bright_mask=imdilate(bright_mask,strel('sphere',1));
%      
% %# s = 3D array
% msk = true(3,3,3);
% msk(2,2,2) = false;
% %# assign, to every voxel, the maximum of its neighbors
% maxy = imdilate(disty,msk);
% seeds = disty > maxy;
% %seeds = seeds .* bright_mask;
% seeds = imdilate(seeds, strel('sphere',4));
% blub=imgaussfilt3((act-imgaussfilt3(act,ratio)),.1*ratio);
% for i=1:size(blub,3)
%     blub(:,:,i)=adpmedian(blub(:,:,i),9);
% end
% blub=reshape(imadjust(blub(:)/max(blub(:))),size(blub));
% test=imgaussfilt3(nuc,.1*ratio)-blub;%test=nuc+blub1-1*blub2;
% test(test<0)=0;
% nuc2=reshape(imadjust(test(:)/max(test(:))),size(test));
% 
% %nuc2=imgaussfilt3(nuc,3/10*ratio)-abs(act-imgaussfilt3(act,2*ratio));
% nuc2(nuc2<0)=0;
% nuc2=reshape(imadjust(nuc2(:)/max(nuc2(:))),size(nuc2));

fgm_act= act>.7*graythresh(act) & Enmask;
for ws=1:12
    for z=1:size(fgm_act,3)
        fgm_act(:,:,z) = fgm_act(:,:,z).*adaptivethreshold(act(:,:,z),2*ws^2*disk_radius);%im2bw(nuc(:,:,i),blub);%imbinarize(nuc(:,:,i));
    %fgm(:,:,i) = bradley(nuc(:,:,i),[35 35], 0);
    end
end

nuc2=nuc-2*act.*fgm_act;
nuc2(nuc2<0)=0;
nuc2=nuc2/max(nuc2(:));
fgm= nuc2>.6*graythresh(nuc2) & Enmask;
for ws=1:10
    for z=1:size(fgm,3)
        fgm(:,:,z) = fgm(:,:,z).*adaptivethreshold(nuc2(:,:,z),2*ws^2*disk_radius);%im2bw(nuc(:,:,i),blub);%imbinarize(nuc(:,:,i));
    %fgm(:,:,i) = bradley(nuc(:,:,i),[35 35], 0);
    end
end

if false  %% important for 436 spheroids but slightly disturbs other segmentations 
border=fgm & ~ imerode(imerode(fgm, strel('sphere',1)),strel('disk',4)) & ~imerode(Enmask,strel('sphere',12));
border_rest=fgm &  imerode(imerode(fgm, strel('sphere',1)),strel('disk',4)) & ~imerode(Enmask,strel('sphere',12));

border_rest2=imdilate(imdilate(border_rest, strel('sphere',1)),strel('disk',2)) & fgm;
fgm=logical(fgm - border + border_rest2);
end

%fgm=fgm.*(nuc>blub);
%fgm = imdilate(fgm, strel('sphere',1));
for z=1:size(fgm,3)
    fgm(:,:,z)=imfill_up2_size(fgm(:,:,z),disk_radius^2);
	%fgm(:,:,z) = imfill(fgm(:,:,z),'holes');
end

fgm = imerode(fgm, strel('sphere',2));
fgm = imdilate(fgm, strel('sphere',2));
for z=1:size(fgm,3)
    fgm(:,:,z)=imfill_up2_size(fgm(:,:,z),disk_radius^2);
	%fgm(:,:,z) = imfill(fgm(:,:,z),'holes');
end

fgm = imerode(fgm, strel('sphere',1));
%fgm = imdilate(fgm, strel('sphere',1));
CN=bwconncomp(fgm,6);
fgm = zeros(size(nuc));
for i=1:length(CN.PixelIdxList)
    if length(CN.PixelIdxList{i})>.2*disk_radius^2
        fgm(CN.PixelIdxList{i})=1;
    end
end

fgm = imdilate(fgm, strel('sphere',1));
fgm = imerode(fgm, strel('sphere',1));
fgm = imdilate(fgm, strel('sphere',1));
fgm = imerode(fgm, strel('sphere',1));
% %fgm = imdilate(fgm, strel('sphere',2));
% for i=1:3
%     disty=bwdist(~fgm);
%     fgm=disty>1;
% end

%fgm=fgm.*BW1;

disty_start=bwdistsc(~fgm,1./ratio);
disty=imgaussfilt3(disty_start,2*ratio);

%     lehisto=imhist(reshape(nuc,[size(nuc,1)*size(nuc,3),size(nuc,2)]),100);
%     level=otsuthresh(lehisto);
%     bright_mask=nuc > level;
%     bright_mask=imdilate(bright_mask,strel('sphere',1));
     
%# s = 3D array
msk = true(3,3,3);
msk(2,2,2) = false;
%# assign, to every voxel, the maximum of its neighbors
maxy = imdilate(disty,msk);
seeds = disty >= maxy & disty>0;
%seeds = seeds .* bright_mask;
seeds = imdilate(seeds, strel('sphere',2));

%clear disty_start maxy disty CN border_rest border
%% testing nuc segmentation

%% nucsy

% [Gx_nuc,Gy_nuc,Gz_nuc] = imgradientxyz(nuc);
% 
% gradmag_nuc = Gx_nuc.^2 + Gy_nuc.^2 + Gz_nuc.^2;% + Gx_nuc.^2 + Gy_nuc.^2 + Gz_nuc.^2;
% for k=1:size(gradmag_nuc,3)
%     gradmag_nuc(:,:,k)=gradmag_nuc(:,:,k)-imopen(gradmag_nuc(:,:,k),strel('disk',5));
% end
% gradmag_nuc = imgaussfilt3(gradmag_nuc,1);
% gradmag_nuc = reshape(imadjust(gradmag_nuc(:)/max(gradmag_nuc(:))),size(gradmag_nuc));
% 
% ws=25;
% 
% nuc_norm=double(imgaussfilt3(nuc,ws));
% nuc_norm=nuc_norm/max(nuc_norm(:));
% 
% gradmag_nuc=gradmag_nuc.*nuc_norm;
% gradmag_nuc = gradmag_nuc/max(gradmag_nuc(:));
% 
% lehisto=imhist(reshape(gradmag_nuc,[size(nuc,1)*size(nuc,3),size(nuc,2)]),100);
% level=otsuthresh(lehisto);
% test=gradmag_nuc > level;
%imshow3D(test);

[seeds_com,seeds2]=find_nuc_com(seeds);
%seeds2=imdilate(seeds2,strel('sphere',1));
seeds_com=idx_vec2num(seeds_com,size(nuc));
%[~,Enmask2]=find_nuc_com(~Enmask);
% 


% IM=nuc-imgaussfilt3(nuc,10)+act-imgaussfilt3(act,10); 
% 
%      lehisto=imhist(reshape(IM,[size(nuc,1)*size(nuc,3),size(nuc,2)]),100);
%      level=otsuthresh(lehisto);
% outside = IM<level;
% outside = ~bwareaopen(~outside, 25);
% outside = bwareaopen(outside, 25);
% 
% outside3 = ~imerode(imdilate(imdilate(~outside+seeds,strel('sphere',6)),strel('disk',8)),strel('disk',7));
% 
% for z=1:size(outside2,3)
% outside3(:,:,z)=~imfill_up2_size(~outside2(:,:,z),1000);
% end

% 
% gradmag2 = imimposemin(gradmag_nuc,outside2 | seeds2); % use this if you estimate the background
% 
% L = watershed(gradmag2);
% I = zeros(size(gradmag2));
% I(imdilate(L == 0, ones(3, 3))) = 1;
% 
% outlines= I > 0.5;
% 
% CC_shapes=bwconncomp(~outlines);
% shapes_nuc=zeros(size(nuc));
% for c=2:length(CC_shapes.PixelIdxList)
%     if ismember(1,ismember(seeds_com, CC_shapes.PixelIdxList{1,c}))==1 && length(CC_shapes.PixelIdxList{1,c}) < 5525000 && length(CC_shapes.PixelIdxList{1,c}) > 500
%         shapes_nuc(CC_shapes.PixelIdxList{1,c})=1;
%     end
% end
% 
% [~,seeds3]=find_nuc_com(shapes_nuc);
%imshow3D(shapes_nuc)



if false % just to copy + paste for debugging purpous
    overlay_nuc=zeros(size(nuc,1),size(nuc,2),size(nuc,3),3);
    overlay_nuc(:,:,:,1)=act;%~imerode(Enmask,strel('sphere',6));
    overlay_nuc(:,:,:,2)=fgm;%seeds3;
    overlay_nuc(:,:,:,3)=nuc;
    imshow3D(overlay_nuc);
end


%% act

[Gx,Gy,Gz] = imgradientxyz(act);

gradmag_act = Gx.^2 + Gy.^2 + Gz.^2;% + Gx_nuc.^2 + Gy_nuc.^2 + Gz_nuc.^2;
for k=1:size(gradmag_act,3)
    gradmag_act(:,:,k)=gradmag_act(:,:,k)-imopen(gradmag_act(:,:,k),strel('disk',5));
%    gradmag_act(:,:,k)= imgaussfilt(gradmag_act(:,:,k),1);
end
gradmag_act = imgaussfilt3(gradmag_act,0.4*1./ratio);
gradmag_act = reshape(imadjust(gradmag_act(:)/max(gradmag_act(:))),size(gradmag_act));

gradmag_act(gradmag_act<0)=0;
gradmag_act=(gradmag_act+.15).*(act+.4)-.1*nuc; 
% for k=1:size(gradmag_act,3)
% gradmag_act(:,:,k)=(gradmag_act(:,:,k)+.2).*(act(:,:,k)+.5)-.2*nuc(:,:,k);  %     gradmag = (imgaussfilt(gradmag,1.5)+.2).*(act+.5);
% end
%gradmag_act=act-.1*nuc;
%gradmag_act(gradmag_act<0)=0;
%gradmag_act=gradmag_act-min(gradmag_act(:));
gradmag_act=gradmag_act/max(gradmag_act(:));
%figure; imshow3D(gradmag_act);


clear Gx Gy Gz

gradmag2 = imimposemin(gradmag_act, ~Enmask|seeds2);
L = watershed(gradmag2);
I3 = zeros(size(gradmag2));
I3(imdilate(L == 0, ones(3, 3))) = 1;

outlines=zeros(size(I3));
outlines_smooth=zeros(size(I3));
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
end
% 

holes=ones(size(outlines_smooth))-outlines_smooth;
bigholes = bwareaopen(holes, lower_bound_cell_area);
smallholes = holes & ~bigholes;
outlines_smooth=outlines_smooth | smallholes;
% outlines_smooth=bwmorph(outlines_smooth,'shrink',Inf);
% outlines_smooth=bwmorph(outlines_smooth,'dilate');
% outlines_smooth=bwmorph(outlines_smooth,'erode');
%figure; imshow3D(outlines_smooth)

holes=ones(size(outlines_smooth))-(outlines_smooth &Enmask);
if shrink_bound
    outlines_shrink=shrink_edge3D(outlines_smooth);
    holes=ones(size(outlines_shrink))-(outlines_shrink &Enmask);
end

CC=bwconncomp(holes,6);
clear c_size
for c=1:length(CC.PixelIdxList)
    c_size(c)=length(CC.PixelIdxList{1,c});
end
CC.PixelIdxList(c_size==max(c_size))=[];
c_size(c_size==max(c_size))=[];
cells=CC.PixelIdxList(c_size >100 & c_size < 25000);%& c_size <5000
shapes=zeros(size(outlines_shrink));
shapes_labeled=zeros(size(outlines_shrink));
for c=1:length(cells)
    shapes(cells{1,c})=1;
    shapes_labeled(cells{1,c})=c;
end

if get_shape_isometry
    k=struct;
    k.shapes=shapes;
    k.Scale=ratio;
    k.Description=name;
	% find the regions
	k.Regions		= bwconncomp( k.shapes, 6 );

	% remove die äußere (die größte) region weil sie kein ZELLE ist sondern
	% nur der backgound und den wollen wir uns ab sofort vom Hals halten
	Sizes			= cellfun( @numel, k.Regions.PixelIdxList );
	[~, Outside]	= max(Sizes);
	k.Regions.PixelIdxList	= {k.Regions.PixelIdxList{[1:Outside-1,Outside+1:end]}};
	k.Regions.NumObjects	= k.Regions.NumObjects-1;
	
	% the label matrix
	k.LabelMatrix			= labelmatrix(k.Regions);

	% find region centroids and volume
	k.Regions.Props	= regionprops( k.Regions, 'Centroid', 'Area', 'BoundingBox' );
    k = Knoedel_FindKnoedelCenter( k );
    k = Knoedel_Filter( k ); 
    k = Knoedel_ComputeShapeTensor3D( k );
    Centroids		= vertcat(k.Regions.Props.Centroid);
    relCentroids	= Centroids - repmat( k.KnoedelDetection.FitCenter, size(Centroids, 1), 1 );
    relCentroids	= relCentroids .* repmat( k.Scale, size(Centroids, 1), 1 );
    distances		= norm3d( relCentroids );

    good = intersect( k.KnoedelDetection.RemainingParticles, k.NotTooSmall);
    
    d_good=distances(good);
    s_good=k.ShapeTensor.ElShape(good)';
%     [~,bin]=histc(distances(good),linspace(min(distances(good)),max(distances(good)),10));
%     ploty=[];
%     for i=1:max(bin)-1
%         ploty(i,:)=[min(d_good)+(i-.5)/10*(max(d_good)-min(d_good)),median(s_good(bin==i)),std(s_good(bin==i)),sum(bin==i)];  
%     end
% figure;
% plot(ploty(:,1),ploty(:,2))

end
if false
    [Knoedel_com,Knoedel_radius]=estimate_Knoedel_centre(Enmask,ratio);  %% Knoedel_radius in �m
    [boxes, labels]=imBoundingBox(shapes_labeled);
    h = waitbar(0,'3D aspect ratio');
    for c=1:length(cells)
        particular_cell=shapes_labeled==c;
        s = regionprops(particular_cell,'centroid');
        cells{5,c}=s.Centroid;
        particular_cell=particular_cell(max(boxes(c,3)-.5,1):min(boxes(c,4)+.5,size(particular_cell,1)),max(boxes(c,1)-.5,1):min(boxes(c,2)+.5,size(particular_cell,2)),max(boxes(c,5)-.5,1):min(boxes(c,6)+.5,size(particular_cell,3)));
        [ar,dG,p_el_shape]=get_axis_of_shape_in_box(particular_cell,ratio);
        cells{2,c}=boxes;
        cells{3,c}=dG;
        cells{4,c}=ar;   
        cells{7,c}=p_el_shape;
        waitbar(c/length(cells),h);
        ars(1,c)=ar;
        cells{6,c}=sum(((cells{5,c}-Knoedel_com').*ratio).^2)^.5;
    end
    
    shape_for_distance=[];
    for c=1:length(cells)
        shape_for_distance(1,c)=sum(((cells{5,c}-Knoedel_com').*ratio).^2)^.5;
        shape_for_distance(2,c)=cells{4,c};
        shape_for_distance(3,c)=cells{7,c};
    end
[~,bin]=histc(shape_for_distance(1,:),linspace(min(shape_for_distance(1,:)),max(shape_for_distance(1,:)),20));
for i=1:max(bin)-1
   ploty(i,:)=[min(shape_for_distance(1,:))+(i-.5)/20*(max(shape_for_distance(1,:))-min(shape_for_distance(1,:))),median(shape_for_distance(3,bin==i)),std(shape_for_distance(3,bin==i)),sum(bin==i)];  
end
figure;
plot(ploty(:,1),ploty(:,2))
end


% shapes=shapes.*bwdist(~Enmask)>20;
% shapes = imdilate(shapes, strel('sphere',3));
% shapes = imerode(shapes, strel('sphere',3));
% holes=holes+shapes;

%shapes=zeros(size(shapes));

if false
overlay=zeros(size(act,1),size(act,2),size(act,3),3);
%overlay_nuc=zeros(size(act,1),size(act,2),size(act,3),3);
overlay(:,:,:,1)=double(imdilate(outlines_shrink,strel('disk',1)));
overlay(:,:,:,2)=double(act)/255;
overlay(:,:,:,3)=double(nuc)/255;
figure;
imshow3D(overlay);
end
% if multi_knuedel>1 && shrink_bound
%     overlay(:,:,:,1)=outlines_shrink;

[outlines_shrink2,seeds_sorted,shapes_notsmall]=seed_selection(act,nuc,Enmask,outlines_shrink,seeds,ratio);

nuc_d=nuc;
act_d=act;
nuc=uint8(round(255*nuc));
act=uint8(round(255*act));
nuc_org=uint8(round(255*nuc_org));
act_org=uint8(round(255*act_org));
outlines=logical(outlines);
Enmask=logical(Enmask);
shapes=logical(shapes);
seeds_sorted=logical(seeds_sorted);
shapes_notsmall=logical(shapes_notsmall);
if isstruct(k)
    save([pathy, name, '_maybe_improved.mat'],'nuc_org','act_org','nuc','act','outlines','outlines_smooth','outlines_shrink','Enmask','seeds','shapes','cells','k','ratio','outlines_shrink2','seeds_sorted','shapes_notsmall'); 
else
    save([pathy, name, '_maybe_improved.mat'],'nuc_org','act_org','nuc','act','outlines','outlines_smooth','outlines_shrink','Enmask','seeds','shapes','cells','ratio','outlines_shrink2','seeds_sorted','shapes_notsmall'); 
end
% end
%end


% mov = immovie(permute(overlay,[1 2 4 3]));
% movie2avi(mov, '10A_shapes.avi');
% for i=1:size(overlay,3)
% %imwrite(squeeze(overlay(:,:,i,:)),['overlay_10A_',num2str(i),'.png']);
% imwrite(squeeze(overlay(:, :, i,:)), 'overlay_stack.tiff', 'WriteMode', 'append');
% end

clear shape_for_distance  disty_start Eim holes I I3 IM Ix Iy L lehisto maxy nuc_norm

if false
shapes_r=zeros(size(act,1),size(act,2),size(act,3));
shapes_g=zeros(size(act,1),size(act,2),size(act,3));
shapes_b=zeros(size(act,1),size(act,2),size(act,3));
for c=1:length(cells)
    if mod(c,12)==0 || mod(c,12)==1 || mod(c,12)==2 || mod(c,12)==11 || mod(c,12)==10
    shapes_b(cells{c})=1;
    end
    if mod(c,12)==4 || mod(c,12)==9 
    shapes_b(cells{c})=0.5;
    end
    if mod(c,12)==2 || mod(c,12)==3 || mod(c,12)==4 || mod(c,12)==5 || mod(c,12)==6
    shapes_r(cells{c})=1;
    end
    if mod(c,12)==1 || mod(c,12)==7 
    shapes_r(cells{c})=0.5;
    end
    if mod(c,12)==6 || mod(c,12)==7 || mod(c,12)==8 || mod(c,12)==9 || mod(c,12)==10
    shapes_g(cells{c})=1;
    end
    if mod(c,12)==5 || mod(c,12)==11 
    shapes_g(cells{c})=0.5;
    end
end
 shapesrgb=zeros(size(act,1),size(act,2),size(act,3),3);
 shapesrgb(:,:,:,1)=shapes_r;
 shapesrgb(:,:,:,2)=shapes_g;
 shapesrgb(:,:,:,3)=shapes_b;
end

waitbar(kn/length(org));
end
close(h) 
%save([name,'.mat'],'overlay','outlines','outlines_shrink','outlines_smooth', 'seeds','shapes','shapes_nuc','-v7.3')

% setpref('Internet','SMTP_Server','server1.rz.uni-leipzig.de');
% setpref('Internet','E_mail','juergen.lippoldt@uni-leipzig.de');
% sendmail('juergen.lippoldt@uni-leipzig.de','Matlab done')



%% ploty
if false 
   % plane 230 
    overlay=zeros(size(act,1),size(act,2),size(act,3),3);
%overlay_nuc=zeros(size(act,1),size(act,2),size(act,3),3);
overlay(:,:,:,1)=double(imdilate(outlines_shrink2,strel('disk',1)));
overlay(:,:,:,2)=double(act)/255;
overlay(:,:,:,3)=double(nuc)/255;
figure;
imshow3D(overlay);
overlay_seg=squeeze(overlay(:,:,230,:));
imwrite(overlay_seg,'start.tif')
    
end
