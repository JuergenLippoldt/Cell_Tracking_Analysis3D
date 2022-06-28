
fname = 'C:\Users\sm69w\Documents\MATLAB\overlay_14th tumor MCA 106 5h-post-OP 5min.tif';
info = imfinfo(fname);
num_images = numel(info);
for n = 1:num_images-1
%for n=1+50*21:51*21
    t=1+floor((n-1)/21);
    z=n-(t-1)*21;
    l=double(imread(fname, n, 'Info', info));
    imorg(:,:,z,t) = l/max(l(:));
  %  imorg(:,:,z) = l/max(l(:));
end

%imshow3D(imorg(:,:,:,1))

ratio=[1-2/15,1-2/15,4];
disk_radius=round(12*.58/ratio(1));
upper_bound_nuc_area=round(20000*(.58/ratio(1))^2 );

for t=1:size(imorg,4)
   	nuc_org=squeeze(imorg(:,:,:,t));
    for z=1:size(nuc_org,3)
        nuc_t_sl(:,:,z)=nuc_org(:,:,z)-imopen(nuc_org(:,:,z),strel('disk',2*disk_radius));
        nuc_t_sl(:,:,z)=adpmedian(nuc_t_sl(:,:,z), 3); 
    end
    nuc_t_sl = reshape(imadjust(nuc_t_sl(:)/max(nuc_t_sl(:))),size(nuc_t_sl));
    nuc = nuc_t_sl;
    
    int_mask=[];
    for z=1:1:size(nuc,3)
        int_mask(:,:,z)=imclose(nuc(:,:,z),strel('disk',3*disk_radius));
    end
    
    int_mask=int_mask/max(int_mask(:));
    %int_mask = reshape(double(int_mask(:))/max(double(int_mask(:))),size(int_mask));
    int_mask = imdilate(imdilate(int_mask >.4*graythresh(int_mask(int_mask>0)),strel('sphere',1)),strel('disk',3));    
    
    Eim = entropyfilt(nuc);
    %G = fspecial('gaussian',[10 10],3);

    ws=round(disk_radius/2);
    sig_norm=imgaussfilt3(nuc,ws);
    outside2=imerode(sig_norm<.2*graythresh(sig_norm),strel('sphere',1));
    Enmask2=imerode(imdilate(~outside2,strel('disk',3)),strel('disk',1));
    for z=1:size(Enmask2,3)
        Enmask2(:,:,z)=imfill_up2_size(Enmask2(:,:,z),disk_radius^2);
    end
    Enmask=logical(Enmask2.*int_mask);
    
    fgm= nuc>.1*graythresh(nuc(Enmask)) & Enmask;
    for ws=1:10
        for z=1:size(fgm,3)
            fgm(:,:,z) = fgm(:,:,z).*adaptivethreshold(nuc(:,:,z),3*ws^2*disk_radius);%im2bw(nuc(:,:,i),blub);%imbinarize(nuc(:,:,i));
        end
    end
    
    for z=1:size(fgm,3)
    fgm(:,:,z)=imfill_up2_size(fgm(:,:,z),disk_radius^2);
	%fgm(:,:,z) = imfill(fgm(:,:,z),'holes');
    end

    fgm = imerode(fgm, strel('sphere',1));

    
    fgm = imdilate(fgm, strel('sphere',1));
    for z=1:size(fgm,3)
        fgm(:,:,z)=imfill_up2_size(fgm(:,:,z),disk_radius^2);
	%fgm(:,:,z) = imfill(fgm(:,:,z),'holes');
    end
    
    fgm = bwareaopen(fgm,10,6);
    
   	fgm = imerode(fgm, strel('sphere',1));
    fgm = imdilate(fgm, strel('sphere',1));
    bw_nucs(:,:,:,t)=fgm;    
    overlay=0.5*nuc+0.5*fgm;
    for z=1:size(overlay,3)
        imwrite(overlay(:,:,z),'overlay_nuc2_Cropped 3 of 14th-106.tif','WriteMode','append')
    end
end
