
fname = 'G:\Fred_spidi\2019-07-10 14th tumor MCA 106\scene8\14th tumor MCA 106 5h-post-OP 5min';
cd(fname)

imdir=dir('*tif');

for n = 1:length(imdir)
    t=1+floor((n-1)/21);
    z=n-(t-1)*21;
    l=rgb2gray(double(imread(imdir(n).name))/255);
    %l=imread(fname, n, 'Info', info);
    imorg(:,:,z,t) = l/max(l(:));
end

%imshow3D(imorg(:,:,:,1))  rgb2gray(RGB)

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
    
    for z=1:size(nuc_org,3)
    imwrite(nuc_t_sl(:,:,z), 'nuc_Cropped 3 of 14th-106.tiff','WriteMode','append')
    end

end

%% save nucs with meta data as ome.tif

    
if false
bfsave(nuc, 'nuc_Cropped 3 of 14th-106.ome.ome.tiff');

metadata = createMinimalOMEXMLMetadata(nuc);
pixelSize = ome.units.quantity.Length(java.lang.Double(.05), ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeX(pixelSize, 0.8667);
metadata.setPixelsPhysicalSizeY(pixelSize, 0.8667);
pixelSizeZ = ome.units.quantity.Length(java.lang.Double(.2), ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeZ(pixelSizeZ, 4);
bfsave(plane, 'metadata.ome.tiff', 'metadata', metadata);
end

