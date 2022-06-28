function [im_filled]=imfill_up2_size(im,size)
    im_filled=im;
    CN=bwconncomp(~im,4);
    if CN.NumObjects>0
        for i=1:CN.NumObjects
            if length(CN.PixelIdxList{i})<size
                im_filled(CN.PixelIdxList{i})=1;
            end
        end
    end

end