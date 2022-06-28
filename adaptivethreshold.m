function bw=adaptivethreshold(IM,ws,tm,fac)
%ADAPTIVETHRESHOLD An adaptive thresholding algorithm that seperates the
%foreground from the background with nonuniform illumination.
%  bw=adaptivethreshold(IM,ws,C) outputs a binary image bw with the local 
%   threshold mean-C or median-C to the image IM.
%  ws is the local window size.
%  tm is 0 or 1, a switch between mean and median. tm=0 mean(default); tm=1 median.
%  fac is a factor for grayscale
%
%  Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn)
%  at Tsinghua University, Beijing, China.
%
%  For more information, please see
%  http://homepages.inf.ed.ac.uk/rbf/HIPR2/adpthrsh.htm

if (nargin<2)
    error('You must provide the image IM, the window size ws, and C.');
elseif (nargin==2)
    tm=0;
    fac=0.8;
elseif (nargin==3) && (tm==0 || tm==1)
    fac=0.8;
elseif (tm~=0 && tm~=1)
    error('tm must be 0 or 1.');
end

IM=mat2gray(IM);

if tm==0
    if length(size(IM))==2
        mIM=averagefilter(IM, [ws ws], 'replicate');%imfilter(IM,fspecial('average',ws),'replicate');
    elseif length(size(IM))==3 && size(IM,3)>3
        h = fspecial3('average',ws);
        mIM=imfilter(IM,h, 'replicate');
        %mIM=fspecial(IM,ws(1)/4,'FilterSize',2*ceil(ws./4)+1);    
    end
else
    mIM=medfilt2(IM,[ws ws]);
end
sIM=mIM-IM; 
thresh=fac*graythresh(sIM); 
bw=sIM > thresh;
bw=logical(-(bw-1));