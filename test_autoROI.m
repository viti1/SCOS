rec = 'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\Vika Head\Fiber400um\Objx20conj_SDS1cm_Gain0_Exp22ms_FR20Hz\';
files = dir([rec '\*.tiff']);

im = double(imread([rec filesep files(1).name]));
my_imagesc(im);

relativeIforCircle = 0.3;

% Find center of mass
x = 1 : size(im, 2); % Columns.
y = 1 : size(im, 1); % Rows.
[X, Y] = meshgrid(x, y);
meanIm = mean(im(:));
centerOfMassX = mean(im(:) .* X(:)) / meanIm;
centerOfMassY = mean(im(:) .* Y(:)) / meanIm;
% figure; imshowpair(CreateCircleMask(size(im),5,centerOfMassY,centerOfMassX),im)

% Find MaxI
maxI = mean(im(CreateCircleMask(size(im),5,centerOfMassY,centerOfMassX)));

% Find Expected Radius
numOfPixelsInSpot = nnz(im(:) > maxI*relativeIforCircle);
approxRadius = sqrt(numOfPixelsInSpot/pi());

[centersBright, radiiBright] = imfindcircles(im > maxI*relativeIforCircle,round(approxRadius*[0.9 1.1]),'ObjectPolarity','bright','Sensitivity',0.99);
my_imagesc(im); viscircles(centersBright, radiiBright,'Color','r','EnhanceVisibility',false,'LineWidth',1);
