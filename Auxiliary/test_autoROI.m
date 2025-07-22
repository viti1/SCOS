% rec = 'C:\Users\tarlevi\OneDrive - Bar-Ilan University - Students\PhD Research\SCOS\Records\Vika Head\Fiber400um\Objx20conj_SDS1cm_Gain0_Exp22ms_FR20Hz\';
rec = 'E:\SCOS\Records\test 7 channels\t1_white_paper';
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
approxRadius = 90; %sqrt(numOfPixelsInSpot/pi());

%%
[centersBright, radiiBright, metric] = imfindcircles(im ,round(approxRadius*[0.95 1.05]),'ObjectPolarity','bright','Sensitivity',0.993); %,'EdgeThreshold',0.3);
my_imagesc(im); viscircles(centersBright, radiiBright,'Color','b','EnhanceVisibility',false,'LineWidth',1);
%%
N = numel(radiiBright);
distances = nan(N);
checked_rows = [];
reducedCenters = nan([0 2]);
reducedRadii = [];
reducedMetric = [];
for i = 1:N    
    if ismember(i,checked_rows)
        continue;
    end
    same_circle_rows = i;
    for j = 1:N
        if j==i; continue; end
        distances(i,j) = norm(diff(centersBright([i,j],:)));
        if distances(i,j) < approxRadius
            same_circle_rows = [same_circle_rows j];
        end
    end 
    checked_rows = [ checked_rows , same_circle_rows ];

    [~, max_ind] =  max( metric(same_circle_rows) ) ; 
    max_row = same_circle_rows(max_ind);
    reducedCenters(end+1,:) = centersBright(max_row,:);
    reducedRadii(end+1)   = radiiBright(max_row);
    reducedMetric(end+1)  = metric(max_row);
end

hold on
my_imagesc(im); viscircles(reducedCenters, reducedRadii,'Color','y','EnhanceVisibility',false,'LineWidth',1);


