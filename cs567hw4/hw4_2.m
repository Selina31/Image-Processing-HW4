% HW4 Q2

% The feature I chose is the largest distance among the top 1% brightest 
% pixels as an adaptive threshold in the gray-scaled image since it varies
% from image to image. The top 1% brightest pixel was extracted by sorting 
% the value of the gray scaled images in descending order.
% After I retrieved data of the brightest pixels, I have to distinguish 
% the exudates to from the optic disc. Since it is quite difficult to 
% distinguish them by color, I decided to observe of the distribution of 
% the pixels. Brightest pixels of a healthy retina tends to cluster in one,
% while the ones with exudates usually spread across the image. 

clear;
folder = './distributed';

% prevent file error
if ~isfolder(folder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

filePattern = fullfile(folder, '*.ppm');
imgs = dir(filePattern);
Ms = zeros(length(imgs),1);
disp('code is running');
% extracting the feature and calculate scores
for k = 1:length(imgs)
    baseFileName = imgs(k).name;
    fullFileName = fullfile(folder, baseFileName);
    oimg = imread(fullFileName);
    img = double(oimg);
    
    r = mean(img(:,:,1));
    g = mean(img(:,:,2));
    b = mean(img(:,:,3));
    r = r*2.5;
    img = r/(r+g+b)*img(:,:,1)+g/(r+g+b)*img(:,:,2)+b/(r+g+b)*img(:,:,3);
    
    % sorting values of images in descending order
    V = sort(img(:), 'descend');
    top75 = V(ceil(end/10)*7.5);
    imin = top75;
    
    % top 1%
    top1 = V(ceil(end/10*.1));
    img(img<top1)=0;
    [x,y] = find(img ~= 0);
    x = [x,y];
    
    % prevent crash
    if size(x,1)>5000
      errorMessage = sprintf('Insufficient memory');
      uiwait(warndlg(errorMessage));
      return;
    end
    
    D = pdist(x); % distance of all pairs
    Ms(k) = max(D); % max distance of each image
    imagesc(oimg)
    colormap(gray)
end

healthy = Ms<size(img,2)*.4; % threshold changes as image size changes
scatter(1:length(Ms),Ms); % bottom left and top right is correct
hold on;
line([18.5 18.5], get(gca, 'ylim')); 
line(get(gca, 'xlim'), [size(img,2)*.4 size(img,2)*.4]);
score = sum(healthy(1:18)==1)+sum(healthy(19:36)==0); % correct cases
score