clc;
clear;

i=imread('coke.jpg');
j=imresize(i,[1000, 1000]); % resizing taken image  
k=rgb2gray(j); % converting rgb image to gray image


figure
subplot(1,2,1)
imshow(k); % displaying objective image
title('Objective image');

x=imread('pepsi.jpg'); % image to be hidden


y=imresize(x,[1000, 1000]); % resizing hidden image
z=im2bw(y); % converting rbg to binary 


subplot(1,2,2);
imshow(z) % displaying image to be hidden
title('image to be hidden');

z=double(z); % increasing range to double

r=double(k-mod(k,2)); % removal of LSB bits 
l=uint8(r+z); % adding LSB bit from image to be hidden



figure
imshow(l)
title('Invisble watermarked Image'); 

%creating logistic maps m1 and m2
m1 = zeros(1,256);
m2 = zeros(1,256);
r    = 3.99; 
 m1(1) = 0.11;
 m2(1) = 0.23;
 N    = 256;
 for ii = 1:N-1
    m1(ii+1) = r*m1(ii)*(1 - m1(ii));
     m2(ii+1) = r*m2(ii)*(1 - m2(ii));
 end
 
 m1 = m1(129:256);
 m2 = m2(129:256);
 
%creating circulant matrices c1 and c2 from the logistic maps 
%function c = circulantMatrix(m)
c1 = zeros(128,128);
c1(1,:) = m1;
c2 = zeros(128,128);
c2(1,:) = m2;

for i = 2:128
   c1(i,:) = circshift(c1(i-1,:),1); 
end

for i = 2:128
   c2(i,:) = circshift(c2(i-1,:),1); 
end
   
%reducing the relevance among columns of circulant matrices
M = 128;
N = 128;
lambda = 2;

%c1(1,1) = lambda * c1(M,N);

for i = 2:M
    c1(i,1) = lambda * c1(i-1,N);
    c2(i,1) = lambda * c2(i-1,N);
end

for j = 2:N
    for i = 2:M
        c1(i,j) = c1(i-1,j-1);
        c2(i,j) = c2(i-1,j-1);
    end
end



%Read an image
I = l;
I = imresize(I, 0.5);
Idouble = im2double(I);

I1=I(1:size(I,1)/2,1:size(I,2)/2,:);
I2=I(size(I,1)/2+1:size(I,1),1:size(I,2)/2,:);
I3=I(1:size(I,1)/2,size(I,2)/2+1:size(I,2),:);
I4=I(size(I,1)/2+1:size(I,1),size(I,2)/2+1:size(I,2),:);

%I1 = rgb2gray(I1);
%I2 = rgb2gray(I2);
%I3 = rgb2gray(I3);
%I4 = rgb2gray(I4);
%I = rgb2gray(I);

%first block
%dct1 = c1 * I1 * transpose(c1);
d1 = dct2(I1);
d2 = dct2(I2);
d3 = dct2(I3);
d4 = dct2(I4);


d1(abs(d1) < 10) = 0;
d2(abs(d2) < 10) = 0;
d3(abs(d3) < 10) = 0;
d4(abs(d4) < 10) = 0;


L = ([d1 d3; d2 d4]);
L1 = ([I1 I3]);
L2 = ([I2 I4]);
  im = L;
% Get the dimensions of the image.  numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(im);
if numberOfColorBands > 1
	imgray = rgb2gray(im); % Convert to gray level.
end


% Display the original gray scale image.
subplot(1, 7, 2);
imshow(I);
title('Original Image');
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf, 'Name', 'image displays', 'NumberTitle', 'Off')
%%%%%%%%%%%%%%
compressedImage = im;

subplot(1, 7, 3);
imshow(log(abs(compressedImage)),[])
title('Compressed Image');

%SCRAMBLING BEGINS

% Get the order to scramble them in 
scrambleOrder = randperm(rows*columns);
% Scramble according to the scrambling order.
im = im(scrambleOrder);

% Reshape into a 2D image
scrambledImage = reshape(im, [rows, columns]);

% Display the scrambled gray scale image.
subplot(1, 7, 4);
imshow(log(abs(scrambledImage)),[]);
title('Scrambled Image');

%To reconstruct the original image

reconstruct = zeros(rows*columns, 2);
reconstruct(:, 1) = 1 : (rows*columns);
reconstruct(:, 2) = scrambleOrder;


% Sort this to find out where each scrambled location needs to be sent to.
newOrder = sortrows(reconstruct, 2);

% Extract just column 1, which is the order we need.
newOrder = newOrder(:,1);
% Unscramble according to the reconstruct order.
im = im(newOrder);
% Reshape into a 2D image
unscrambledImage = reshape(im, [rows, columns]);

% Display the original gray scale image.
subplot(1, 7, 5);
imshow(log(abs(unscrambledImage)), []);
title('Unscrambled Image');


I5=unscrambledImage(1:size(unscrambledImage,1)/2,1:size(unscrambledImage,2)/2,:);
I6=unscrambledImage(size(unscrambledImage,1)/2+1:size(unscrambledImage,1),1:size(unscrambledImage,2)/2,:);
I7=unscrambledImage(1:size(unscrambledImage,1)/2,size(unscrambledImage,2)/2+1:size(unscrambledImage,2),:);
I8=unscrambledImage(size(unscrambledImage,1)/2+1:size(unscrambledImage,1),size(unscrambledImage,2)/2+1:size(unscrambledImage,2),:);


%first block
%idct1 = transpose(c1) * I1 * c1;
K1 = idct2(d1);
K2 = idct2(d2);
K3 = idct2(d3);
K4 = idct2(d4);


L22 = ([K1 K3; K2 K4]);


% Display the original image.
subplot(1, 7, 6);
imshow(L22, []);
title('Uncompressed Image');

%detection of hidden image

h=mod(l,2);
p=zeros(1000,1000);

for x=1:1000
    for y=1:1000
        if(h(x,y)==1)
            p(x,y)=255;
        end
    end
end

s=im2bw(p);
figure
imshow(s); % hidden image
title('Recovered hidden image')

whos
