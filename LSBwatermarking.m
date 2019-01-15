clc;
clear all;
close all;

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