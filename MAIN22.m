%% 直方图移位嵌入实验仿真
clc;clear all;close all;clear;

% 本文加解密结果
im = double(imread('date/woman.tif'));
cover = double(imread('date/bridge.bmp'));
 
% im = double(imread('date/lena_gray_512.tif'));
% cover = double(rgb2gray(imread('date/4.2.06.tiff')));
%  
% im = double(imread('date/4.2.01.tiff'));
% cover = double(imread('date/4.2.03.tiff'));
% 
% im = double(rgb2gray(imread('date/4.2.07.tiff')));
% cover = double(imread('date/house.tiff'));
% 
% im = double(imread('date/ZeldaColor.bmp'));
% cover = double(imread('date/5.3.01.tiff'));

key = '8d5ab8ba5340fce4420829ad5d12a0e45dacb0858544163d04c1d02b73e3697d';
cr=0.5;

emb_type = 1;  % 空域
tt=11;

% emb_type = 2;  %小波
% tt=11;

[M,N,K] = size(cover);
[m,n,k] = size(im);
im_shape = [m,n,k];


[ cip,dnkey,mmr,TT,tcip,rC ] = encryption( im,cover,key,cr,emb_type,tt );
[ rim,rtcip ] = dencryption( cip,cover,dnkey,mmr,TT,cr,emb_type,im_shape,rC );


[mssim_cip] = mssim(cover,cip,M-8);
[ssim_cip, ~] = ssim(uint8(cover),uint8(cip));
[psnr_cip, ~] = psnr(cover,double(uint8(cip)),255);

[mssim_rim] = mssim(im,rim,m-8);
[ssim_rim, ~] = ssim(uint8(im),uint8(rim));
[psnr_rim, ~] = psnr(im,double(uint8(rim)),255);
 

figure(1)
imshow(uint8(im))
figure(2)
imshow(uint8(cover))
figure(3)
imshow(uint8(tcip))
figure(4)
imshow(uint8(cip))
figure(5)
imshow(uint8(rim))
figure(6)
imshow(uint8(30*(cip-cover)))
figure(7)
imshow(uint8(30*(rim-im)))