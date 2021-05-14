%%% color-greyscale mutimodal image fusion (functional-anatomical)

clear
% clc

addpath('utilities');

%% fusion problem

fusion_mods = 'T2-PET';
% fusion_mods = 'T2-TC';
% fusion_mods = 'T2-TI';
% fusion_mods = 'Gad-PET';

%% parameters

opts.k = 5; % maximum nnonzero entries in sparse vectors
opts.rho = 10; % optimization penalty term
opts.plot = false; % plot decomposition components

%% loading input images

I1rgb = double(imread(['Source_Images\' fusion_mods '_A.png']))/255;
I1ycbcr = rgb2ycbcr(I1rgb);
I1 = I1ycbcr(:,:,1);
I2 = double(imread(['Source_Images\' fusion_mods '_B.png']))/255;
if size(I2,3)>1, I2 = rgb2gray(I2); end

%% performing decomposition and fusion

n = 32; b = 8;
D0 = DCT(n,b);  % initializing the dictionaries with DCT matrices

tic;
[~,~,Ie1,Ie2,D1,D2,A1,A2] = perform_Corr_Ind_Decomp(I1,I2,D0,D0,opts); % Decomposition
[IF, IF_int] = Fuse_color(Ie2,Ie1,D2,D1,A2,A1,I1ycbcr); % Fusion
toc; % runtime

%% results
F = uint8(IF*255);
imwrite(F,['Results\' fusion_mods '_F.png']);

figure(23)
subplot 131
imshow(I1rgb,[])
xlabel('I_1')
subplot 132
imshow(I2,[])
xlabel('I_2')
subplot 133
imshow(IF,[])
xlabel('I^F')

%%  dictionary atoms
% ID1 = displayPatches(D1);
% ID2 = displayPatches(D2);
% 
% figure(37)
% subplot 121
% imshow(ID1)
% xlabel('D1')
% subplot 122
% imshow(ID2)
% xlabel('D2')




