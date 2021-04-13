%%% greyscale-greyscale mutimodal image fusion (anatomical-anatomical)

clear
% clc
addpath('utilities');

%% fusion problem

fusion_mods = 'CT-T2';
% fusion_mods = 'T1-T2';

%% parameters

opts.k = 5;
opts.rho = 10; % penalty parameter
opts.plot = true; % plot decomposition parameters

%% loading input images

I1 = double(imread(['Source_Images\' fusion_mods '_A.png']))/255;
if size(I1,3)>1, I1 = rgb2gray(I1); end
I2 = double(imread(['Source_Images\' fusion_mods '_B.png']))/255;
if size(I2,3)>1, I2 = rgb2gray(I2); end

%% performing decomposition and fusion
n = 36; % number of atoms in dictionaries
b = 8; % patch size
D0 = DCT(n,b);  % initializing the dictionaries with DCT matrices

tic;
[~,~,Ie1,Ie2,D1,D2,A1,A2] = perform_Corr_Ind_Decomp(I1,I2,D0,D0,opts); % Deocomposition
IF = Fuse_grey(Ie1,Ie2,D1,D2,A1,A2); % Fusion
toc % runtime

%% results

F = uint8(IF*255);
imwrite(F,['Results\' fusion_mods '_F.png']);

figure(23)
subplot 131
imshow(I1,[])
xlabel('I_1')
subplot 132
imshow(I2,[])
xlabel('I_2')
subplot 133
imshow(IF,[])
xlabel('I^F')

%%  dictionary atoms
ID1 = displayPatches(D1);
ID2 = displayPatches(D2);

figure(37)
subplot 121
imshow(ID1)
xlabel('D1')
subplot 122
imshow(ID2)
xlabel('D2')