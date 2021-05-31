%%% greyscale-greyscale mutimodal image fusion (Infrared - Visible)

clear
% clc
addpath('utilities');

%% parameters

opts.T = 3; % number of nonzer entries in sparse vectors
opts.rho = 5; % penalty parameter
opts.print_res = true; % print decomposition results

%% loading input images

I1 = imresize(double(imread('Infrared_Visible_Images\IR2.png'))/255,1);
if size(I1,3)>1, I1 = rgb2gray(I1); end
I2 = imresize(double(imread('Infrared_Visible_Images\VIS2.png'))/255,1);
if size(I2,3)>1, I2 = rgb2gray(I2); end

%% performing decomposition and fusion
n = 16; % number of atoms in dictionaries
b = 8; % patch size
D0 = DCT(n,b);  % initializing the dictionaries with DCT matrices

tic;
[~,~,Ie1,Ie2,D1,D2,A1,A2] = perform_Corr_Ind_Decomp(I1,I2,D0,D0,opts); % Deocomposition
IF = Fuse_grey(Ie1,Ie2,D1,D2,A1,A2); % Fusion
toc % runtime

%% results

F = uint8(IF*255);
imwrite(F,'Results\IV_res.png');

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
