function [IF, IF_y] = Fuse_color(Ie1,Ie2,D1,D2,A1,A2,I1ycbcr)
% fusion of multimodal images when one of them is in color and the other is in
% greyscale (anatomical-functional image fusion)


% applying max-abs fusion rule
A1(abs(A1)<abs(A2)) = 0; 
A2(abs(A2)<=abs(A1)) = 0;

% reconstructing the fused correlated-component
IzF = mexCombinePatches([D1 D2]*[A1;A2],zeros(size(Ie1)),sqrt(size(D1,1)),0,1);

IF_y = IzF + Ie1 + Ie2; % the luminance layer of the fused image

% clamping the fused image (standardization)
IF_y(IF_y>1)=1;
IF_y(IF_y<0)=0;

IF_ycbcr(:,:,2:3) = I1ycbcr(:,:,2:3); % the color layers from the functional image
IF_ycbcr(:,:,1) = IF_y; % the luminace layer is the fused image
IF = ycbcr2rgb(IF_ycbcr); % converting to rgb


