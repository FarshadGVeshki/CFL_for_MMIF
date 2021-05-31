function IF = Fuse_grey(Ie1,Ie2,D1,D2,A1,A2)
% fusion of multimodal images when when both inputs are in
% greyscale (anatomical-anatomical image fusion)


% applying max-abs fusion rule
A1(abs(A1)<abs(A2)) = 0; 
A2(abs(A2)<=abs(A1)) = 0;

% reconstructing the fused correlated-component
IzF = mexCombinePatches([D1 D2]*[A1;A2],zeros(size(Ie1)),sqrt(size(D1,1)),0,1); 

IF = IzF + Ie1 + Ie2;

% clamping the fused image (standardization)
IF(IF>1)=1;
IF(IF<0)=0;