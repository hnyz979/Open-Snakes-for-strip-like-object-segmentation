function [Eextern, I_deri] = ExternalForceImage2D_1(I,Wline, Wedge, Wterm,Sigma,Verbose)
% 从图像计算外力
% 
% 输入： 
%  I        原图
%  Sigma	计算图像导数的因子 
%  Wline    文章中4式的三个权重值
%  Wedge	文章中4式的三个权重值
%  Wterm	文章中4式的三个权重值
% 输出：
%  Eextern	外能量函数，文章中的Gext

% 1 求出来图像对x或者y的各阶导数。
Ix=ImageDerivatives2D(I,Sigma,'x');
Iy=ImageDerivatives2D(I,Sigma,'y');
Ixx=ImageDerivatives2D(I,Sigma,'xx');
Ixy=ImageDerivatives2D(I,Sigma,'xy');
Iyy=ImageDerivatives2D(I,Sigma,'yy');
threshhold=1e-6;
Ix(abs(Ix)<threshhold)=0;
% 即Ix(find(Ix<threshhold))=0;把那些小于1e-6的都弄成0，然后通过那个Eterm(isnan(Eterm))=0把Eterm矩阵中相应元素弄成0。
% 否则的话，Eterm里会存在非常大的数值，如1e+11，对计算非常不利。
% 之所以取1e-6，是因为画出了(Ix.^2 + Iy.^2).^(3/2)的图，发现有意义的数值都在1e-5以上，而比他小的都是噪声。
Iy(abs(Iy)<threshhold)=0;
Ixx(abs(Ixx)<threshhold)=0;
Ixy(abs(Ixy)<threshhold)=0;
Iyy(abs(Iyy)<threshhold)=0;
% 2 计算Eline（4式中的G_l），其实就是原图。这个好像没啥用。
Eline = I;
if(Verbose),  figure,imshow(Eline,[]),title('Eline'); end%这两行画Eline的二维分布图，调试后可以注释掉。
% 3 计算Eterm（4式中的G_t）
Eterm = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((Ix.^2 + Iy.^2).^(3/2));
Eterm(isnan(Eterm))=0;
Eterm_threshhold=1e+0;
Eterm(abs(Eterm)>Eterm_threshhold)=0;
%Eterm真不能太大，应该把那个Eterm_threshhold顶多设成10（甚至设成1），防止噪声造成分母过小，对蛇的运动产生影响。
if(Verbose), figure,imshow(Eterm,[]),title('Eterm'); end%这两行画Eterm的二维分布图，调试后可以注释掉。
% 4 计算Eedge（4式中的G_e）
[grady,gradx] = gradient(I);%用另一种方法求Ix和Iy，看看效果如何。
Eedge = sqrt(gradx.^2 + grady.^2); 
if(Verbose), figure,imshow(Eedge,[]),title('Eedge'); end%%这两行画Eedge的二维分布图，调试后可以注释掉。
I_deri=Ix;
% 5 加和得到Eextern（文章中的Gext） = (Wline*Eline - Wedge*Eedge -Wterm * Eterm); 
Eextern= (Wline*Eline - Wedge*Eedge + Wterm * Eterm); 