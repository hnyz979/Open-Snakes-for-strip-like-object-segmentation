function J=ImageDerivatives2D(I,sigma,type)
% 基于高斯分布，求二维图像对于x或y方向的一阶、二阶导数。
% 输入：
%   I       输入图像（2D）
%   sigma	计算图像导数的因子
%   type	五种，分别为'x', 'y', 'xx', 'xy', 'yy'，即计算导数的方向
% 输出：
%   J       图像导数

% 1 生成对图像求导用到的高斯卷积核。
[x,y]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
% 生成的卷积核DGauss是一个矩阵，类似于各种边缘提取算子的作用。
switch(type)
    case 'x'
        DGauss=-(x./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
    case 'y'
        DGauss=-(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
    case 'xx'
        DGauss = 1/(2*pi*sigma^4) * (x.^2/sigma^2 - 1) .* exp(-(x.^2 + y.^2)/(2*sigma^2));
    case {'xy','yx'}
        DGauss = 1/(2*pi*sigma^6) * (x .* y)           .* exp(-(x.^2 + y.^2)/(2*sigma^2));
    case 'yy'
        DGauss = 1/(2*pi*sigma^4) * (y.^2/sigma^2 - 1) .* exp(-(x.^2 + y.^2)/(2*sigma^2));
end
% 2 通过用原图与高斯核卷积（imfilter函数）的方法，求出图像灰度的导数。
J = imfilter(I,DGauss,'conv','symmetric');