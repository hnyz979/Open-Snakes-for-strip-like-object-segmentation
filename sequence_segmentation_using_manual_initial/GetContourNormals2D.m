function N=GetContourNormals2D(P)
% 利用轮廓上每个点的若干个相邻点，计算轮廓上该点的法线方向。
% 输入：
%  P    2*N的轮廓点，不必须首尾相接
% 输出：
%  N    2*N的轮廓法线方向

% 要用几个临近点？
a=4;
% 搞出来轮廓每个点的x和y坐标。
xt=P(:,1); yt=P(:,2);

% 轮廓的导数，即，点坐标错位相减。
n=length(xt);
f=(1:n)+a; f(f>n)=f(f>n)-n;
b=(1:n)-a; b(b<1)=b(b<1)+n;
dx=xt(f)-xt(b);
dy=yt(f)-yt(b);

% 归一化，让法线向量的长度为1，得到法向量。
l=sqrt(dx.^2+dy.^2);
nx = -dy./l; % 
ny =  dx./l;
N(:,1)=nx; N(:,2)=ny;