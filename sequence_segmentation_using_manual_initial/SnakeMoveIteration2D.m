function P=SnakeMoveIteration2D(B,P,Fext,gamma,kappa,delta,eliminate_x_force)
% 用外部力和内部力，计算轮廓(蛇)的运动。
% 输入：
%   B       内力矩阵
%   P       轮廓点，2*N的矩阵
%   Fext    外力矩阵（文章中的5式）
%   gamma	步长，即这些力作用在轮廓上的时间，时间越长，运动就越多
%   kappa	外部力的权重，越大，外力对轮廓的影响越大
%   delta	气球力的权重
% 输出：
%   P       输出轮廓点，2*N的矩阵

% 把超出边界的轮廓点，弄回边界处
P(:,1)=min(max(P(:,1),1),size(Fext,2));
P(:,2)=min(max(P(:,2),1),size(Fext,1));

% 1 轮廓点上的外力。
Fext1=zeros(size(P,1),1);
Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(:,1),P(:,2));%把原来计算出来的外力矩阵的y分量，插值到每个轮廓点上
if (eliminate_x_force),Fext1(:,2)=kappa*interp2(Fext(:,:,2),P(:,1),P(:,2));end
% 类似地，如果eliminate_x_force=1，即保留受力的x分量，就把它也插值到每个轮廓点上
% 但是，为了简化计算，我们只让轮廓在y方向运动，所以这个Fext1(:,2)（x方向受的力）保持为初始值0。
Fext1(isnan(Fext1))=0;  % 把nan都变成0。

% 2 轮廓点的气球力（就是轮廓法向量乘以一个因子delta）。
N=GetContourNormals2D(P);
Fext2=delta*N;
Fext2(isnan(Fext2))=0;

% 3 按照两个力，乘以作用时间步长，得到距离。
% 有个地方搞反了，P(:,1)是图像中的x坐标，但是，在这个程序中，Fext1(i,1)是第i点y方向受的力，
% 而Fext2(i,1)是第i点y方向的气球力，所以加的时候，会是扭了一下。。。
if (eliminate_x_force), ssx = P(:,1) + gamma*Fext1(:,2) + gamma*Fext2(:,2); end
ssy = P(:,2) + gamma*Fext1(:,1) + gamma*Fext2(:,1); % 这个是y方向受力。

% 4 处理边界条件
len=size(P,1);
expand_len=10;  % 边界要多管几个点？
% 先是x方向
if (eliminate_x_force)%eliminate_x_force=1表示保留x方向受力
    ssx_expanded(expand_len+1:len+expand_len)=ssx;%增广P向量，目测应该是加了图像力之后，而不是之前。
    ssx_expanded(expand_len)=ssx(1)-1;
%     for i=1:expand_len-1,ssx_expanded(i)=ssx(expand_len-i);end%镜面反射边界条件
    for i=1:expand_len-1,ssx_expanded(i)=ssx(expand_len);end%相等边界条件
    ssx_expanded(len+expand_len+1)=ssx(len)+1;
%     for i=len+expand_len+2:len+expand_len*2,ssx_expanded(i)=ssx(len-(i-(len+expand_len+1)));end%镜面反射边界条件
    for i=len+expand_len+2:len+expand_len*2,ssx_expanded(i)=ssx(len);end%相等边界条件
    P_expand(:,1)=B*ssx_expanded'; 
    P(:,1)= P_expand(expand_len+1:len+expand_len,1);
end
% 然后是y方向
ssy_expanded(expand_len+1:len+expand_len)=ssy;%处理边界条件（轮廓跑到图像边缘了怎么办）。
% ssy_expanded(expand_len)=ssy(1)-1;%镜面反射边界条件用，否则注释掉
% for i=1:expand_len-1,ssy_expanded(i)=ssy(expand_len-i);end%镜面反射边界条件用，否则注释掉
% ssy_expanded(len+expand_len+1)=ssy(len)+1;%镜面反射边界条件用，否则注释掉
% for i=len+expand_len+2:len+expand_len*2,ssy_expanded(i)=ssy(len-(i-(len+expand_len+1)));end%镜面反射边界条件用，否则注释掉
% for i=1:expand_len,ssy_expanded(i)=ssy(expand_len+1);end%相等边界条件用，否则注释掉
% for i=len+expand_len+1:len+expand_len*2,ssy_expanded(i)=ssy(len);end%相等边界条件用，否则注释掉
ssy_dot_start=zeros(3,1);%斜率相等边界条件用，否则注释掉。这个是初始的三个点的斜率，都置零。
for i=1:3, ssy_dot_start(i)=ssy(i+1)-ssy(i);end%斜率相等边界条件用，否则注释掉。这个是ssy中的前三个点（即边界）的斜率。
mean_ssy_dot_start=mean(ssy_dot_start);%斜率相等边界条件用，否则注释掉。前三个点斜率平均值。
for i=expand_len:-1:1, ssy_expanded(i)=ssy_expanded(i+1)-mean_ssy_dot_start/(expand_len-i+1);end%斜率相等边界条件用，否则注释掉。
% 扩展的边界中的力的大小，让他越远（索引i越小），力就越小。
ssy_dot_end=zeros(3,1);%斜率相等边界条件用，否则注释掉。类似于前面，现在处理后三个点的斜率。
for i=1:3, ssy_dot_end(i)=ssy(len-3+i)-ssy(len-4+i);end%斜率相等边界条件用，否则注释掉。
mean_ssy_dot_end=mean(ssy_dot_end);%斜率相等边界条件用，否则注释掉。
for i=len+expand_len+1:len+expand_len*2, ssy_expanded(i)=ssy_expanded(i-1)+mean_ssy_dot_end/(i-len-expand_len);end%斜率相等边界条件用，否则注释掉。
% 类似地，最后的边界，也是让他越远（索引i越大），力就越小。

% 5 考虑内力B的影响，计算出来新的轮廓位置。
P_expand(:,2)=B*ssy_expanded'; 
P(:,2)= P_expand(expand_len+1:len+expand_len,2);  % 把补的边界删除掉。
% 处理边界上的点（前5个和后5个），保持按照那个平均斜率走，不让他们“飘起来”
% 实际情况下，如果不加这个的话，图像边缘的一些很强的灰度梯度，会把蛇结果带走，特别烦。索性加了个土办法。
P_dot_start=zeros(5,1);
for i=1:5, P_dot_start(i)=P(i+5,2)-P(i+4,2); end%乘完B后，第6~10个点的导数
mean_P_dot_start=mean(P_dot_start);
for i=5:-1:1, P(i,2)=P(i+1,2)-mean_P_dot_start;end%
P_dot_end=zeros(5,1);
for i=1:5, P_dot_end(i)=P(len-10+i,2)-P(len-11+i,2); end%乘完B后，导数第10~6个点的导数
mean_P_dot_end=mean(P_dot_end);
for i=1:5, P(len-5+i,2)=P(len-6+i,2)+mean_P_dot_end;end%
% 如果超出边界了，就让他回到边界上。
P(:,1)=min(max(P(:,1),1),size(Fext,2));
P(:,2)=min(max(P(:,2),1),size(Fext,1));