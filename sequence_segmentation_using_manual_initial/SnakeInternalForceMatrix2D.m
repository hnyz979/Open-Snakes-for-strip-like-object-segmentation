function B=SnakeInternalForceMatrix2D(nPoints,alpha,beta)
% 计算蛇的内力矩阵（保持蛇的光滑性）
% 输入：
%   nPoints     蛇上有几个轮廓点？
%   alpha       文中4式的alpha（曲线一阶导数对应的力）
%   beta        文中4式的beta（曲线二阶导数对应的力）
%   gamma       时间步长
% 输出：
%   B           内力矩阵

% 文中5式中，后面5项的系数（分母都不要）。忽然发现，文中倒数第二项好像多了个2。
b(1)=beta;
b(2)=-(alpha + 4*beta);
b(3)=(2*alpha + 6 *beta);
b(4)=b(2);
b(5)=b(1);

% 构造M矩阵（见Everything You Always Wanted To Know About Snakes第20页6.6式）。那个文章，有几个地方也有问题。。。
nPoints=nPoints+20;%因为蛇不连续，所以给A和B矩阵增广，前后各增广10行10列。
% 注意，现在本函数返回的那个B矩阵，就是增广后的B矩阵。到时候在SnakeMoveIteration2D里再把那个P向量前后做个镜面反射的边界条件。
A=b(1)*circshift(eye(nPoints),2);%M矩阵中右上角的两个p，和最下面的那一斜条p
A=A+b(2)*circshift(eye(nPoints),1);%右上角的一个q，和倒数第二靠下的那一斜条q
A=A+b(3)*circshift(eye(nPoints),0);%中间对角线上的r，但是r应该是1+2*alpha + 6 *beta，现在没加1，一会儿再加。
A=A+b(4)*circshift(eye(nPoints),-1);
A=A+b(5)*circshift(eye(nPoints),-2);
A=A + eye(nPoints);
% 现在完全按照01中的来干，就把原来的B=inv(gamma * A + eye(nPoints));或者B=inv(A + gamma.* eye(nPoints));中的gamma去掉了。
B=inv(A);