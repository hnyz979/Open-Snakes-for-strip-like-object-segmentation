clear
close all  
clc
% -----------------------0、读入首帧的手动分割结果。-----------------------
root_dir = ['B:\02-动脉血管追踪数据（给师弟的补充文章）\matlab程序\上传github的版本（基于备份25）\' ...
    'Open-Snakes-for-strip-like-object-segmentation\data\'];
program_dir = ['B:\02-动脉血管追踪数据（给师弟的补充文章）\matlab程序\上传github的版本（基于备份25）\' ...
    'Open-Snakes-for-strip-like-object-segmentation\sequence_segmentation_using_manual_initial'];
% 以上分别是存放数据的目录和存放程序的目录。
sequence_name = '1';
load([root_dir sequence_name '\LI.mat']);  % load出来之后，LI的初始值就是red_manu
load([root_dir sequence_name '\MA.mat']);  % load出来之后，MA的初始值就是blue_manu
I=imread([root_dir sequence_name '\000001'],'bmp');
rmax_index=fix(red_manu);%变成整形。。
bmax_index=fix(blue_manu);
xlen=length(red_manu);
if(size(I,3)==3), I=rgb2gray(I); end%把I变成灰度图。
figure, imshow(I), hold on;
plot(1:xlen,rmax_index,'r.',1:xlen,bmax_index,'b.');
title('第1幅图')
hold off
Red=[1:xlen; rmax_index']';%LI初始控制点集合Red，现在是个2*xlen的矩阵。
Blue=[1:xlen; bmax_index']';%MA初始控制点集合Blue，现在是个2*xlen的矩阵。
% -----------------------1、判断首帧灰度（噪声）、斜率等情况，作为后续帧处理的基准。-----------------------
grayscale_rmax_up_mean=zeros(1,xlen);
grayscale_rmax_down_mean=zeros(1,xlen);
grayscale_bmax_up_mean=zeros(1,xlen);
grayscale_bmax_down_mean=zeros(1,xlen);grayscale_bmax_mid_std=zeros(1,xlen);
inside_IMC_is_bright=zeros(1,xlen);%看看是否有MA上方比下方还亮的，即，可能不太正常的。
inside_IMC_is_dark=zeros(1,xlen);%看看是否有LI下方很暗的，即，可能不太正常的。
IMW=bmax_index-rmax_index;%初始轮廓的IM厚度。是个1*361的向量，用于设置邻域和贴近修正阈值。
% 1.1 首先通过灰度值，检测噪声情况。
% 1.1.1 观察LI上方或MA下方空间内的灰度值，通过噪声情况，修正外邻域（LI上方或者MA下方）的大小。
neighbor_out_r=min(max(4,fix(IMW*0.7)),9);%根据IMW值选取外邻域大小初始值。
neighbor_out_b=min(max(4,fix(IMW*0.7)),9);
for i=1:xlen
    r_small=I(rmax_index(i)-1:rmax_index(i)+1,max(1,i-2):min(i+2,xlen));
    b_small=I(bmax_index(i)-1:bmax_index(i)+1,max(1,i-2):min(i+2,xlen));
    for j=rmax_index(i)-4:-1:rmax_index(i)-9%在初始LI上方搜索，从上方4个点搜索到9个点。
        if i==104
            pause(0.1)  % 这个可以用来逐个检查每一个点的情况。
        end
        if I(j,i)>=mean2(r_small)+15%如果有一个点的灰度，比起LI处的平均灰度大了15，就认为这一点很可能是噪声了。
            neighbor_out_r(i)=rmax_index(i)-j-1;%那么，记录下这个点的位置。对于这个i，LI的上邻域就取rmax_index(i)-j。
            break
        end
    end
    for k=bmax_index(i)+4:bmax_index(i)+9%初始MA下方搜索，从下方4个点搜索到下方9个点。
        if I(k,i)<=mean2(b_small)-15%如果有一个点的灰度，比起MA处的平均灰度小了15，就认为这一点很可能是不对了。比如，MA下方的、上白下黑的地方。
            neighbor_out_b(i)=k-bmax_index(i)-1;
            break
        end
    end
end
% 1.1.2 确定外邻域（LI上方和MA下方）和内邻域（LI下方和MA上方）的平均灰度（噪声水平）。
neighbor_in=min(max(4,fix(IMW/3)),6);%内邻域，如果IM很窄的话，就仍然取4，否则就加宽一些
for i=1:xlen
    r_up=I(rmax_index(i)-neighbor_out_r(i)-1:rmax_index(i)-1,max(1,i-2):min(i+2,xlen));
    r_down=I(rmax_index(i)+1:rmax_index(i)+neighbor_in(i)+1,max(1,i-2):min(i+2,xlen));%类似上面，只不过这儿是加了1
    b_up=I(bmax_index(i)-neighbor_in(i)-1:bmax_index(i)-1,max(1,i-2):min(i+2,xlen));%同r_up
    b_down=I(bmax_index(i)+1:bmax_index(i)+neighbor_out_b(i)+1,max(1,i-2):min(i+2,xlen));%同r_down
    b_mid=I(bmax_index(i)-neighbor_in(i):bmax_index(i)+neighbor_in(i),max(1,i-2):min(i+2,xlen));
    grayscale_rmax_up_mean(i)=mean2(r_up);%求初始LI上方邻域内的平均灰度。
    % 例如，初始LI上一点为(3,202)，则就求图像中从(1,198)开始到(5,202)结束的这25个点的平均灰度。
    grayscale_rmax_down_mean(i)=mean2(r_down);%求初始LI下方邻域内的平均灰度。
    % 例如，初始LI上一点为(3,202)，则就求图像中从(1,202)开始到(5,206)结束的这25个点的平均灰度。
    grayscale_bmax_up_mean(i)=mean2(b_up);%求初始MA上方邻域内的平均灰度。
    grayscale_bmax_down_mean(i)=mean2(b_down);%求初始MA下方邻域内的平均灰度。
    grayscale_bmax_mid_std(i)=std2(b_mid);%求初始MA下方邻域内的灰度标准差。
    grayscale_inside=I(rmax_index(i)+1:bmax_index(i),i);
    if sum(grayscale_inside>grayscale_bmax_down_mean(i))>=10
        % 如果这一列，内中膜之间有10个以上的点的灰度大于下方平均灰度，就做个标记，不要轻易加能量值，使得MA往上移动推。
        % 此处对应文章P4的两个exceptions，有改动。
        inside_IMC_is_bright(i)=1;
    end
    if sum(grayscale_inside<=10)>=10
        % 如果这一列，内中膜之间有10个以上的点的灰度小于10，就做个标记，不要轻易加能量值，使得LI往下移动。
        inside_IMC_is_dark(i)=1;
    end
end
exists_noise_temp=find(grayscale_rmax_up_mean>=10);%LI上方较大区域灰度值大于10，初步认为有噪声。
temp=0;
for k=1:size(exists_noise_temp,2)
    % 用这个循环，找出exists_noise_temp中的点的10-邻域内所有点。例如，如果exists_noise_temp中有一点67，那么从57~77都认为可疑有噪声。
    qu=ismember(max(1,exists_noise_temp(k)-2):min(exists_noise_temp(k)+2,xlen),exists_noise_temp);
    % 判断第k个数左右两边的数是否也在exists_noise_temp里，也就是说，判断是否连着5个点都有噪声。
    if min(abs(qu))==1
        % 如果qu里的5个数都是1的话，就认为这一带噪声大，才找出10邻域内的点。
        sus=max(1,exists_noise_temp(k)-10):min(exists_noise_temp(k)+10,xlen);
    else
        sus=exists_noise_temp(k);
    end
    temp=[temp sus];
end
temp(temp==0)=[];%删掉第一个0
exists_noise=unique(temp);
medium_noise=find(grayscale_rmax_up_mean>=20 & grayscale_rmax_up_mean<30);%LI上方灰度值大于20且小于30，认为噪声中等严重。
severe_noise=find(grayscale_rmax_up_mean>=30);%LI上方灰度值大于30，认为噪声比较严重。
%根据噪声严重程度，决定蛇函数中灰度达到多少才加Grayscale constraint energy。
exists_fuzzy_temp=find(grayscale_bmax_mid_std<=20);%
temp=0;
for k=1:size(exists_fuzzy_temp,2)
    %类似于前面，处理MA附近邻域。
    qu=ismember(max(1,exists_fuzzy_temp(k)-2):min(exists_fuzzy_temp(k)+2,xlen),exists_fuzzy_temp);
    if min(abs(qu))==1
        sus=max(1,exists_fuzzy_temp(k)-10):min(exists_fuzzy_temp(k)+10,xlen);
    else
        sus=exists_fuzzy_temp(k);
    end
    temp=[temp sus];
end
temp(temp==0)=[];%删掉第一个0
exists_fuzzy=unique(temp);
medium_fuzzy=find(grayscale_bmax_mid_std<=15 & grayscale_bmax_mid_std>10);
severe_fuzzy=find(grayscale_bmax_mid_std<=10);
lightflag=0;lightflag1=0;noise_flag0=0;noise_flag1=0;noise_flag2=0;fuzzy_flag1=0;
if size(exists_noise_temp,2)>15%用的是exists_noise_temp，即没扩展之前的
    noise_flag0=1;
    disp('本序列LI上方有噪声')
end
if size(find(grayscale_bmax_up_mean>=150),2)>200%如果有200个以上的点，他的MA上方灰度值大于150。对1055L1用的。
    lightflag=1;%snake2D函数中，MA向上鼓包往下推的时候会有一句，“if lightflag==1, vol_times=vol_times*4;end”
    disp('本序列MA上方很亮')%已测试，只有1055L1MA上方很亮
end
if size(find(grayscale_bmax_down_mean>=200),2)>200%如果有200个以上的点，他的MA下方灰度值大于200。
    lightflag1=1;%
    disp('本序列MA下方很亮')%已测试，39R，47L/R，48R，49L/R，55的四个，59R，60L/R，62L/R，63L/R，64L，66L1满足此条件。
end
if size(exists_noise_temp,2)>200 && size(medium_noise,2)>=30%50R
    noise_flag1=1;
    disp('本序列LI上方噪声严重')
end
if size(severe_noise,2)>size(medium_noise,2) && size(severe_noise,2)>=100%
    noise_flag2=1;
    disp('本序列LI上方噪声严重，且严重噪声比一般噪声还多')
end
if size(exists_fuzzy_temp,2)>=30%用的是exists_fuzzy_temp，即没扩展之前的
    fuzzy_flag1=1;
    disp('本序列MA有30个点以上模糊')
end
% 1.2 检查初始轮廓中的导数情况。对应文中Derivative constraint energy。
red_dot1=Red(:,2);blue_dot1=Blue(:,2);%初始化，其实都弄成0也可以的。
for i=1:xlen-1, red_dot1(i)=Red(i+1,2)-Red(i,2);blue_dot1(i)=Blue(i+1,2)-Blue(i,2);end%计算导数
red_dot1(xlen)=red_dot1(xlen-1);blue_dot1(xlen)=blue_dot1(xlen-1);
% -----------------------2、读入第2帧图，用蛇算法演化后，设置文章中H-inf的参数。-----------------------
% 2.1 演化得到第2帧图的LI和MA。
I2=imread([root_dir sequence_name '\000002'],'bmp');
if(size(I2,3)==3), I2=rgb2gray(I2); end
sigma_filt=1;
I2 = gaussian_filter(I2,sigma_filt);
I2(:,1)=I2(:,2);I2(:,xlen)=I2(:,xlen-1);%滤波的时候，第一列和最后一列会出问题，特别黑。现在改掉。
Options=struct;
Options.Iterations=20;
Options.Verbose=false;
i=1;
fast_r=0;fast_b=0;%LI和MA是否运动非常快的flag。先初始化为0，在某个序列中，一旦这两个东西被置1，就不会再返回为0了。
[Red1,Blue1,J,fast_r,fast_b]=Snake2D(i,I2,Red,Blue,Options,neighbor_in,neighbor_out_r,neighbor_out_b, ...
    grayscale_rmax_up_mean,grayscale_rmax_down_mean,grayscale_bmax_up_mean,grayscale_bmax_down_mean, ...
    inside_IMC_is_bright,inside_IMC_is_dark,severe_noise,medium_noise,exists_noise, ...
    severe_fuzzy,medium_fuzzy,exists_fuzzy,IMW,red_dot1,blue_dot1,fast_r,fast_b,...
    lightflag,lightflag1,noise_flag0,noise_flag1,noise_flag2,fuzzy_flag1);%用蛇确定下一幅图的轮廓
% 2.2 求解状态空间方程。
Options.Verbose=false;
q=1;
r=1;%测量噪声协方差R越小（趋于零）， 测量变量（蛇的结果）的权重越来越大
matrixSize=size(Red,1);
A = eye(matrixSize,matrixSize);
H = eye(matrixSize,matrixSize);
Q = q*eye(matrixSize,matrixSize);
R = r*eye(matrixSize,matrixSize);
% 用karmen及H inf的时候，按理说，只取Red(:,2)和Blue(:,2)就可以了，因为他们的第一维是不变的，就是1：len。
Red_new=Red; Blue_new=Blue;
u=0;G=0;
L=eye(matrixSize,matrixSize);S=eye(matrixSize,matrixSize);I=eye(matrixSize,matrixSize);
theta=0.1;
[Red_new(:,2),K_red] = Hinf(Red(:,2),Red1(:,2),u,A,G,H,L,R,Q,S,I,theta);
[Blue_new(:,2),K_blue] = Hinf(Blue(:,2),Blue1(:,2),u,A,G,H,L,R,Q,S,I,theta);
%实验表明，如果θ=0.5，那么，K_red和K_blue都是单位阵，这意味着下一次的LI和MA将只和蛇的结果有关，而和上次的结果无关。见02那本书的11.102式。
figure, imshow(I2,[]), hold on; 
plot(Red(:,1),Red(:,2),'m.');
plot(Red_new(:,1),Red_new(:,2),'r.');
plot(Blue(:,1),Blue(:,2),'c.');
plot(Blue_new(:,1),Blue_new(:,2),'b.');hold off;
title('第2幅图')
% -----------------------3、逐一读入后续图像。现在想想，第2帧和后续帧可以放在一起做。先不花时间修改了。-----------------------
srcDir=[root_dir sequence_name]; %选择那些图像的目录
cd(srcDir);%进入该目录
allnames=struct2cell(dir('*.bmp')); %该目录下的所有bmp图像
[~,len]=size(allnames);%len是图像张数
str1='第';str2='幅图';%画图时写题目用的
%初始化一个3维的Red_3D，第一维和第二维与前面相同，第三维则表示它是第k张图的结果。Blue与之相同。
Red_3D=zeros(size(Red,1),size(Red,2),len);
Blue_3D=zeros(size(Blue,1),size(Blue,2),len);
Red_3D(:,:,1)=[Red(:,1),Red(:,2)];%000001图中，原始手动画的结果
Red_3D(:,:,2)=[Red_new(:,1),Red_new(:,2)];%000002图中，蛇+卡尔曼的结果。
Blue_3D(:,:,1)=[Blue(:,1),Blue(:,2)];
Blue_3D(:,:,2)=[Blue_new(:,1),Blue_new(:,2)];
red_dot2=Red_new(:,2);blue_dot2=Blue_new(:,2);%初始化，其实都弄成0也可以的。
for i=1:xlen-1, red_dot2(i)=Red_new(i+1,2)-Red_new(i,2);blue_dot2(i)=Blue_new(i+1,2)-Blue_new(i,2);end%计算导数
red_dot2(xlen)=red_dot2(xlen-1);blue_dot2(xlen)=blue_dot2(xlen-1);
%盖楼，每一层都是一个采样时刻，即一个k。
%从第3张开始，这样的话，Red_3D(:,:,i)就对应于第i幅图的LI的横纵坐标。
figure,
for i=3:min(100,len)
    name=allnames{1,i};
    I=imread(name);%此时I是第i张图读入进来的矩阵。
    if(size(I,3)==3), I=rgb2gray(I); end
    cd(program_dir)%刚刚是存放图片的文件夹，现在要回到存放程序的文件夹。
    I = gaussian_filter(I,sigma_filt);
    I(:,1)=I(:,2);I(:,xlen)=I(:,xlen-1);%滤波的时候，第一列和最后一列会出问题，特别黑。现在改掉。
    if i==30 %想要检查那一张图，就让i等于几，然后加断点。断点后单步调试，进入Snake2D程序即可观察此时的Eext等矩阵。
        pause(0.1);
    end
    [Red_3D(:,:,i),Blue_3D(:,:,i),J,fast_r,fast_b]=Snake2D(i,I,Red_3D(:,:,i-1),Blue_3D(:,:,i-1),Options, ... 
        neighbor_in,neighbor_out_r,neighbor_out_b,grayscale_rmax_up_mean,grayscale_rmax_down_mean,grayscale_bmax_up_mean, ...
        grayscale_bmax_down_mean,inside_IMC_is_bright,inside_IMC_is_dark,severe_noise,medium_noise,exists_noise, ...
        severe_fuzzy,medium_fuzzy,exists_fuzzy,IMW,red_dot2,blue_dot2,fast_r,fast_b,lightflag,lightflag1,...
        noise_flag0,noise_flag1,noise_flag2,fuzzy_flag1);
    Red_3D(:,2,i) = A * Red_3D(:,2,i-1) + G * u + A * K_red * (Red_3D(:,2,i) - H * Red_3D(:,2,i-1));
    Blue_3D(:,2,i) = A * Blue_3D(:,2,i-1) + G * u + A * K_blue * (Blue_3D(:,2,i) - H * Blue_3D(:,2,i-1));
    %上面一行，等号右边Red_3D(:,2,i-1)是上一幅图的结果，Red_3D(:,2,i)是这一幅图的观测值（蛇的结果），K_red是11.89式中的K值，假设为不变的。
    %然后等号左边的Red_3D(:,2,i)是这一幅图的结果（蛇且H inf后的结果）。
    imshow(I,[])
    hold on; 
    plot(Red_3D(:,1,i),Red_3D(:,2,i),'r.');%Red_3D(:,1,i)仍然是1:len那些自然数，Red_3D(:,2,i)是第k幅图处理过的LI坐标。
    plot(Blue_3D(:,1,i),Blue_3D(:,2,i),'b.'); 
    str=[str1 num2str(i) str2];
    title(str)
    hold off
    pause(0.01);
    cd (srcDir)
end
% save B:\02-动脉血管追踪数据（给师弟的补充文章）\matlab程序\上传github的版本（基于备份25）\Red_1Dsnake_Hinf Red_3D
% save B:\02-动脉血管追踪数据（给师弟的补充文章）\matlab程序\上传github的版本（基于备份25）\Blue_1Dsnake_Hinf Blue_3D