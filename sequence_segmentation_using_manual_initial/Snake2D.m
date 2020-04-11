function [P,Q,J1,fast_r,fast_b]=Snake2D(fr,I,P,Q,Options,neighbor_in,neighbor_out_r,neighbor_out_b,...
    grayscale_rmax_up_mean,grayscale_rmax_down_mean,grayscale_bmax_up_mean,grayscale_bmax_down_mean,...
    inside_IMC_is_bright,inside_IMC_is_dark,severe_noise,medium_noise,exists_noise,...
    severe_fuzzy,medium_fuzzy,exists_fuzzy,IMW,red_dot_ref,blue_dot_ref,fast_r,fast_b,...
    lightflag,lightflag1,noise_flag0,noise_flag1,noise_flag2,fuzzy_flag1)
% ----------------------------------------------输入----------------------------------------------
%   I    当前帧图像矩阵。
%   P   上一时刻的后验估计（LI、红色线），作为这一时刻LI蛇的初值。
%   Q   上一时刻的后验估计（MA、蓝色线），作为这一时刻MA蛇的初值。
%   Options   滤波、控制选项。
%   neighbor_in/neighbor_out_r/neighbor_out_b   首帧中的内邻域、LI外邻域、MA外邻域。
%   grayscale_rmax_up_mean/grayscale_rmax_down_mean   初始LI上/下方邻域内的平均灰度。
%   grayscale_bmax_up_mean/grayscale_bmax_down_mean   初始MA上/下方邻域内的平均灰度。
%   inside_IMC_is_bright/inside_IMC_is_dark   初始LI/MA之间区域内，灰度是否过高或过低。
%   severe_noise/medium_noise/exists_noise    初始LI上方，是否有严重噪声、中等噪声、轻度噪声。
%                                             用于确定何时加入Grayscale constraint energy。
%   severe_fuzzy,medium_fuzzy,exists_fuzzy    初始MA下方，是否有严重噪声、中等噪声、轻度噪声。
%                                             用于确定何时加入Grayscale constraint energy。
%   IMW   初始的内中膜厚度（LI-MA间距）
%   red_dot_ref/blue_dot_ref   初始的LI和MA每个点的斜率，用于Derivative constraint energy。
%   fast_r/fast_b   是否发现LI/MA在某些帧中运动速度较快。
%   lightflag/lightflag1/noise_flag0/noise_flag1/noise_flag2/fuzzy_flag1   是否过于亮或者有噪声。
%                                       注意，这个是之前的一些设置，几年后重看，觉得这样反而可能导致过拟合，先暂且保留。
% ----------------------------------------------输出----------------------------------------------
%   P : 这一时刻LI（红）蛇的结果，即这一时刻的观测值
%   Q : 这一时刻MA（蓝）蛇的结果，即这一时刻的观测值
%   J : Binary image with the segmented region
%   fast_r/fast_b   是否发现LI/MA在某些帧中运动速度较快。

% -----------------------0、处理输入。-----------------------
defaultoptions=struct('Verbose',false,'Wline',0.9,'Wedge',10,'Wterm',2,'Sigma1',1,'Sigma2',1,'Alpha',4,'Beta',2,'Delta',0.1,'Gamma',0.0025,'Kappa',1,'Iterations',200,'GIterations',0,'Mu',0.2,'Sigma3',1);
if(~exist('Options','var'))
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);%
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end%
    if(length(tags)~=length(fieldnames(Options))) 
        warning('snake:unknownoption','unknown options found');
    end
end
% 图像变成double类型。
I = double(I);
[ylen, xlen]=size(I);
% 1 开始计算蛇的能量，用文章中的4式的第三个，处理输入的图像I，返回Eext（外能量函数，文章中的Gext）
% 和J（图像对y的导数，不过MATLAB里x和y是颠倒的，有点迷惑）。
[Eext,J] = ExternalForceImage2D_1(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1,Options.Verbose);
if(Options.Verbose)  
    figure,imshow(Eext,[]),title('Eextern')  
end
if(Options.Verbose)  
    figure,imshow(Eext,[]),title('Eextern with panelty')  
end
% 2 在偏离上一幅图较大的地方，加上较大的惩罚函数。不是很重要，文中没写。
Estereo=stereo(ylen,xlen,P,Q);
Eext=Eext+0.001*Estereo;
if(nargout>1)
     J1=DrawSegmentedArea2D(P,size(I));
end   
% 3 开始修正Eext，分别对LI和MA修正。
% 3.1 初始化两条蛇的Eext、Ggc或Gdc值、各种阈值。
volcano=200;  % 文章中的Ggc或者Gdc值。
Eext_red=Eext;
Eext_blue=Eext;
eliminate_x_force=0;
yi_red=fix(P(:,2));%LI红线
yi_blue=fix(Q(:,2));%MA蓝线。
red_dot=P(:,2);blue_dot=Q(:,2);%初始化，其实都弄成0也可以的。
for i=1:xlen
    % 一些初步的处理。
    % 首先分别为LI和MA设置一个“禁区”，在禁区内增大能量值，不要让二者的图像梯度互相影响。
    Eext_red(yi_blue(i)-fix(IMW(i)/2):ylen,i)=200;%把Eext_red中，原来LI-MA的中间线以下（靠近MA线）的部分，增大能量值，不要让LI蛇跑到这里来。
    Eext_blue(1:yi_red(i)+fix(IMW(i)/2),i)=200;%类似，把Eext_blue中，LI-MA的中间线上方的部分，增大能量值。
    % 然后，计算LI和MA每个点的导数，后面会将其用于Derivative constraint energy修正。
    if i~=xlen
        red_dot(i)=P(i+1,2)-P(i,2);%计算导数
        blue_dot(i)=Q(i+1,2)-Q(i,2);
    else
        red_dot(i)=red_dot(i-1);blue_dot(i)=blue_dot(i-1);
    end
end
%红蓝线上下推阈值，初始值都设为第一幅图的相应值。
push_r_up_thr=grayscale_rmax_up_mean;
push_r_down_thr=grayscale_rmax_down_mean;
push_b_up_thr=grayscale_bmax_up_mean;
push_b_down_thr=grayscale_bmax_down_mean;
% 下面四个，表示超出阈值（push_r_up_thr）多少，才需要加入Ggc。这个15 50 40 50什么的，都是看了好多个序列，才试出来的。
% 四年后再看，不但费时费力，而且还有可能过拟合，当时花这么多时间真心血亏。当时是不懂机器学习，现在再做，
% 如果能把它和强化学习结合起来，自动训练这些阈值，甚至自动设定别的阈值，应该会比现有工作强得多。
margin_r_push_up=15*ones(1,xlen);
margin_b_push_up=50*ones(1,xlen);
margin_r_push_down=40*ones(1,xlen);
margin_b_push_down=50*ones(1,xlen);
% 下面Iymax什么的，是用来找到LI/MA附近合理的灰度梯度最大值位置的。这个最大值位置可以用于辅助判断是否要加入Grayscale constraint energy。
% 最终用到的是Iymax_real_r和Iymax_real_b，其他的都是临时变量。
Iymax_r=zeros(1,xlen);Iymax_b=zeros(1,xlen);
Iymax_temp_r=zeros(1,xlen);Iymax_temp_b=zeros(1,xlen);
Iymax_real_r=zeros(1,xlen);Iymax_real_b=zeros(1,xlen);
Iymax_real_temp_r=zeros(1,xlen);Iymax_real_temp_b=zeros(1,xlen);
% 下面四个，是LI/MA灰度梯度最大值位置的、上下邻域内的、灰度均值。初始化为0。
mean_up_r=zeros(1,xlen);mean_up_b=zeros(1,xlen);
mean_down_r=zeros(1,xlen);mean_down_b=zeros(1,xlen);
% 下面四个，是输入的LI/MA位置的、上下邻域内的、灰度均值。初始化为0。
mean_red_up=zeros(1,xlen);mean_blue_up=zeros(1,xlen);
mean_red_down=zeros(1,xlen);mean_blue_down=zeros(1,xlen);
% 3.2 对于每一列，分别为LI/MA找到原来位置附近（search_l个像素点内）的、符合灰度要求的、灰度梯度最大值点的y坐标。
% 目的是根据灰度梯度最大值位置的信息，辅助判断是否要加入Grayscale constraint energy。
% 目前这些条件都是根据经验选取的，好处是看了很多的图，应该有信心比较鲁棒；但坏处是调参真心太枯燥。
% 限于篇幅，文中没有写这个辅助判断准则。
invalid_r_up=zeros(1,xlen);invalid_b_up=zeros(1,xlen);invalid_r_down=zeros(1,xlen);invalid_b_down=zeros(1,xlen);
for i=1:xlen
    search_l=min(8,fix(IMW(i)/2));
    Ji_r=J(yi_red(i)-search_l:yi_red(i)+search_l,i);  % LI
    Ji_b=J(yi_blue(i)-search_l:yi_blue(i)+search_l,i);  % MA
    Ji_rs=sort(Ji_r);%排序，排序后的Ji_rs中，最后一个是Ji_r中的最大值，即灰度梯度最大值。
    Ji_bs=sort(Ji_b);
    Jilen=size(Ji_r,1);
    for j=Jilen:-1:1%对于LI，从灰度梯度最大的开始试探，就看他是不是满足灰度条件
        k=find(Ji_r==Ji_rs(j));%Ji_rs(j)是最大值。k就是最大值的值的位置。
        if size(k,1)>1%如果不止一个最大值，那么，放弃这个限制条件。
            %因为这种情况，一般来讲都是LI往下跑得很快，跑到黑色区域里去了，所以就让Iymax_r(i)>search_l，让他直接下推。
            Iymax_r(i)=search_l+1;
            break;
        end
        Iymax_temp_r(i)=k;
        Iymax_real_temp_r(i)=Iymax_temp_r(i)+yi_red(i)-search_l;
        mean_up_r_temp=mean2(I(Iymax_real_temp_r(i)-neighbor_out_r(i)-1:Iymax_real_temp_r(i)-1,max(1,i-2):min(i+2,xlen)));
        mean_down_r_temp=mean2(I(Iymax_real_temp_r(i)+1:Iymax_real_temp_r(i)+neighbor_in(i)+1,max(1,i-2):min(i+2,xlen)));
        if abs(mean_up_r_temp-grayscale_rmax_up_mean(i))<abs(mean_up_r_temp-grayscale_bmax_up_mean(i)) ...
                && abs(mean_down_r_temp-grayscale_rmax_down_mean(i))<abs(mean_down_r_temp-grayscale_bmax_down_mean(i))...
                && abs(Q(i,2)-Iymax_real_temp_r(i))>abs(P(i,2)-Iymax_real_temp_r(i))
            % 这个if是判断这个灰度梯度最大值点，是否合理（如果不满足一些灰度条件，就可能不合理）
            Iymax_r(i)=k;%如果满足灰度条件，说明Ji_r中的第k个数确实就是LI的灰度梯度最大处。
            Iymax_real_r(i)=Iymax_r(i)+yi_red(i)-search_l;
            mean_up_r(i)=mean_up_r_temp;
            mean_down_r(i)=mean_down_r_temp;
            break;%如果满足灰度条件的话，那么，跳出for j=Jilen:-1:1这个循环。
        end
        %如果不满足灰度条件，那么j就减了1，找第二大的数，重新看是否满足灰度条件。   
        if j==1
            invalid_r_up(i)=1;invalid_r_down(i)=1;
            Iymax_r(i)=search_l+1;%就是搜索区域中间的那个点，相当于恢复初始位置了。
            Iymax_real_r(i)=P(i,2);%也恢复初始位置，本来应该是yi_red(i)+1的，但我们还是让他恢复成P好了。
        end%如果找到最后一个都没找到满足上面那个if的点，那就说明灰度梯度最大值的方法无效，这个点就不用他来判断了。
    end
    for j=Jilen:-1:1%类似上面，现在处理MA。
        k=find(Ji_b==Ji_bs(j));%Ji_bs(j)是最大值。这句话是说，在Ji_b中找到等于Ji_bs(j)，也就是最大值的值的位置。
        Iymax_temp_b(i)=k;
        Iymax_real_temp_b(i)=Iymax_temp_b(i)+yi_blue(i)-search_l;
        mean_up_b_temp=mean2(I(Iymax_real_temp_b(i)-neighbor_in(i)-1:Iymax_real_temp_b(i)-1,max(1,i-2):min(i+2,xlen)));
        mean_down_b_temp=mean2(I(Iymax_real_temp_b(i)+1:Iymax_real_temp_b(i)+neighbor_out_b(i)+1,max(1,i-2):min(i+2,xlen)));
        if abs(mean_up_b_temp-grayscale_bmax_up_mean(i))<abs(mean_up_b_temp-grayscale_rmax_up_mean(i)) ...
                && abs(mean_down_b_temp-grayscale_bmax_down_mean(i))<abs(mean_down_b_temp-grayscale_rmax_down_mean(i))...
                && abs(P(i,2)-Iymax_real_temp_b(i))>abs(Q(i,2)-Iymax_real_temp_b(i)) %20170502加了这个条件，要求距离红线更远。
            Iymax_b(i)=k;%如果灰度和蓝线的更接近，说明Ji_r中的第k个数确实就是蓝线的灰度梯度最大处。
            Iymax_real_b(i)=Iymax_b(i)+yi_blue(i)-search_l;
            mean_up_b(i)=mean_up_b_temp;
            mean_down_b(i)=mean_down_b_temp;
            break;%如果满足灰度条件的话，那么，跳出for j=Jilen:-1:1这个循环。
        end
        %如果不满足灰度条件，那么j就减了1，找第二大的数，重新看是否满足灰度条件。 
        if j==1
            invalid_b_up(i)=1;invalid_b_down(i)=1;
            Iymax_b(i)=search_l+1;
            Iymax_real_b(i)=Q(i,2);
        end%如果找到最后一个都没找到满足上面那个if的点，那就不找了，就让它等于蓝线原来位置吧。
    end
end
% 以下，修正灰度梯度最大值位置中的“鼓包”。
Iymax_real_r=smooth(Iymax_real_r,39,'sgolay',2);%图像滤波，体现大多数点的运动情况，而忽略掉少数错误点的运动情况。
Iymax_real_b=smooth(Iymax_real_b,39,'sgolay',2);
Iymax_r=smooth(Iymax_r,39,'sgolay',2);
Iymax_b=smooth(Iymax_b,39,'sgolay',2);
if sum(abs(red_dot_ref)>0.1)<=0.5*xlen && sum(abs(red_dot_ref)>0.2)<=0.25*xlen ...
        && ~(sum(red_dot_ref>0)>=0.8*xlen||sum(red_dot_ref<0)>=0.8*xlen)
    [Iymax_real_r,red_alert]=repair_bumpupanddown(I,Iymax_real_r,0.1,0.2,30,'red',[]);
else
    red_alert=[0 0 0 0];
end
if sum(abs(blue_dot_ref)>0.1)<=0.5*xlen && sum(abs(blue_dot_ref)>0.2)<=0.25*xlen
    [Iymax_real_b,~]=repair_bumpupanddown(I,Iymax_real_b,0.1,0.2,30,'blue',red_alert);
end
% 3.3 观察当前LI/MA曲线的斜率，以及移动情况，记录下来。
red_dot_front_is_neg=red_dot(1:7)<0;%判断红线前7个点斜率是否都是负数。如果这个red_dot_is_neg里的7个数都是1（斜率是个绝对值较大的负数），就说明前7个点的数值一直在下降（在matlab的图中是位置上升）。
blue_dot_front_is_neg=blue_dot(1:7)<0;%判断蓝线前7个点斜率是否都是负数。
red_dot_front_is_pos=red_dot(1:7)>0;%判断红线前7个点斜率是否都是正数。如果是，说明红线前7个点的数值一直在上升（在matlab的图中是位置下降），说明有可能红线前面往下（在matlab的图中是往上）翘了。
blue_dot_front_is_pos=blue_dot(1:7)>0;%判断蓝线前7个点斜率是否都是正数。
red_dot_tail_is_neg=red_dot(xlen-7:xlen)<0;%判断红线后7个点斜率是否都是负数。如果是，说明红线后7个点的数值一直在下降（在matlab的图中是位置上升），有可能在往下（在matlab的图中是往上）翘。
blue_dot_tail_is_neg=blue_dot(xlen-7:xlen)<0;%判断蓝线后7个点斜率是否都是负数。
red_dot_tail_is_pos=red_dot(xlen-7:xlen)>0;%判断红线后7个点斜率是否都是正数。如果是，说明红线后7个点的数值一直在上升（在matlab的图中是位置下降），有可能在往上（在matlab的图中是往下）翘。
blue_dot_tail_is_pos=blue_dot(xlen-7:xlen)>0;%判断蓝线后7个点斜率是否都是正数。
movement_r=Iymax_real_r-yi_red-1;%刚才那个大for循环里，都有那个if限制，限制红蓝线不要跑太远。这样大多数可以，但是对于1040R这样的动作比较大的就不行，所以这儿来个补充。
movement_b=Iymax_real_b-yi_blue-1;%
movement_r_count_down=sum(movement_r>5);%红线灰度梯度最大处 在 原来红线位置 下方 的点的个数
movement_r_count_down1=sum(movement_r>0);
movement_r_count_down2=sum(movement_r>7);
movement_b_count_down=sum(movement_b>5);
movement_b_count_down1=sum(movement_b>0);
movement_b_count_down2=sum(movement_b>7);
% 3.4 根据上面的灰度、斜率、运动幅度的数值、LI/MA蛇的距离，确定要把Ggc和Gdc加在何处、以调整Eext。
blue_suspect=0;
if (~eliminate_x_force)%只对红蓝线只在y方向运动的时候才有效
    for i=1:xlen
        % 首先，计算输入的LI/MA位置的、上下邻域内的、灰度均值。
        mean_red_up(i)=mean2(I(yi_red(i)-neighbor_out_r(i):yi_red(i),max(1,i-2):min(i+2,xlen)));% 上邻域，yi_red(i)不用-1，因为取了fix。
        mean_red_down(i)=mean2(I(yi_red(i)+1:yi_red(i)+1+neighbor_in(i),max(1,i-2):min(i+2,xlen)));% 下邻域，yi_red(i)要+1，因为fix得到的坐标是LI线上方的信息。
        mean_blue_up(i)=mean2(I(yi_blue(i)-neighbor_in(i):yi_blue(i),max(1,i-2):min(i+2,xlen)));
        mean_blue_down(i)=mean2(I(yi_blue(i)+1:yi_blue(i)+1+neighbor_out_b(i),max(1,i-2):min(i+2,xlen)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%思路：红线和蓝线对于噪声的处理方式并不一样。首先，如果红线没有噪声、蓝线比较清晰，%%%%%%%%
        %%%%%%就都设置比较大的margin_r或者margin_b，尽量不要外加Ggc去影响他们。然后，如果LI线%%%%%%%%%
        %%%%%%有噪声，那么，我们知道应该让他容易向下运动、不容易向上运动，所以，可以直接改变阈值%%%%%%%
        %%%%%%让他容易向下运动，而难以向上运动；然而，如果蓝线不清晰的话，我们并不确定到底应该让他%%%%%%
        %%%%%%往上运动还是往下，所以要减小上推阈值，增大下推阈值，让他既容易向上，又容易向下，而%%%%%%%%
        %%%%%%到底是上还是下，靠灰度梯度最大值的位置（Iymax_b(i)与search_l+1的比较）的判断来决%%%%%%%%
        %%%%%%定往上还是往下。也就是说，灰度值本身的影响被削弱，而灰度梯度位置决定运动的方向。 %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3.4.1加入Ggc让LI向上运动：初衷是，如果红线上方灰度值太大（文中是差异值超过了20%），就加Ggc。
        % 后续处理中，加入了噪声、灰度梯度等的辅助判断条件。
        if ~ismember(i,exists_noise)
            % 通过判断噪声情况，改变加入Ggc的阈值。
            if noise_flag0==1%如果整体有噪声的话，那么，就算这个点没噪声，也把margin_r_push_up(i)增大一些，别乱往上推。
                margin_r_push_up(i)=20;
            end%
            if mean_red_up(i)>20%如果现在上方平均值大了（而一开始的是被判断为没有噪声），那么就很可疑，有可能是一开始没噪声，现在有了。
                margin_r_push_up(i)=30;%其实就是别让他推了。除非真的mean_red_up(i)能大于30。。
            end
            %如果需要，可以加一下，如果总体有噪声和没有噪声，怎么办。
        end
        if(ismember(i,exists_noise) && ~ismember(i,medium_noise) && ~ismember(i,severe_noise))
            push_r_up_thr(i)=grayscale_rmax_down_mean(i);%对于第i个点，有噪声，就用down的灰度来判断红线是否往上推。
            %margin_r_push_up(i)保留初始值10;
            if noise_flag0==1%如果整体有噪声的话，那么，这个点把margin_r_push_up(i)增大一些，别乱往上推。
                margin_r_push_up(i)=20;
            end
        end
        if(ismember(i,medium_noise))
            push_r_up_thr(i)=grayscale_rmax_down_mean(i);
            %margin_r_push_up(i)保留初始值10;
            if noise_flag0==1%如果整体有噪声的话，那么，这个点把margin_r_push_up(i)增大一些，别乱往上推。
                margin_r_push_up(i)=20;
            end
        end
        if(ismember(i,severe_noise))
            push_r_up_thr(i)=grayscale_rmax_down_mean(i);
            %margin_r_push_up(i)保留初始值10;
            if noise_flag0==1%如果整体有噪声的话，那么，这个点把margin_r_push_up(i)增大一些，别乱往上推。
                margin_r_push_up(i)=20;
            end
        end
        if mean_red_up(i)>=push_r_up_thr(i)+margin_r_push_up(i) && (invalid_r_up(i)==1 || Iymax_r(i)<=search_l && mean_down_r(i)>mean_up_r(i))...  %|| (~ismember(i,exists_noise) && sum(I(yi_red(i),max(1,i-5):min(i+5,xlen))>=15)>=10 && I(yi_red(i),i)>=15 && I(yi_red(i)-1,i)>=15 && I(yi_red(i)-2,i)>=15) ...
                || (~ismember(i,exists_noise) && movement_r(i)<-2 && I(yi_red(i),i)>=15 && I(yi_red(i)-1,i)>=15 && I(yi_red(i)-2,i)>=15) ...
                || grayscale_rmax_up_mean(i)<=20 && min(I(yi_red(i),max(1,i-10):min(i+10,xlen)))>=30 && min(I(yi_red(i)-2:yi_red(i),i))>=30 ...
                || movement_r(i)<-3.5 && I(yi_red(i),i)>=50 && I(yi_red(i)-1,i)>=50 && I(yi_red(i)-2,i)>=50 ...
                || (mean_red_up(i)>=10&&mean_red_down(i)<10&&inside_IMC_is_dark(i)==1)
            % 有一个概念，这个grayscale_rmax_up_mean(i)+margin_r_push_up(i)取得越大，
            % 就越难使得这三个都大于这个数，也就是越难通过Ggc让蛇往上运动。
            Eext_red(yi_red(i):yi_red(i)+7,i)=Eext_red(yi_red(i):yi_red(i)+7,i)+volcano;
            Eext_red(yi_red(i):yi_red(i)+7,max(1,i-1):min(i+1,xlen))=Eext_red(yi_red(i):yi_red(i)+7,max(1,i-1):min(i+1,xlen))+volcano/2;
            Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+1.5*volcano;
        end
        % 3.4.2加入Ggc让LI向下运动：初衷是，如果红线下方灰度值太小（文中是差异值超过了20%），就下推。
        % 后续处理中，类似于前面，也是加入了噪声、灰度梯度等的辅助判断条件。
        if(~ismember(i,exists_noise))
            if noise_flag0==0%如果总体没有噪声，那么，下推阈值进一步减小，让他不要下推。
                push_r_down_thr(i)=5;
                margin_r_push_down(i)=0;%如果没噪声的话，除非mean_red_down(i)<push_r_down_thr(i)-margin_r_push_down(i)=10才下推。
            end
            if noise_flag0==1%如果总体有噪声，阈值情况。
                if ~ismember(i,[1:20 xlen-19:xlen])  % 是否为最前或者最后几个点（应该是当初观察这几个点的灰度与众不同）
                    push_r_down_thr(i)=15;
                else
                    push_r_down_thr(i)=10;
                end
                margin_r_push_down(i)=0;%如果没噪声的话，除非mean_red_down(i)<push_r_down_thr(i)-margin_r_push_down(i)=10才下推。
            end
        end
        if(ismember(i,exists_noise) && ~ismember(i,medium_noise) && ~ismember(i,severe_noise))%
            if grayscale_rmax_down_mean(i)>=30
                push_r_down_thr(i)=grayscale_rmax_up_mean(i);
                %如果下方太亮，那就用grayscale_rmax_up_mean去推；否则用默认的grayscale_rmax_down_mean
                %（对应于，如果上方太亮，即有噪声，就用下面的推。然后呢，这儿和上推不一样的是，如果上方有噪声，要让他略容易下推，所以margin_r_push_down是在变的。）
            end
            margin_r_push_down(i)=15;%下方太亮（大于30）的，就是原来上方的减15；否则就是原来下方的（小于30了）减掉15，差不多是15以下。因为你的噪声阈值定的是15，所以现在是，如果下方小于15左右了，那就下推吧。
            %其他情况，invalid_r(i)预留为0。
            %可以预留一下，如果noise_flag0,noise_flag1==1，该怎么办.
        end
        if(ismember(i,medium_noise))
            if grayscale_rmax_down_mean(i)>=30
                push_r_down_thr(i)=grayscale_rmax_up_mean(i);
                %如果下方太亮，那就用grayscale_rmax_up_mean去推；否则用默认的grayscale_rmax_down_mean
            end
            margin_r_push_down(i)=10;%如果是中等噪声，那么，下推阈值改成10。即，下方太亮的，就是上方的减10，否则是下方的减10，大概20以下，容易下推一些。
            %其他情况，invalid_r(i)预留为0。
            %可以预留一下，如果noise_flag0,noise_flag1==1，该怎么办.
        end
        if(ismember(i,severe_noise))
            if grayscale_rmax_down_mean(i)>=30
                push_r_down_thr(i)=grayscale_rmax_up_mean(i);
                %如果下方太亮，那就用grayscale_rmax_up_mean去推；否则用默认的grayscale_rmax_down_mean
            end
            margin_r_push_down(i)=5;%严重噪声
            %其他情况，invalid_r(i)预留为0。
            %可以预留一下，如果noise_flag0,noise_flag1==1，该怎么办.
        end
        if mean_red_down(i)<=max(5, push_r_down_thr(i)-margin_r_push_down(i)) &&~(mean_red_down(i)<=15&&mean_red_up(i)<=15&&max(J(yi_red(i):yi_blue(i)-fix(IMW(i)/2),i))<3)... %
                || mean_red_down(i)<=grayscale_rmax_down_mean(i) && movement_r_count_down1>xlen*(300/361) && movement_r(i)>2 && mean(movement_r)>1.5 ...%170525把300改成了xlen*(300/361) %&& ~(min(I(yi_red(i),max(1,i-20):min(i+20,xlen)))>20)
                && (invalid_r_down(i)==1 || Iymax_r(i)>=search_l && mean_down_r(i)>mean_up_r(i))...
                || (I(yi_red(i),i)<=5 && I(yi_red(i)+1,i)<=5 && I(yi_red(i)+2,i)<=5) ...
                && inside_IMC_is_dark(i)==0
        %红线往下运动的各种条件总结。
            Eext_red(yi_red(i)-7:yi_red(i),i)=Eext_red(yi_red(i)-7:yi_red(i),i)+volcano;
            Eext_red(yi_red(i)-7:yi_red(i),max(1,i-1):min(i+1,xlen))=Eext_red(yi_red(i)-7:yi_red(i),max(1,i-1):min(i+1,xlen))+volcano/2; 
            Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+1.5*volcano;
        end
        %3.4.3加入Ggc让MA向上运动：每个点都确定是否清晰。越清晰，就越没必要用用外加的力。
        if(~ismember(i,exists_fuzzy))
            push_b_up_thr(i)=max(210,grayscale_bmax_up_mean(i));%原来说的也是，只要蓝线清晰，就不要往上推了就得了。
            %margin_b_push_up(i)保留初始值，即50;
        end
        if(ismember(i,exists_fuzzy) && ~ismember(i,medium_fuzzy) &&~ismember(i,severe_fuzzy))
            %不清晰的部分，动作阈值继续改小，让他略容易往上运动。
            margin_b_push_up(i)=40;
        end
        if(ismember(i,medium_fuzzy))
            %push_b_up_thr(i)保留初始值，即grayscale_bmax_up_mean(i)
            margin_b_push_up(i)=35;
        end
        if(ismember(i,severe_fuzzy))
            %push_b_up_thr(i)保留初始值，即grayscale_bmax_up_mean(i)
            margin_b_push_up(i)=30;
        end
        if mean_blue_up(i)>=push_b_up_thr(i)+margin_b_push_up(i) && (invalid_b_up(i)==1 || Iymax_b(i)<=search_l && mean_down_b(i)>mean_up_b(i))
        %如果蓝线上方三个像素点平均值有点过于的白（大于push_b_up_thr(i)+margin_b），就认为蓝线跑到了白色的组织里面。
        %push_b_up_thr(i)+margin_b_push_up(i)取得越大，就越难往上推。
        %Iymax_b(i)<search_l+1表示，灰度梯度最大值点要在原来的蓝线蛇位置上方（<search_l+1），才往上推。mean_down_b(i)>mean_up_b(i)不解释。
            Eext_blue(yi_blue(i):yi_blue(i)+7,i)=Eext_blue(yi_blue(i):yi_blue(i)+7,i)+volcano;
            Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))=Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))+volcano/2; %惠及周围
            Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+1.5*volcano;
        end
        if mean_blue_up(i)>mean_blue_down(i) && ~(max(I(yi_blue(i)+1:yi_blue(i)+10,i))>=max(I(yi_blue(i)-neighbor_in(i):yi_blue(i),i))) ...%20170502
                && ~(grayscale_bmax_up_mean(i)>grayscale_bmax_down_mean(i)) ...
                && ~(sum(grayscale_bmax_up_mean(max(1,i-50):min(xlen,i+50))>150)>10) ...
                && inside_IMC_is_bright(i)==0 ...
                && max(abs(blue_dot(max(1,i-20):min(i+20,xlen))))>=0.05 ... 
                &&~(sum(blue_dot(max(1,i-20):i)<0)>10&&sum(blue_dot(i:min(i+20,xlen))>0)>10)
            Eext_blue(yi_blue(i):yi_blue(i)+7,i)=Eext_blue(yi_blue(i):yi_blue(i)+7,i)+2*volcano;
            Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))=Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))+1.5*volcano; %惠及周围
            Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+4*volcano;
        end
        %3.4.4加入Ggc让MA向下运动：类似于上面，通过灰度斜率等条件，判断阈值的变化。
        if(~ismember(i,exists_fuzzy))
            if grayscale_bmax_down_mean(i)>=180%应该有一大堆大于180的，就像红线下方也有一大堆大于30的。
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %如果下方太亮，那就用grayscale_bmax_up_mean去推；否则用默认的grayscale_bmax_down_mean
                %和红线下推很像。
                margin_b_push_down(i)=10;%如果这样的话，把margin_b_push_down(i)改成10，毕竟已经用上面的灰度去推了。。。
            end
            if mean_blue_down(i)<=110 && std2(I(yi_blue(i)-neighbor_in(i):yi_blue(i),max(1,i-2):min(i+2,xlen)))<=20%如果蓝线一开始是清晰的，但是现在不清晰了，那就比较可疑。类似于红线上推，只不过多加了个灰度值比较小的条件。
                blue_suspect=1;
            end%如果蓝线可疑，这张图的蓝线步长除以2。。
        end
        if(ismember(i,exists_fuzzy) && ~ismember(i,medium_fuzzy) &&~ismember(i,severe_fuzzy))%不清晰的部分，动作阈值继续改小。
            if grayscale_bmax_down_mean(i)>=180
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %如果下方太亮，那就用grayscale_bmax_up_mean去推；否则用默认的grayscale_bmax_down_mean
                %和红线下推很像。
                margin_b_push_down(i)=10;%如果这样的话，把margin_b_push_down(i)改成10，毕竟已经用上面的灰度去推了。。。
            end
        end
        if(ismember(i,medium_fuzzy))
            if grayscale_bmax_down_mean(i)>=180
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %如果下方太亮，那就用grayscale_bmax_up_mean去推；否则用默认的grayscale_bmax_down_mean
                %和红线下推很像。
                margin_b_push_down(i)=10;%如果这样的话，把margin_b_push_down(i)改成10，毕竟已经用上面的灰度去推了。。。其实10都未必行，弄不好得-10。。。
            end
        end
        if(ismember(i,severe_fuzzy))
            if grayscale_bmax_down_mean(i)>=180
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %如果下方太亮，那就用grayscale_bmax_up_mean去推；否则用默认的grayscale_bmax_down_mean
                %和红线下推很像。
                margin_b_push_down(i)=10;%如果这样的话，把margin_b_push_down(i)改成10，毕竟已经用上面的灰度去推了。。。
            end
        end
        if (mean_blue_down(i)<=push_b_down_thr(i)-margin_b_push_down(i) || mean_blue_down(i)<=grayscale_bmax_down_mean(i) && movement_b_count_down1>xlen*(300/361) && movement_b(i)>2 && mean(movement_b)>1.5)...
                && (invalid_b_down(i)==1 || Iymax_b(i)>=search_l && mean_down_b(i)>mean_up_b(i))
            Eext_blue(yi_blue(i)-7:yi_blue(i),i)=Eext_blue(yi_blue(i)-7:yi_blue(i),i)+volcano;
            Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-1):min(i+1,xlen))=Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-1):min(i+1,xlen))+volcano/2;
            Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+1.5*volcano;
        end
        % 以上，用灰度判断是否加入Ggc的部分结束。下面，加入红蓝线贴近的修正程序。
        % 3.4.5如果LI/MA只差了2像素点以内，就可能有问题。根据局部灰度复核一下，如果有问题，就加入Ggc调整蛇位置。
        tol=min(IMW(i),6);
        if (Q(i,2)-P(i,2)<=tol)
            yi_mid=fix((yi_red(i)+yi_blue(i))/2);
            I_local=I(yi_mid-2:yi_mid+4,max(P(i,1)-2,1):min(P(i,1)+2,xlen));%这是为了防止索引出现0或者大于xlen的情况。
            mean_local=mean(mean(I_local));%局部平均灰度
            volcano1=200*(tol/(Q(i,2)-P(i,2)));
            if abs(mean_local-grayscale_rmax_down_mean(i))<10 || abs(mean_local-(grayscale_bmax_up_mean(i)/2+grayscale_bmax_down_mean(i)/2))>10 && mean_blue_down(i)<=grayscale_bmax_down_mean(i)
                %此if处理蓝线向上的情况
                if (ismember(i,1:7) && min(abs(blue_dot_front_is_pos))==1) || ...
                        (ismember(i,xlen-7:xlen) && min(abs(blue_dot_tail_is_neg))==1)
                    volcano1=volcano1*2;%如果是蓝线左边或者右边往上翘，加大Ggc往下驱赶蛇。
                end
                if Q(i,2)-P(i,2)<=tol*0.8%如果太小了，那就直接生把蓝线往下拉。
                    Q(i,2)=Q(i,2)+tol/2;
                    disp('blue pushed up manually because it is getting too close to red')
                    disp(fr)
                end
                Eext_blue(yi_blue(i)-7:yi_blue(i),i)=Eext_blue(yi_blue(i)-7:yi_blue(i),i)+2*volcano1;
                Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+4*volcano1;%这次不要惠及周围了，因为下一个点还会循环到这儿来的。你要是惠及周围了，相当于加了很多次火山。
            end
            if abs(mean_local-grayscale_bmax_up_mean(i))<10 || abs(mean_local-(grayscale_rmax_up_mean(i)/2+grayscale_rmax_down_mean(i)/2))>10 && mean_red_up(i)>max(10,grayscale_rmax_up_mean(i))
                %此if处理红线向下的情况
                if (ismember(i,1:7) && min(abs(red_dot_front_is_neg))==1) || ...
                        (ismember(i,xlen-7:xlen) && min(abs(red_dot_tail_is_pos))==1)
                    volcano1=volcano1*2;%如果是红线左边或者右边往下耷拉，妈的，加大火山往上推。不信tmd治不了他了。
                end
                if Q(i,2)-P(i,2)<=tol*0.8%如果太小了，那就直接生把红线往上拉.不过，先拉一半再加火山吧，如果拉太多了，怕出事儿。
                    P(i,2)=P(i,2)-tol/2;
                    disp('red pushed up manually because it is getting too close to blue')
                    disp(fr)
                end
                if P(i,2)-yi_red(i)>0.8, trick=1;else trick=0;end%如果红线正好是***.9的话，那么，就在他下面的一个点再开始加火山。否则可能反而给往下推了。1060R第88张。
                %蓝线不用加这个东西，因为是往下推，在上面加火山，本来那个yi_blue就把小数点抹掉了。
                Eext_red(yi_red(i)+trick:yi_red(i)+trick+7,i)=Eext_red(yi_red(i)+trick:yi_red(i)+trick+7,i)+2*volcano1;
                Eext_red(yi_red(i)+trick,i)=Eext_red(yi_red(i)+trick,i)+4*volcano1;
            end
        end
        % 3.4.6加入鼓包/误差修正程序
        % 用后面的图和第一张图每一点的斜率去对比，如果差太多了，那就加Gdc驱赶。实际做的时候，第一幅图的斜率不是0就是1，所以用第二幅图的作为判定标准。
        if fr~=1
            vol_times=10;%主要是那个斜率平均值太小了，一般是0.几甚至0.0几，所以这儿加大一些Gdc的倍数。
            mean5_abs_red_dot_ref=mean(abs(red_dot_ref(max(1,i-5):min(i+5,xlen))));%第i点附近的、第2张图的红蓝线斜率绝对值的均值。
            mean5_abs_blue_dot_ref=mean(abs(blue_dot_ref(max(1,i-5):min(i+5,xlen))));
            mean10_abs_red_dot_ref=mean(abs(red_dot_ref(max(1,i-10):min(i+10,xlen))));%第i点附近的、第2张图的红蓝线斜率绝对值的均值。
            mean10_abs_blue_dot_ref=mean(abs(blue_dot_ref(max(1,i-10):min(i+10,xlen))));
            if (i>7 && i<xlen-7)%前7个点和后7个点暂且不管，后面单独考虑。
                if (red_dot(i)<-0.01 && red_dot(i+1)>0.01)%LI线出现极小值（图中表示向上的鼓包）。注意，此时不能直接就给他加Gdc，得再判断斜率绝对值，见下一个if。
                    %本来是0，但为了防止有些1*10^-4这样的值来混淆视听，后来把0全都改成0.01了。
                    mean5_abs_red_dot=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_red_dot=mean(abs(red_dot(max(1,i-10):min(i+10,xlen))));
                    if (mean5_abs_red_dot>mean5_abs_red_dot_ref+0.05 && sum(Iymax_real_r(i-5:i+5)>P(i-5:i+5,2))>=8) ...
                            ||(mean10_abs_red_dot>mean10_abs_red_dot_ref+0.05 && sum(Iymax_real_r(max(1,i-10):min(i+10,xlen))>P(max(1,i-10):min(i+10,xlen),2))>=16) ...
                            || sum(Iymax_real_r(i-5:i+5)>P(i-5:i+5,2)+7)>=8
                        % 用鼓包附近的斜率（用了两种不同邻域内的斜率）、灰度梯度最大值点，这些因素的综合考虑，判断是否增加Gdc。
                        % 这些手工特性，需要考虑能否用自动特征代替。
                        for try_loca=1:30%这个for循环是在找鼓包的长度
                            if mean(abs(red_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(red_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;%i-forced_area:i+forced_area这个范围内，加Gdc。
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))-red_dot_ref(max(1,i-5):min(i+5,xlen))));
                        volcano1=volcano*vol_times*scalarco;%加Gdc，这个平均灰度绝对值越大，Gdc就越大。
                        Eext_red(yi_red(i)-7:yi_red(i),max(1,i-forced_area):min(i+forced_area,xlen))=Eext_red(yi_red(i)-7:yi_red(i),max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
                if (red_dot(i)>0.01 && red_dot(i+1)<-0.01)%LI线出现极大值（图中表示向下的鼓包）
                    mean5_abs_red_dot=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_red_dot=mean(abs(red_dot(max(1,i-10):min(i+10,xlen))));
                    if (mean5_abs_red_dot>mean5_abs_red_dot_ref+0.05 && sum(Iymax_real_r(i-5:i+5)<P(i-5:i+5,2))>=8) ...
                            ||(mean10_abs_red_dot>mean10_abs_red_dot_ref+0.05 && sum(Iymax_real_r(max(1,i-10):min(i+10,xlen))<P(max(1,i-10):min(i+10,xlen),2))>=8)...
                            || sum(Iymax_real_r(i-5:i+5)<P(i-5:i+5,2)-7)>=8 && lightflag==0
                        %逻辑和上面一样。
                        for try_loca=1:30
                            if mean(abs(red_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(red_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;%i-forced_area:i+forced_area这个范围内，加Gdc。
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))-red_dot_ref(max(1,i-5):min(i+5,xlen))));
                        volcano1=volcano*vol_times*scalarco;
                        Eext_red(yi_red(i):yi_red(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))=Eext_red(yi_red(i):yi_red(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
                if (blue_dot(i)<-0.01 && blue_dot(i+1)>0.01)%MA线出现极小值（图中表示向上的鼓包）。
                    mean5_abs_blue_dot=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_blue_dot=mean(abs(blue_dot(max(1,i-10):min(i+10,xlen))));
                    if lightflag==1, mean5_abs_blue_dot=mean5_abs_blue_dot+0.1;mean10_abs_blue_dot=mean10_abs_blue_dot+0.1;end%我就是想让1055L1有往上的鼓包，就给我推下去。但是这个还不好，因为有的推多了。先不管了！！！！
                    if (mean5_abs_blue_dot>mean5_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(i-5:i+5)>Q(i-5:i+5,2))>=8) ...
                            || sum(Iymax_real_b(i-5:i+5)>Q(i-5:i+5,2)+7)>=8 ...
                            || (mean10_abs_blue_dot>mean10_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(max(1,i-10):min(i+10,xlen))>Q(max(1,i-10):min(i+10,xlen),2))>=16)
                        %同上,用鼓包附近的斜率（用了两种不同邻域内的斜率）、灰度梯度最大值点，这些因素的综合考虑，判断是否增加Gdc。
                        for try_loca=1:30%应该不会有超过30个像素点的鼓包了。。。
                            if mean(abs(blue_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(blue_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;%i-forced_area:i+forced_area这个范围内，加Gdc。
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))-blue_dot_ref(max(1,i-5):min(i+5,xlen))));
                        if lightflag==1, vol_times=vol_times*4;end
                        volcano1=volcano*vol_times*scalarco;
                        Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-forced_area):min(i+forced_area,xlen))=Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
                if (blue_dot(i)>0.01 && blue_dot(i+1)<-0.01)%MA线出现极大值（图中表示向下的鼓包）。
                    mean5_abs_blue_dot=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_blue_dot=mean(abs(blue_dot(max(1,i-10):min(i+10,xlen))));
                    if (mean5_abs_blue_dot>mean5_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(i-5:i+5)<Q(i-5:i+5,2))>=8)  ...
                            ||(mean10_abs_blue_dot>mean10_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(max(1,i-10):min(i+10,xlen))<Q(max(1,i-10):min(i+10,xlen),2))>=16)...
                            || sum(Iymax_real_b(i-5:i+5)<Q(i-5:i+5,2)-7)>=8 && lightflag==0
                    % 用鼓包附近的斜率（用了两种不同邻域内的斜率）、灰度梯度最大值点，这些因素的综合考虑，判断是否增加Gdc。
                        for try_loca=1:30
                            if mean(abs(blue_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(blue_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))-blue_dot_ref(max(1,i-5):min(i+5,xlen))));
                        volcano1=volcano*vol_times*scalarco;
                        Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))=Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
            end
            v1=5;
            % 下面是类似的处理鼓包的程序，只不过是考虑前面7个点和最后7个点的。
            if ismember(i,1:7)
                m_rd_f=mean(red_dot(1:7));m_rd_f_r=mean(red_dot_ref(1:7));
                m_bd_f=mean(blue_dot(1:7));m_bd_f_r=mean(blue_dot_ref(1:7));
                if min(abs(red_dot_front_is_neg))==1 && abs(m_rd_f)>0.1 && ~(m_rd_f_r<-0.1) || (m_rd_f_r<-0.2 && abs(m_rd_f)>=abs(m_rd_f_r)*1.5)
                    %那个if条件的分析：如果前7个点红线斜率都是负的，且绝对值平均值大于0.1，就说明红线前面7个点可能往上（在matlab的图中是往下）翘，需要用Gdc修正。
                    %但是，如果第二幅图前7个点斜率平均值也小于-0.1（就是斜率都是负的，且绝对值大于0.1），那么，就认为是正常的了。
                    volcano1=volcano*v1*max(0,mean(red_dot_ref(1:7)-red_dot(1:7)));%如果第一幅图比较平缓（斜率没负那么多）。
                    if sum(Iymax_real_r(1:7)<P(1:7,2)-5)==8,volcano1=volcano1*2;end%如果发现，前8个点中每一个点的红线梯度最大值处都在原来红线上方5以上，那么，Gdc加倍。
                    if movement_r_count_down1<=xlen*(180/361) %如果至少一半的点的位移小于0（往上移了）。
                        %这个if是为了防止有的时候LI/MA线整体往下移，但是也满足“min(abs(red_dot_front_is_neg))==1...”这个if条件，
                        %而导致既在上面加Gdc，又在下面加Gdc的情况（这会使得蛇在超大能量的趋势下开始开始乱动）。
                        %movement_r_count1是“LI线梯度最大值处”在LI线下方的点数。如果这个点数小于一半，就证明在这张图中，
                        %LI线整体没有发生明显向下的移动，才放心地加Gdc。
                        Eext_red(yi_red(i):yi_red(i)+4,i)=Eext_red(yi_red(i):yi_red(i)+4,i)+volcano1;%下方加Gdc。
                        Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                    end
                end
                if min(abs(blue_dot_front_is_neg))==1 && abs(m_bd_f)>0.1 && ~(m_bd_f_r<-0.1) || (m_bd_f_r<-0.2 && abs(m_bd_f)>=abs(m_bd_f_r)*1.5)%如果这样，就说明蓝线前面7个点往上（在matlab的图中是往下）翘了
                    volcano1=volcano*v1*max(0,mean(blue_dot_ref(1:7)-blue_dot(1:7)));%
                    if sum(Iymax_real_b(1:7)<Q(1:7,2)-5)==8,volcano1=volcano1*2;end%
                    if movement_b_count_down1<=xlen*(180/361) 
                        %和上面类似，如果至少一半的点的位移小于0（往上移了），即MA线没有明显地向下移动，才放心地加Gdc。
                        Eext_blue(yi_blue(i):yi_blue(i)+4,i)=Eext_blue(yi_blue(i):yi_blue(i)+4,i)+volcano1;
                        Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                    end
                end
                if min(abs(red_dot_front_is_pos))==1 && abs(m_rd_f)>0.1 && ~(m_rd_f_r>0.1) || (m_rd_f_r>0.2 && abs(m_rd_f)>=abs(m_rd_f_r)*1.5)%如果这样，就说明红线前面7个点往下（在matlab的图中是往上）翘了
                    volcano1=volcano*v1*max(0,mean(red_dot(1:7)-red_dot_ref(1:7)));%这个时候应该反过来减。
                    if sum(Iymax_real_r(1:7)>P(1:7,2)+5)==8,volcano1=volcano1*2;end%
                    %这个不再加if movement_r_count1<=180这样的条件了，主要是没发现LI/MA蛇明显快速地往上移的。
                    Eext_red(yi_red(i)-4:yi_red(i),i)=Eext_red(yi_red(i)-4:yi_red(i),i)+volcano1;
                    Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                end
                if min(abs(blue_dot_front_is_pos))==1 && abs(m_bd_f)>0.1 && ~(m_bd_f_r>0.1) || (m_bd_f_r>0.2 && abs(m_bd_f)>=abs(m_bd_f_r)*1.5)%如果这样，就说明红线前面7个点往下（在matlab的图中是往上）翘了
                    volcano1=volcano*v1*max(0,mean(blue_dot(1:7)-blue_dot_ref(1:7)));%这个时候应该反过来减。
                    if sum(Iymax_real_b(1:7)>Q(1:7,2)+5)==8,volcano1=volcano1*2;end%
                    %这个不再加if movement_b_count1<=180这样的条件了，主要是没发现LI/MA蛇明显快速地往上移的。
                    Eext_blue(yi_blue(i)-4:yi_blue(i),i)=Eext_blue(yi_blue(i)-4:yi_blue(i),i)+volcano1;
                    Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                end
            end
            if ismember(i,xlen-7:xlen)  % 这个是最后7个点的情况，和最前7个点的情况一样的。
                m_rd_l=mean(red_dot(xlen-7:xlen));m_rd_l_r=mean(red_dot_ref(xlen-7:xlen));
                m_bd_l=mean(blue_dot(xlen-7:xlen));m_bd_l_r=mean(blue_dot_ref(xlen-7:xlen));
                if min(abs(red_dot_tail_is_neg))==1 && abs(m_rd_l)>0.1 && ~(m_rd_l_r<-0.1) || (m_rd_l_r<-0.2 && abs(m_rd_l)>=abs(m_rd_l_r)*1.5)%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %如果满足：①最后7个点斜率都小于0，且绝对值平均值大于0.1；②第二幅图后7个点斜率平均值并不小于-0.1，
                    %或者，这一幅图后7个点斜率平均值，大于第二幅图后7个点斜率平均值的1.5倍，且第二幅图后7个点斜率平均值小于-0.2（比较斜）
                    %就说明LI线后面7个点在matlab的图中往上翘（坐标偏小）了
                    volcano1=volcano*v1*max(0,mean(red_dot_ref(xlen-7:xlen)-red_dot(xlen-7:xlen)));%调整Gdc值
                    if sum(Iymax_real_r(xlen-7:xlen)>P(xlen-7:xlen,2)+5)==8,volcano1=volcano1*2;end
                    %如果发现，最后8个点中每一个点的红线梯度最大值处都在原来红线下方5以下，那么，Gdc加倍。
                    Eext_red(yi_red(i)-4:yi_red(i),i)=Eext_red(yi_red(i)-4:yi_red(i),i)+volcano1;
                    Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                end
                if min(abs(blue_dot_tail_is_neg))==1 && abs(m_bd_l)>0.1 && ~(m_bd_l_r<-0.1) || (m_bd_l_r<-0.2 && abs(m_bd_l)>=abs(m_bd_l_r)*1.5)
                    %如果这样，就说明蓝线后面7个点在matlab的图中往上翘（坐标偏小）了
                    volcano1=volcano*v1*max(0,mean(blue_dot_ref(xlen-7:xlen)-blue_dot(xlen-7:xlen)));
                    if sum(Iymax_real_b(xlen-7:xlen)>Q(xlen-7:xlen,2)+5)==8,volcano1=volcano1*2;end
                    Eext_blue(yi_blue(i)-4:yi_blue(i),i)=Eext_blue(yi_blue(i)-4:yi_blue(i),i)+volcano1;
                    Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                end
                if min(abs(red_dot_tail_is_pos))==1 && abs(m_rd_l)>0.1 && ~(m_rd_l_r>0.1) || (m_rd_l_r>0.2 && abs(m_rd_l)>=abs(m_rd_l_r)*1.5)
                    volcano1=volcano*v1*max(0,mean(red_dot(xlen-7:xlen)-red_dot_ref(xlen-7:xlen)));
                    if sum(Iymax_real_r(xlen-7:xlen)<P(xlen-7:xlen,2)-5)==8,volcano1=volcano1*2;end%
                    if movement_r_count_down1<=xlen*(180/361) %确定没发生LI线整体往下移的情况
                        Eext_red(yi_red(i):yi_red(i)+4,i)=Eext_red(yi_red(i):yi_red(i)+4,i)+volcano1;
                        Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                    end
                end
                if min(abs(blue_dot_tail_is_pos))==1 && abs(m_bd_l)>0.1 && ~(m_bd_l_r>0.1) || m_bd_l_r>0.2 && abs(m_bd_l)>=abs(m_bd_l_r)*1.5
                    %如果这样，就说明蓝线后面7个点往上（在matlab的图中是往下）翘了
                    volcano1=volcano*v1*max(0,mean(blue_dot(xlen-7:xlen)-blue_dot_ref(xlen-7:xlen)));%调整Gdc值
                    if sum(Iymax_real_b(xlen-7:xlen)<Q(xlen-7:xlen,2)-5)==8,volcano1=volcano1*2;end%
                    if movement_b_count_down1<=xlen*(180/361) %确定没发生MA线整体往下移的情况
                        Eext_blue(yi_blue(i):yi_blue(i)+4,i)=Eext_blue(yi_blue(i):yi_blue(i)+4,i)+volcano1;
                        Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                    end
                end
            end
        end
    end
    % 3.4.7 受到上面启发，加入整体运动程序
    fl_r=0;fl_b=0;r=0;%每一张图都置零，如果这一张图动作大，才置一。
    % 先是处理LI向下移动的情况。
    if movement_r_count_down>xlen*(180/361) && movement_r_count_down2>xlen*(100/361)  
        % 上面的180和100也是试出来的。就是如果有180个都往下移动了（5个像素点），且有100个往下移动较大（7个像素点），才加入Ggc。
        P(:,2)=P(:,2)+abs(mean(movement_r));%另外一个思路，直接生往下推
        for i=1:xlen%%
            Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)=Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)+volcano;
            Eext_red(fix(P(i,2)),i)=Eext_red(fix(P(i,2)),i)+volcano;
        end
        disp('red moves down fast'),disp(fr)
        fl_r=1;fast_r=1;
    else
        if movement_r_count_down1>xlen*(330/361) && fast_r==0 && abs(mean(movement_r))>3 && ...
                size(severe_noise,2)<xlen*(80/361) || ...
                (noise_flag1==1 && movement_b_count_down1>xlen*(330/361) &&...
                fast_b==0 && abs(mean(movement_b))>3 && abs(mean(movement_r))>1)
            %如果没有上面往下动得那么剧烈，但是往下动的点数多，也加Ggc往下驱赶蛇。
            for i=1:xlen  % 或许可以更精细，不同的点，根据运动状态，不同地去加Ggc，不过后续没时间试了。
                if (~ismember(i,exists_noise) && I(yi_red(i),i)<10) || (ismember(i,exists_noise))
                    %如果没有噪声且LI附近灰度小，才加Ggc。主要是觉得LI老是误动作，蓝线倒是不至于。
                    if movement_r(i)>5, P(i,2)=P(i,2)+3;disp('red pushed down manually'),disp(fr),end
                    Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)=Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)+volcano*max(0,abs(mean(movement_r))-1)/8;
                    Eext_red(fix(P(i,2)),i)=Eext_red(fix(P(i,2)),i)+volcano*max(0,abs(mean(movement_r))-1)/4;
                end
            end
            disp('red moves down'),disp(abs(mean(movement_r))),disp(fr)
            r=1;
        end
    end
    % 类似于上面思路，处理MA。阈值会不太一样。
    if movement_b_count_down>xlen*(180/361)  && movement_b_count_down2>xlen*(50/361)%
        Q(:,2)=Q(:,2)+abs(mean(movement_b));%
        for i=1:xlen
            Eext_blue(max(1,fix(Q(i,2))-4):fix(Q(i,2)),i)=Eext_blue(max(1,fix(Q(i,2))-4):fix(Q(i,2)),i)+volcano;
            Eext_blue(fix(Q(i,2)),i)=Eext_blue(fix(Q(i,2)),i)+volcano;
        end
        disp('blue moves down fast'),disp(fr)
        fl_b=1;fast_b=1;
    else
        if movement_b_count_down1>xlen*(330/361) && fast_b==0 && abs(mean(movement_b))>3 || (r==1 && abs(mean(movement_b))>1)
            %MA没有LI那么容易往下降很多，其实并不是没降，而是看不出来。所以，这儿和LI的判定标准不太一样。两个条件之一：
            %①平均下降了3个；②平均下降了1个，但是LI的在下降“red moves down”即r==1
            for i=1:xlen
                if movement_b(i)>6 && grayscale_bmax_down_mean(i)>200, Q(i,2)=Q(i,2)+3;disp('blue pushed down manually'),disp(fr),end%
                %类似地，找到加Ggc的条件。
                Eext_blue(fix(Q(i,2))-4:fix(Q(i,2)),i)=Eext_blue(fix(Q(i,2))-4:fix(Q(i,2)),i)+volcano*max(0,abs(mean(movement_b))-1)/8;
                Eext_blue(fix(Q(i,2)),i)=Eext_blue(fix(Q(i,2)),i)+volcano*max(0,abs(mean(movement_b))-1)/4;
            end
            disp('blue moves down'),disp(abs(mean(movement_b))),disp(fr)
        end
    end
    % 3.4.8 根据整体运动趋势，设定蛇算法的行进步长。
    % 先是LI
    devi_r=abs(mean(movement_r));
    if fl_r==1%
        Gamma_red=Options.Gamma*min(10,max(2,devi_r));%。
    else
        Gamma_red=Options.Gamma;
    end
    if (noise_flag1==1 && fl_r==0) || (noise_flag2==1 && fl_r==0)%noise_flag1是对整个序列而言的。
        Gamma_red=Options.Gamma/2;
    end
    % MA线。
    devi_b=abs(mean(movement_b));
    if fl_b==1%
        Gamma_blue=Options.Gamma*min(10,max(2,devi_b));
    else
        Gamma_blue=Options.Gamma;
    end
    if lightflag==1, Gamma_blue=Options.Gamma/4;end%
    if (blue_suspect==1 && lightflag==0 && fl_b==0),Gamma_blue=Options.Gamma/3;end%
    if (fuzzy_flag1==1 || lightflag1==1 && blue_suspect==0 && lightflag==0 && fl_r==0)%
        Gamma_blue=Gamma_blue/2;
    end
    % 4 根据修正后的能量矩阵Eext，开始演化蛇的轮廓。%
    % 4.1 计算图像能量导数、以及LI、MA两条蛇的外部力（文章中4式的Gext的导数）。
    f0=1;%计算力的那个重力因子。
    Fx_red=ImageDerivatives2D(Eext_red,Options.Sigma2,'x');  % LI蛇的图像力，文章中的5式，这里力就是图像导数乘以一个因子。
    Fext_red(:,:,1)=-f0*Fx_red*2*Options.Sigma2^2;
    Fext_red(:,:,2)=0;
    Fx_blue=ImageDerivatives2D(Eext_blue,Options.Sigma2,'x');  % MA蛇的图像力
    Fext_blue(:,:,1)=-f0*Fx_blue*2*Options.Sigma2^2;
    Fext_blue(:,:,2)=0;
    if(Options.Verbose)
        [x,y]=ndgrid(1:10:size(Fext_red,1),1:10:size(Fext_red,2));
        figure,imshow(Eext_red,[]),title('Eextern after volcano (red)')
        figure,imshow(I,[]); hold on, quiver(y,x,Fext_red(1:10:end,1:10:end,2),Fext_red(1:10:end,1:10:end,1));        
        title('The external force field after volcano (red)')
        [x,y]=ndgrid(1:10:size(Fext_blue,1),1:10:size(Fext_blue,2));
        figure,imshow(Eext_blue,[]),title('Eextern after volcano (blue)')
        figure,imshow(I,[]); hold on, quiver(y,x,Fext_blue(1:10:end,1:10:end,2),Fext_blue(1:10:end,1:10:end,1));        
        title('The external force field after volcano (blue)')
    end
    % 4.2 计算蛇的内部力（文章中4式的Gint的导数）
    S=SnakeInternalForceMatrix2D(xlen,Options.Alpha,Options.Beta);  % 蛇受到的内力
    for i=1:Options.Iterations
    % 4.3 根据受力，推动蛇的运动。
        P=SnakeMoveIteration2D(S,P,Fext_red,Gamma_red,Options.Kappa,Options.Delta,eliminate_x_force);
        Q=SnakeMoveIteration2D(S,Q,Fext_blue,Gamma_blue,Options.Kappa,Options.Delta,eliminate_x_force);
    end
end