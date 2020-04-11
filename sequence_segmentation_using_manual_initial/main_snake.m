clear
close all  
clc
% -----------------------0��������֡���ֶ��ָ�����-----------------------
root_dir = ['B:\02-����Ѫ��׷�����ݣ���ʦ�ܵĲ������£�\matlab����\�ϴ�github�İ汾�����ڱ���25��\' ...
    'Open-Snakes-for-strip-like-object-segmentation\data\'];
program_dir = ['B:\02-����Ѫ��׷�����ݣ���ʦ�ܵĲ������£�\matlab����\�ϴ�github�İ汾�����ڱ���25��\' ...
    'Open-Snakes-for-strip-like-object-segmentation\sequence_segmentation_using_manual_initial'];
% ���Ϸֱ��Ǵ�����ݵ�Ŀ¼�ʹ�ų����Ŀ¼��
sequence_name = '1';
load([root_dir sequence_name '\LI.mat']);  % load����֮��LI�ĳ�ʼֵ����red_manu
load([root_dir sequence_name '\MA.mat']);  % load����֮��MA�ĳ�ʼֵ����blue_manu
I=imread([root_dir sequence_name '\000001'],'bmp');
rmax_index=fix(red_manu);%������Ρ���
bmax_index=fix(blue_manu);
xlen=length(red_manu);
if(size(I,3)==3), I=rgb2gray(I); end%��I��ɻҶ�ͼ��
figure, imshow(I), hold on;
plot(1:xlen,rmax_index,'r.',1:xlen,bmax_index,'b.');
title('��1��ͼ')
hold off
Red=[1:xlen; rmax_index']';%LI��ʼ���Ƶ㼯��Red�������Ǹ�2*xlen�ľ���
Blue=[1:xlen; bmax_index']';%MA��ʼ���Ƶ㼯��Blue�������Ǹ�2*xlen�ľ���
% -----------------------1���ж���֡�Ҷȣ���������б�ʵ��������Ϊ����֡����Ļ�׼��-----------------------
grayscale_rmax_up_mean=zeros(1,xlen);
grayscale_rmax_down_mean=zeros(1,xlen);
grayscale_bmax_up_mean=zeros(1,xlen);
grayscale_bmax_down_mean=zeros(1,xlen);grayscale_bmax_mid_std=zeros(1,xlen);
inside_IMC_is_bright=zeros(1,xlen);%�����Ƿ���MA�Ϸ����·������ģ��������ܲ�̫�����ġ�
inside_IMC_is_dark=zeros(1,xlen);%�����Ƿ���LI�·��ܰ��ģ��������ܲ�̫�����ġ�
IMW=bmax_index-rmax_index;%��ʼ������IM��ȡ��Ǹ�1*361�������������������������������ֵ��
% 1.1 ����ͨ���Ҷ�ֵ��������������
% 1.1.1 �۲�LI�Ϸ���MA�·��ռ��ڵĻҶ�ֵ��ͨ���������������������LI�Ϸ�����MA�·����Ĵ�С��
neighbor_out_r=min(max(4,fix(IMW*0.7)),9);%����IMWֵѡȡ�������С��ʼֵ��
neighbor_out_b=min(max(4,fix(IMW*0.7)),9);
for i=1:xlen
    r_small=I(rmax_index(i)-1:rmax_index(i)+1,max(1,i-2):min(i+2,xlen));
    b_small=I(bmax_index(i)-1:bmax_index(i)+1,max(1,i-2):min(i+2,xlen));
    for j=rmax_index(i)-4:-1:rmax_index(i)-9%�ڳ�ʼLI�Ϸ����������Ϸ�4����������9���㡣
        if i==104
            pause(0.1)  % �����������������ÿһ����������
        end
        if I(j,i)>=mean2(r_small)+15%�����һ����ĻҶȣ�����LI����ƽ���Ҷȴ���15������Ϊ��һ��ܿ����������ˡ�
            neighbor_out_r(i)=rmax_index(i)-j-1;%��ô����¼��������λ�á��������i��LI���������ȡrmax_index(i)-j��
            break
        end
    end
    for k=bmax_index(i)+4:bmax_index(i)+9%��ʼMA�·����������·�4�����������·�9���㡣
        if I(k,i)<=mean2(b_small)-15%�����һ����ĻҶȣ�����MA����ƽ���Ҷ�С��15������Ϊ��һ��ܿ����ǲ����ˡ����磬MA�·��ġ��ϰ��ºڵĵط���
            neighbor_out_b(i)=k-bmax_index(i)-1;
            break
        end
    end
end
% 1.1.2 ȷ��������LI�Ϸ���MA�·�����������LI�·���MA�Ϸ�����ƽ���Ҷȣ�����ˮƽ����
neighbor_in=min(max(4,fix(IMW/3)),6);%���������IM��խ�Ļ�������Ȼȡ4������ͼӿ�һЩ
for i=1:xlen
    r_up=I(rmax_index(i)-neighbor_out_r(i)-1:rmax_index(i)-1,max(1,i-2):min(i+2,xlen));
    r_down=I(rmax_index(i)+1:rmax_index(i)+neighbor_in(i)+1,max(1,i-2):min(i+2,xlen));%�������棬ֻ��������Ǽ���1
    b_up=I(bmax_index(i)-neighbor_in(i)-1:bmax_index(i)-1,max(1,i-2):min(i+2,xlen));%ͬr_up
    b_down=I(bmax_index(i)+1:bmax_index(i)+neighbor_out_b(i)+1,max(1,i-2):min(i+2,xlen));%ͬr_down
    b_mid=I(bmax_index(i)-neighbor_in(i):bmax_index(i)+neighbor_in(i),max(1,i-2):min(i+2,xlen));
    grayscale_rmax_up_mean(i)=mean2(r_up);%���ʼLI�Ϸ������ڵ�ƽ���Ҷȡ�
    % ���磬��ʼLI��һ��Ϊ(3,202)�������ͼ���д�(1,198)��ʼ��(5,202)��������25�����ƽ���Ҷȡ�
    grayscale_rmax_down_mean(i)=mean2(r_down);%���ʼLI�·������ڵ�ƽ���Ҷȡ�
    % ���磬��ʼLI��һ��Ϊ(3,202)�������ͼ���д�(1,202)��ʼ��(5,206)��������25�����ƽ���Ҷȡ�
    grayscale_bmax_up_mean(i)=mean2(b_up);%���ʼMA�Ϸ������ڵ�ƽ���Ҷȡ�
    grayscale_bmax_down_mean(i)=mean2(b_down);%���ʼMA�·������ڵ�ƽ���Ҷȡ�
    grayscale_bmax_mid_std(i)=std2(b_mid);%���ʼMA�·������ڵĻҶȱ�׼�
    grayscale_inside=I(rmax_index(i)+1:bmax_index(i),i);
    if sum(grayscale_inside>grayscale_bmax_down_mean(i))>=10
        % �����һ�У�����Ĥ֮����10�����ϵĵ�ĻҶȴ����·�ƽ���Ҷȣ���������ǣ���Ҫ���׼�����ֵ��ʹ��MA�����ƶ��ơ�
        % �˴���Ӧ����P4������exceptions���иĶ���
        inside_IMC_is_bright(i)=1;
    end
    if sum(grayscale_inside<=10)>=10
        % �����һ�У�����Ĥ֮����10�����ϵĵ�ĻҶ�С��10����������ǣ���Ҫ���׼�����ֵ��ʹ��LI�����ƶ���
        inside_IMC_is_dark(i)=1;
    end
end
exists_noise_temp=find(grayscale_rmax_up_mean>=10);%LI�Ϸ��ϴ�����Ҷ�ֵ����10��������Ϊ��������
temp=0;
for k=1:size(exists_noise_temp,2)
    % �����ѭ�����ҳ�exists_noise_temp�еĵ��10-���������е㡣���磬���exists_noise_temp����һ��67����ô��57~77����Ϊ������������
    qu=ismember(max(1,exists_noise_temp(k)-2):min(exists_noise_temp(k)+2,xlen),exists_noise_temp);
    % �жϵ�k�����������ߵ����Ƿ�Ҳ��exists_noise_temp�Ҳ����˵���ж��Ƿ�����5���㶼��������
    if min(abs(qu))==1
        % ���qu���5��������1�Ļ�������Ϊ��һ�������󣬲��ҳ�10�����ڵĵ㡣
        sus=max(1,exists_noise_temp(k)-10):min(exists_noise_temp(k)+10,xlen);
    else
        sus=exists_noise_temp(k);
    end
    temp=[temp sus];
end
temp(temp==0)=[];%ɾ����һ��0
exists_noise=unique(temp);
medium_noise=find(grayscale_rmax_up_mean>=20 & grayscale_rmax_up_mean<30);%LI�Ϸ��Ҷ�ֵ����20��С��30����Ϊ�����е����ء�
severe_noise=find(grayscale_rmax_up_mean>=30);%LI�Ϸ��Ҷ�ֵ����30����Ϊ�����Ƚ����ء�
%�����������س̶ȣ������ߺ����лҶȴﵽ���ٲż�Grayscale constraint energy��
exists_fuzzy_temp=find(grayscale_bmax_mid_std<=20);%
temp=0;
for k=1:size(exists_fuzzy_temp,2)
    %������ǰ�棬����MA��������
    qu=ismember(max(1,exists_fuzzy_temp(k)-2):min(exists_fuzzy_temp(k)+2,xlen),exists_fuzzy_temp);
    if min(abs(qu))==1
        sus=max(1,exists_fuzzy_temp(k)-10):min(exists_fuzzy_temp(k)+10,xlen);
    else
        sus=exists_fuzzy_temp(k);
    end
    temp=[temp sus];
end
temp(temp==0)=[];%ɾ����һ��0
exists_fuzzy=unique(temp);
medium_fuzzy=find(grayscale_bmax_mid_std<=15 & grayscale_bmax_mid_std>10);
severe_fuzzy=find(grayscale_bmax_mid_std<=10);
lightflag=0;lightflag1=0;noise_flag0=0;noise_flag1=0;noise_flag2=0;fuzzy_flag1=0;
if size(exists_noise_temp,2)>15%�õ���exists_noise_temp����û��չ֮ǰ��
    noise_flag0=1;
    disp('������LI�Ϸ�������')
end
if size(find(grayscale_bmax_up_mean>=150),2)>200%�����200�����ϵĵ㣬����MA�Ϸ��Ҷ�ֵ����150����1055L1�õġ�
    lightflag=1;%snake2D�����У�MA���Ϲİ������Ƶ�ʱ�����һ�䣬��if lightflag==1, vol_times=vol_times*4;end��
    disp('������MA�Ϸ�����')%�Ѳ��ԣ�ֻ��1055L1MA�Ϸ�����
end
if size(find(grayscale_bmax_down_mean>=200),2)>200%�����200�����ϵĵ㣬����MA�·��Ҷ�ֵ����200��
    lightflag1=1;%
    disp('������MA�·�����')%�Ѳ��ԣ�39R��47L/R��48R��49L/R��55���ĸ���59R��60L/R��62L/R��63L/R��64L��66L1�����������
end
if size(exists_noise_temp,2)>200 && size(medium_noise,2)>=30%50R
    noise_flag1=1;
    disp('������LI�Ϸ���������')
end
if size(severe_noise,2)>size(medium_noise,2) && size(severe_noise,2)>=100%
    noise_flag2=1;
    disp('������LI�Ϸ��������أ�������������һ����������')
end
if size(exists_fuzzy_temp,2)>=30%�õ���exists_fuzzy_temp����û��չ֮ǰ��
    fuzzy_flag1=1;
    disp('������MA��30��������ģ��')
end
% 1.2 ����ʼ�����еĵ����������Ӧ����Derivative constraint energy��
red_dot1=Red(:,2);blue_dot1=Blue(:,2);%��ʼ������ʵ��Ū��0Ҳ���Եġ�
for i=1:xlen-1, red_dot1(i)=Red(i+1,2)-Red(i,2);blue_dot1(i)=Blue(i+1,2)-Blue(i,2);end%���㵼��
red_dot1(xlen)=red_dot1(xlen-1);blue_dot1(xlen)=blue_dot1(xlen-1);
% -----------------------2�������2֡ͼ�������㷨�ݻ�������������H-inf�Ĳ�����-----------------------
% 2.1 �ݻ��õ���2֡ͼ��LI��MA��
I2=imread([root_dir sequence_name '\000002'],'bmp');
if(size(I2,3)==3), I2=rgb2gray(I2); end
sigma_filt=1;
I2 = gaussian_filter(I2,sigma_filt);
I2(:,1)=I2(:,2);I2(:,xlen)=I2(:,xlen-1);%�˲���ʱ�򣬵�һ�к����һ�л�����⣬�ر�ڡ����ڸĵ���
Options=struct;
Options.Iterations=20;
Options.Verbose=false;
i=1;
fast_r=0;fast_b=0;%LI��MA�Ƿ��˶��ǳ����flag���ȳ�ʼ��Ϊ0����ĳ�������У�һ����������������1���Ͳ����ٷ���Ϊ0�ˡ�
[Red1,Blue1,J,fast_r,fast_b]=Snake2D(i,I2,Red,Blue,Options,neighbor_in,neighbor_out_r,neighbor_out_b, ...
    grayscale_rmax_up_mean,grayscale_rmax_down_mean,grayscale_bmax_up_mean,grayscale_bmax_down_mean, ...
    inside_IMC_is_bright,inside_IMC_is_dark,severe_noise,medium_noise,exists_noise, ...
    severe_fuzzy,medium_fuzzy,exists_fuzzy,IMW,red_dot1,blue_dot1,fast_r,fast_b,...
    lightflag,lightflag1,noise_flag0,noise_flag1,noise_flag2,fuzzy_flag1);%����ȷ����һ��ͼ������
% 2.2 ���״̬�ռ䷽�̡�
Options.Verbose=false;
q=1;
r=1;%��������Э����RԽС�������㣩�� �����������ߵĽ������Ȩ��Խ��Խ��
matrixSize=size(Red,1);
A = eye(matrixSize,matrixSize);
H = eye(matrixSize,matrixSize);
Q = q*eye(matrixSize,matrixSize);
R = r*eye(matrixSize,matrixSize);
% ��karmen��H inf��ʱ�򣬰���˵��ֻȡRed(:,2)��Blue(:,2)�Ϳ����ˣ���Ϊ���ǵĵ�һά�ǲ���ģ�����1��len��
Red_new=Red; Blue_new=Blue;
u=0;G=0;
L=eye(matrixSize,matrixSize);S=eye(matrixSize,matrixSize);I=eye(matrixSize,matrixSize);
theta=0.1;
[Red_new(:,2),K_red] = Hinf(Red(:,2),Red1(:,2),u,A,G,H,L,R,Q,S,I,theta);
[Blue_new(:,2),K_blue] = Hinf(Blue(:,2),Blue1(:,2),u,A,G,H,L,R,Q,S,I,theta);
%ʵ������������=0.5����ô��K_red��K_blue���ǵ�λ������ζ����һ�ε�LI��MA��ֻ���ߵĽ���йأ������ϴεĽ���޹ء���02�Ǳ����11.102ʽ��
figure, imshow(I2,[]), hold on; 
plot(Red(:,1),Red(:,2),'m.');
plot(Red_new(:,1),Red_new(:,2),'r.');
plot(Blue(:,1),Blue(:,2),'c.');
plot(Blue_new(:,1),Blue_new(:,2),'b.');hold off;
title('��2��ͼ')
% -----------------------3����һ�������ͼ���������룬��2֡�ͺ���֡���Է���һ�������Ȳ���ʱ���޸��ˡ�-----------------------
srcDir=[root_dir sequence_name]; %ѡ����Щͼ���Ŀ¼
cd(srcDir);%�����Ŀ¼
allnames=struct2cell(dir('*.bmp')); %��Ŀ¼�µ�����bmpͼ��
[~,len]=size(allnames);%len��ͼ������
str1='��';str2='��ͼ';%��ͼʱд��Ŀ�õ�
%��ʼ��һ��3ά��Red_3D����һά�͵ڶ�ά��ǰ����ͬ������ά���ʾ���ǵ�k��ͼ�Ľ����Blue��֮��ͬ��
Red_3D=zeros(size(Red,1),size(Red,2),len);
Blue_3D=zeros(size(Blue,1),size(Blue,2),len);
Red_3D(:,:,1)=[Red(:,1),Red(:,2)];%000001ͼ�У�ԭʼ�ֶ����Ľ��
Red_3D(:,:,2)=[Red_new(:,1),Red_new(:,2)];%000002ͼ�У���+�������Ľ����
Blue_3D(:,:,1)=[Blue(:,1),Blue(:,2)];
Blue_3D(:,:,2)=[Blue_new(:,1),Blue_new(:,2)];
red_dot2=Red_new(:,2);blue_dot2=Blue_new(:,2);%��ʼ������ʵ��Ū��0Ҳ���Եġ�
for i=1:xlen-1, red_dot2(i)=Red_new(i+1,2)-Red_new(i,2);blue_dot2(i)=Blue_new(i+1,2)-Blue_new(i,2);end%���㵼��
red_dot2(xlen)=red_dot2(xlen-1);blue_dot2(xlen)=blue_dot2(xlen-1);
%��¥��ÿһ�㶼��һ������ʱ�̣���һ��k��
%�ӵ�3�ſ�ʼ�������Ļ���Red_3D(:,:,i)�Ͷ�Ӧ�ڵ�i��ͼ��LI�ĺ������ꡣ
figure,
for i=3:min(100,len)
    name=allnames{1,i};
    I=imread(name);%��ʱI�ǵ�i��ͼ��������ľ���
    if(size(I,3)==3), I=rgb2gray(I); end
    cd(program_dir)%�ո��Ǵ��ͼƬ���ļ��У�����Ҫ�ص���ų�����ļ��С�
    I = gaussian_filter(I,sigma_filt);
    I(:,1)=I(:,2);I(:,xlen)=I(:,xlen-1);%�˲���ʱ�򣬵�һ�к����һ�л�����⣬�ر�ڡ����ڸĵ���
    if i==30 %��Ҫ�����һ��ͼ������i���ڼ���Ȼ��Ӷϵ㡣�ϵ�󵥲����ԣ�����Snake2D���򼴿ɹ۲��ʱ��Eext�Ⱦ���
        pause(0.1);
    end
    [Red_3D(:,:,i),Blue_3D(:,:,i),J,fast_r,fast_b]=Snake2D(i,I,Red_3D(:,:,i-1),Blue_3D(:,:,i-1),Options, ... 
        neighbor_in,neighbor_out_r,neighbor_out_b,grayscale_rmax_up_mean,grayscale_rmax_down_mean,grayscale_bmax_up_mean, ...
        grayscale_bmax_down_mean,inside_IMC_is_bright,inside_IMC_is_dark,severe_noise,medium_noise,exists_noise, ...
        severe_fuzzy,medium_fuzzy,exists_fuzzy,IMW,red_dot2,blue_dot2,fast_r,fast_b,lightflag,lightflag1,...
        noise_flag0,noise_flag1,noise_flag2,fuzzy_flag1);
    Red_3D(:,2,i) = A * Red_3D(:,2,i-1) + G * u + A * K_red * (Red_3D(:,2,i) - H * Red_3D(:,2,i-1));
    Blue_3D(:,2,i) = A * Blue_3D(:,2,i-1) + G * u + A * K_blue * (Blue_3D(:,2,i) - H * Blue_3D(:,2,i-1));
    %����һ�У��Ⱥ��ұ�Red_3D(:,2,i-1)����һ��ͼ�Ľ����Red_3D(:,2,i)����һ��ͼ�Ĺ۲�ֵ���ߵĽ������K_red��11.89ʽ�е�Kֵ������Ϊ����ġ�
    %Ȼ��Ⱥ���ߵ�Red_3D(:,2,i)����һ��ͼ�Ľ��������H inf��Ľ������
    imshow(I,[])
    hold on; 
    plot(Red_3D(:,1,i),Red_3D(:,2,i),'r.');%Red_3D(:,1,i)��Ȼ��1:len��Щ��Ȼ����Red_3D(:,2,i)�ǵ�k��ͼ�������LI���ꡣ
    plot(Blue_3D(:,1,i),Blue_3D(:,2,i),'b.'); 
    str=[str1 num2str(i) str2];
    title(str)
    hold off
    pause(0.01);
    cd (srcDir)
end
% save B:\02-����Ѫ��׷�����ݣ���ʦ�ܵĲ������£�\matlab����\�ϴ�github�İ汾�����ڱ���25��\Red_1Dsnake_Hinf Red_3D
% save B:\02-����Ѫ��׷�����ݣ���ʦ�ܵĲ������£�\matlab����\�ϴ�github�İ汾�����ڱ���25��\Blue_1Dsnake_Hinf Blue_3D