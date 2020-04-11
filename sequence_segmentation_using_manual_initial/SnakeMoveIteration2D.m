function P=SnakeMoveIteration2D(B,P,Fext,gamma,kappa,delta,eliminate_x_force)
% ���ⲿ�����ڲ�������������(��)���˶���
% ���룺
%   B       ��������
%   P       �����㣬2*N�ľ���
%   Fext    �������������е�5ʽ��
%   gamma	����������Щ�������������ϵ�ʱ�䣬ʱ��Խ�����˶���Խ��
%   kappa	�ⲿ����Ȩ�أ�Խ��������������Ӱ��Խ��
%   delta	��������Ȩ��
% �����
%   P       ��������㣬2*N�ľ���

% �ѳ����߽�������㣬Ū�ر߽紦
P(:,1)=min(max(P(:,1),1),size(Fext,2));
P(:,2)=min(max(P(:,2),1),size(Fext,1));

% 1 �������ϵ�������
Fext1=zeros(size(P,1),1);
Fext1(:,1)=kappa*interp2(Fext(:,:,1),P(:,1),P(:,2));%��ԭ��������������������y��������ֵ��ÿ����������
if (eliminate_x_force),Fext1(:,2)=kappa*interp2(Fext(:,:,2),P(:,1),P(:,2));end
% ���Ƶأ����eliminate_x_force=1��������������x�������Ͱ���Ҳ��ֵ��ÿ����������
% ���ǣ�Ϊ�˼򻯼��㣬����ֻ��������y�����˶����������Fext1(:,2)��x�����ܵ���������Ϊ��ʼֵ0��
Fext1(isnan(Fext1))=0;  % ��nan�����0��

% 2 ���������������������������������һ������delta����
N=GetContourNormals2D(P);
Fext2=delta*N;
Fext2(isnan(Fext2))=0;

% 3 ��������������������ʱ�䲽�����õ����롣
% �и��ط��㷴�ˣ�P(:,1)��ͼ���е�x���꣬���ǣ�����������У�Fext1(i,1)�ǵ�i��y�����ܵ�����
% ��Fext2(i,1)�ǵ�i��y����������������Լӵ�ʱ�򣬻���Ť��һ�¡�����
if (eliminate_x_force), ssx = P(:,1) + gamma*Fext1(:,2) + gamma*Fext2(:,2); end
ssy = P(:,2) + gamma*Fext1(:,1) + gamma*Fext2(:,1); % �����y����������

% 4 ����߽�����
len=size(P,1);
expand_len=10;  % �߽�Ҫ��ܼ����㣿
% ����x����
if (eliminate_x_force)%eliminate_x_force=1��ʾ����x��������
    ssx_expanded(expand_len+1:len+expand_len)=ssx;%����P������Ŀ��Ӧ���Ǽ���ͼ����֮�󣬶�����֮ǰ��
    ssx_expanded(expand_len)=ssx(1)-1;
%     for i=1:expand_len-1,ssx_expanded(i)=ssx(expand_len-i);end%���淴��߽�����
    for i=1:expand_len-1,ssx_expanded(i)=ssx(expand_len);end%��ȱ߽�����
    ssx_expanded(len+expand_len+1)=ssx(len)+1;
%     for i=len+expand_len+2:len+expand_len*2,ssx_expanded(i)=ssx(len-(i-(len+expand_len+1)));end%���淴��߽�����
    for i=len+expand_len+2:len+expand_len*2,ssx_expanded(i)=ssx(len);end%��ȱ߽�����
    P_expand(:,1)=B*ssx_expanded'; 
    P(:,1)= P_expand(expand_len+1:len+expand_len,1);
end
% Ȼ����y����
ssy_expanded(expand_len+1:len+expand_len)=ssy;%����߽������������ܵ�ͼ���Ե����ô�죩��
% ssy_expanded(expand_len)=ssy(1)-1;%���淴��߽������ã�����ע�͵�
% for i=1:expand_len-1,ssy_expanded(i)=ssy(expand_len-i);end%���淴��߽������ã�����ע�͵�
% ssy_expanded(len+expand_len+1)=ssy(len)+1;%���淴��߽������ã�����ע�͵�
% for i=len+expand_len+2:len+expand_len*2,ssy_expanded(i)=ssy(len-(i-(len+expand_len+1)));end%���淴��߽������ã�����ע�͵�
% for i=1:expand_len,ssy_expanded(i)=ssy(expand_len+1);end%��ȱ߽������ã�����ע�͵�
% for i=len+expand_len+1:len+expand_len*2,ssy_expanded(i)=ssy(len);end%��ȱ߽������ã�����ע�͵�
ssy_dot_start=zeros(3,1);%б����ȱ߽������ã�����ע�͵�������ǳ�ʼ���������б�ʣ������㡣
for i=1:3, ssy_dot_start(i)=ssy(i+1)-ssy(i);end%б����ȱ߽������ã�����ע�͵��������ssy�е�ǰ�����㣨���߽磩��б�ʡ�
mean_ssy_dot_start=mean(ssy_dot_start);%б����ȱ߽������ã�����ע�͵���ǰ������б��ƽ��ֵ��
for i=expand_len:-1:1, ssy_expanded(i)=ssy_expanded(i+1)-mean_ssy_dot_start/(expand_len-i+1);end%б����ȱ߽������ã�����ע�͵���
% ��չ�ı߽��е����Ĵ�С������ԽԶ������iԽС��������ԽС��
ssy_dot_end=zeros(3,1);%б����ȱ߽������ã�����ע�͵���������ǰ�棬���ڴ�����������б�ʡ�
for i=1:3, ssy_dot_end(i)=ssy(len-3+i)-ssy(len-4+i);end%б����ȱ߽������ã�����ע�͵���
mean_ssy_dot_end=mean(ssy_dot_end);%б����ȱ߽������ã�����ע�͵���
for i=len+expand_len+1:len+expand_len*2, ssy_expanded(i)=ssy_expanded(i-1)+mean_ssy_dot_end/(i-len-expand_len);end%б����ȱ߽������ã�����ע�͵���
% ���Ƶأ����ı߽磬Ҳ������ԽԶ������iԽ�󣩣�����ԽС��

% 5 ��������B��Ӱ�죬��������µ�����λ�á�
P_expand(:,2)=B*ssy_expanded'; 
P(:,2)= P_expand(expand_len+1:len+expand_len,2);  % �Ѳ��ı߽�ɾ������
% ����߽��ϵĵ㣨ǰ5���ͺ�5���������ְ����Ǹ�ƽ��б���ߣ��������ǡ�Ʈ������
% ʵ������£������������Ļ���ͼ���Ե��һЩ��ǿ�ĻҶ��ݶȣ�����߽�����ߣ��ر𷳡����Լ��˸����취��
P_dot_start=zeros(5,1);
for i=1:5, P_dot_start(i)=P(i+5,2)-P(i+4,2); end%����B�󣬵�6~10����ĵ���
mean_P_dot_start=mean(P_dot_start);
for i=5:-1:1, P(i,2)=P(i+1,2)-mean_P_dot_start;end%
P_dot_end=zeros(5,1);
for i=1:5, P_dot_end(i)=P(len-10+i,2)-P(len-11+i,2); end%����B�󣬵�����10~6����ĵ���
mean_P_dot_end=mean(P_dot_end);
for i=1:5, P(len-5+i,2)=P(len-6+i,2)+mean_P_dot_end;end%
% ��������߽��ˣ��������ص��߽��ϡ�
P(:,1)=min(max(P(:,1),1),size(Fext,2));
P(:,2)=min(max(P(:,2),1),size(Fext,1));