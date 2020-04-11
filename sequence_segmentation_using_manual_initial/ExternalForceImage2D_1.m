function [Eextern, I_deri] = ExternalForceImage2D_1(I,Wline, Wedge, Wterm,Sigma,Verbose)
% ��ͼ���������
% 
% ���룺 
%  I        ԭͼ
%  Sigma	����ͼ���������� 
%  Wline    ������4ʽ������Ȩ��ֵ
%  Wedge	������4ʽ������Ȩ��ֵ
%  Wterm	������4ʽ������Ȩ��ֵ
% �����
%  Eextern	�����������������е�Gext

% 1 �����ͼ���x����y�ĸ��׵�����
Ix=ImageDerivatives2D(I,Sigma,'x');
Iy=ImageDerivatives2D(I,Sigma,'y');
Ixx=ImageDerivatives2D(I,Sigma,'xx');
Ixy=ImageDerivatives2D(I,Sigma,'xy');
Iyy=ImageDerivatives2D(I,Sigma,'yy');
threshhold=1e-6;
Ix(abs(Ix)<threshhold)=0;
% ��Ix(find(Ix<threshhold))=0;����ЩС��1e-6�Ķ�Ū��0��Ȼ��ͨ���Ǹ�Eterm(isnan(Eterm))=0��Eterm��������ӦԪ��Ū��0��
% ����Ļ���Eterm�����ڷǳ������ֵ����1e+11���Լ���ǳ�������
% ֮����ȡ1e-6������Ϊ������(Ix.^2 + Iy.^2).^(3/2)��ͼ���������������ֵ����1e-5���ϣ�������С�Ķ���������
Iy(abs(Iy)<threshhold)=0;
Ixx(abs(Ixx)<threshhold)=0;
Ixy(abs(Ixy)<threshhold)=0;
Iyy(abs(Iyy)<threshhold)=0;
% 2 ����Eline��4ʽ�е�G_l������ʵ����ԭͼ���������ûɶ�á�
Eline = I;
if(Verbose),  figure,imshow(Eline,[]),title('Eline'); end%�����л�Eline�Ķ�ά�ֲ�ͼ�����Ժ����ע�͵���
% 3 ����Eterm��4ʽ�е�G_t��
Eterm = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((Ix.^2 + Iy.^2).^(3/2));
Eterm(isnan(Eterm))=0;
Eterm_threshhold=1e+0;
Eterm(abs(Eterm)>Eterm_threshhold)=0;
%Eterm�治��̫��Ӧ�ð��Ǹ�Eterm_threshhold�������10���������1������ֹ������ɷ�ĸ��С�����ߵ��˶�����Ӱ�졣
if(Verbose), figure,imshow(Eterm,[]),title('Eterm'); end%�����л�Eterm�Ķ�ά�ֲ�ͼ�����Ժ����ע�͵���
% 4 ����Eedge��4ʽ�е�G_e��
[grady,gradx] = gradient(I);%����һ�ַ�����Ix��Iy������Ч����Ρ�
Eedge = sqrt(gradx.^2 + grady.^2); 
if(Verbose), figure,imshow(Eedge,[]),title('Eedge'); end%%�����л�Eedge�Ķ�ά�ֲ�ͼ�����Ժ����ע�͵���
I_deri=Ix;
% 5 �Ӻ͵õ�Eextern�������е�Gext�� = (Wline*Eline - Wedge*Eedge -Wterm * Eterm); 
Eextern= (Wline*Eline - Wedge*Eedge + Wterm * Eterm); 