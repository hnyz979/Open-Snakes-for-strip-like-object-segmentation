function J=ImageDerivatives2D(I,sigma,type)
% ���ڸ�˹�ֲ������άͼ�����x��y�����һ�ס����׵�����
% ���룺
%   I       ����ͼ��2D��
%   sigma	����ͼ����������
%   type	���֣��ֱ�Ϊ'x', 'y', 'xx', 'xy', 'yy'�������㵼���ķ���
% �����
%   J       ͼ����

% 1 ���ɶ�ͼ�����õ��ĸ�˹����ˡ�
[x,y]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
% ���ɵľ����DGauss��һ�����������ڸ��ֱ�Ե��ȡ���ӵ����á�
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
% 2 ͨ����ԭͼ���˹�˾����imfilter�������ķ��������ͼ��Ҷȵĵ�����
J = imfilter(I,DGauss,'conv','symmetric');