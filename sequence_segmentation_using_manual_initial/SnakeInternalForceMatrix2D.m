function B=SnakeInternalForceMatrix2D(nPoints,alpha,beta)
% �����ߵ��������󣨱����ߵĹ⻬�ԣ�
% ���룺
%   nPoints     �����м��������㣿
%   alpha       ����4ʽ��alpha������һ�׵�����Ӧ������
%   beta        ����4ʽ��beta�����߶��׵�����Ӧ������
%   gamma       ʱ�䲽��
% �����
%   B           ��������

% ����5ʽ�У�����5���ϵ������ĸ����Ҫ������Ȼ���֣����е����ڶ��������˸�2��
b(1)=beta;
b(2)=-(alpha + 4*beta);
b(3)=(2*alpha + 6 *beta);
b(4)=b(2);
b(5)=b(1);

% ����M���󣨼�Everything You Always Wanted To Know About Snakes��20ҳ6.6ʽ�����Ǹ����£��м����ط�Ҳ�����⡣����
nPoints=nPoints+20;%��Ϊ�߲����������Ը�A��B�������㣬ǰ�������10��10�С�
% ע�⣬���ڱ��������ص��Ǹ�B���󣬾���������B���󡣵�ʱ����SnakeMoveIteration2D���ٰ��Ǹ�P����ǰ���������淴��ı߽�������
A=b(1)*circshift(eye(nPoints),2);%M���������Ͻǵ�����p�������������һб��p
A=A+b(2)*circshift(eye(nPoints),1);%���Ͻǵ�һ��q���͵����ڶ����µ���һб��q
A=A+b(3)*circshift(eye(nPoints),0);%�м�Խ����ϵ�r������rӦ����1+2*alpha + 6 *beta������û��1��һ����ټӡ�
A=A+b(4)*circshift(eye(nPoints),-1);
A=A+b(5)*circshift(eye(nPoints),-2);
A=A + eye(nPoints);
% ������ȫ����01�е����ɣ��Ͱ�ԭ����B=inv(gamma * A + eye(nPoints));����B=inv(A + gamma.* eye(nPoints));�е�gammaȥ���ˡ�
B=inv(A);