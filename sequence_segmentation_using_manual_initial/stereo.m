function Estereo=stereo(ylen,xlen,P,Q)
k1=10;
k2=10;
Estereo=zeros(ylen,xlen);
for i=1:xlen%wid��������wid=10��ʾ��ͼ�еĵ�10�С�
    for j=1:ylen%j����һ�в�ͬλ�ô���y����
        Estereo(j,i)=k1*(j-P(i,2))^2+k2*(j-Q(i,2))^2;%P(i,2)�ǵ�i�к�������λ�õ�y����
    end
end
