function Estereo=stereo(ylen,xlen,P,Q)
k1=10;
k2=10;
Estereo=zeros(ylen,xlen);
for i=1:xlen%wid是列数，wid=10表示该图中的第10列。
    for j=1:ylen%j是这一列不同位置处的y坐标
        Estereo(j,i)=k1*(j-P(i,2))^2+k2*(j-Q(i,2))^2;%P(i,2)是第i列红线所在位置的y坐标
    end
end
