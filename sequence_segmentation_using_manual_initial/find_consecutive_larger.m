function [biggest_count, longest_pos]=find_consecutive_larger(x,a)
%�������Ҵ��ڵ���a���������������˼��Σ��Լ������������Щaֵ������x�����еĵڼ�����x������ľ���a��Ҫ�ҵ�����
%����˵��x=[1 0 0 1 1 1 1 1 1 1 0 0 0 1 0 1 1];Ҫ��1���ֵĴ�����Ӧ�÷���7�Ρ��ӵ�4����ʼ��
len=length(x);
count=[];c=0;
pos=[];
is_reached=zeros(1,len);
for i=1:len
    if is_reached(i)==0
        if x(i)>=a%��������i����a��һ�γ���
            pos=[pos i];%��¼һ����Ρ�������a����ʼλ�á�
            c=1;
            for j=i+1:len%�ӵ�i����������
                if x(j)>=a%���������a���ǾͰ�c+1
                  c=c+1;
                  is_reached(j)=1;
                else%�����Ƚ���c��Ȼ���c���㡣
                  count=[count c];
                  c=0;
                  is_reached(j)=1;%ֱ�Ӵӵ�j+1�������Ϳ����ˣ���Ϊ�Ѿ�ȷ����i������i+1����...��j-1������a�ˣ�����j�����ǡ�
                  %�������ǣ�ֱ��дi=j�ǲ��еģ���Ϊ�������forѭ���i�ָ�����ϴε�i+1�ˡ�
                  break
                end
            end
        end
        is_reached(i)=1;
    end
end
count=[count c];%�����˰����һ������ȥ�����������x1=[1 0 0 0 0 1 0 0 1 0 1 0 1 1 1 1 1 1 1];����û���õ�7�����õ�1��
biggest_count=max(count);
if ~isempty(pos)
    longest_pos=pos(count==biggest_count);
else
    longest_pos=0;
end