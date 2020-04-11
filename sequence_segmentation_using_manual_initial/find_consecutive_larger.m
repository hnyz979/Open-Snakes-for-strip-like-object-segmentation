function [biggest_count, longest_pos]=find_consecutive_larger(x,a)
%本函数找大于等于a的数连续最多出现了几次，以及最多连续的那些a值，是在x数组中的第几个。x是输入的矩阵，a是要找的数。
%比如说，x=[1 0 0 1 1 1 1 1 1 1 0 0 0 1 0 1 1];要找1出现的次数，应该返回7次、从第4个开始。
len=length(x);
count=[];c=0;
pos=[];
is_reached=zeros(1,len);
for i=1:len
    if is_reached(i)==0
        if x(i)>=a%搜索到第i个，a第一次出现
            pos=[pos i];%记录一下这次“连续的a的起始位置”
            c=1;
            for j=i+1:len%从第i个往后搜索
                if x(j)>=a%如果还等于a，那就把c+1
                  c=c+1;
                  is_reached(j)=1;
                else%否则，先结算c，然后把c清零。
                  count=[count c];
                  c=0;
                  is_reached(j)=1;%直接从第j+1个搜索就可以了，因为已经确定第i个，第i+1个，...第j-1个都是a了，而第j个不是。
                  %但问题是，直接写i=j是不行的，因为在外面的for循环里，i又给变成上次的i+1了。
                  break
                end
            end
        end
        is_reached(i)=1;
    end
end
count=[count c];%别忘了把最后一个加上去。否则如果是x1=[1 0 0 0 0 1 0 0 1 0 1 0 1 1 1 1 1 1 1];，就没法得到7，而得到1。
biggest_count=max(count);
if ~isempty(pos)
    longest_pos=pos(count==biggest_count);
else
    longest_pos=0;
end