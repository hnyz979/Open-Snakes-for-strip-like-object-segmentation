function [repaired,red_alert]=repair_bumpupanddown(I,unrepaired,tol1,tol2,wid,mode,red_alert)
%本程序根据灰度和斜率，修正两个方向的鼓包。斜率阈值有两个，默认tol1=0.1，tol2=0.2。wid为修正宽度。
index=unrepaired;len=length(index);
index_diff=zeros(1,len-1);
y=fix(index);grayscale_std=zeros(1,len);grayscale_mean=zeros(1,len);
f=0;
for i=1:len-1
    if i==333
        pause(0.001)
    end
    index_diff(i)=index(i+1)-index(i);
    grayscale_mean(i)=mean2(I(y(i)-2:y(i)+2,max(1,i-2):min(len,i+2)));
    grayscale_std(i)=std2(I(y(i)-2:y(i)+2,max(1,i-2):min(len,i+2)));
end
grayscale_mean(len)=mean2(I(y(len)-2:y(len)+2,len-1:len));
grayscale_std(len)=std2(I(y(len)-2:y(len)+2,len-1:len));
count_gt0_2=length(find(index_diff>0.2));%有多少个大于0.2
count_lt_0_2=length(find(index_diff<-0.2));%有多少个小于-0.2
count_all0_2=count_gt0_2+count_lt_0_2;%二者之和，即多少个绝对值大于0.2的。
red_alert_default=[0 0 0 0];
if strcmp(mode,'red'), red_alert=[0 0 0 0];end%如果是红线的话，就初始化为0。
%如果虽然前面10个都往下斜，且整体并不是往下斜的，但不满足红线灰度条件，所以没有修正的话，那就不修正蓝的。
if strcmp(mode,'blue')%如果是蓝线，而且没有外来的输入，也弄成0
    if ~exist('red_alert','var')
        red_alert=red_alert_default;
    end
end
inclined=0;
if (count_all0_2>=80 && count_gt0_2>=0.75*count_all0_2) || (count_all0_2>=60 && count_gt0_2>=0.85*count_all0_2)
    inclined=1;%如果绝对值大于0.2的有100列，而且，其中有80%以上的都大于0.2，那就认为确实是往下斜的
end
if (count_all0_2>=80 && count_lt_0_2>=0.75*count_all0_2) || (count_all0_2>=60 && count_lt_0_2>=0.85*count_all0_2)
    inclined=-1;%如果绝对值大于0.2的有100列，而且，其中有80%以上的都小于-0.2，那就认为确实是往上斜的。69L改成80列。41R加了“或者大于60列，而其中85%以上都是一个方向的”
end
for i=11:len-11%中间那些列，注意index_diff只有i-1个元素。前10列和后10列单独考虑，见后。
    if i==257
        pause(0.001)
    end
    if (index_diff(i)<-0.001 && index_diff(i+1)>0.001 || ...
            (index_diff(i)<-0 && index_diff(i+1)>0 && index_diff(i-1)<-0.04 && index_diff(i+2)>0.04)) 
        flags_front1=(index_diff(max(1,i-9):i)<-tol1);%看看它前面的点，有多少个斜率是小于-tol1的
        flags_front2=(index_diff(max(1,i-9):i)<-tol2);%看看它前面的点，有多少个斜率是小于-tol2的
        flags_rear1=(index_diff(min(i+1,len-1):min(i+10,len-1))>tol1);%看看它后面的点，有多少个斜率是大于tol1的
        flags_rear2=(index_diff(min(i+1,len-1):min(i+10,len-1))>tol2);%看看它后面的点，有多少个斜率是大于tol2的
        if (sum(flags_front1)>=9&&sum(flags_rear1)>=9) && (sum(flags_front2)>=6&&sum(flags_rear2)>=6)%如果前后各有9个斜率绝对值大于0.1，而且有6个大于0.2的
            %如果满足以上条件，就认为有个鼓包了。
            %先找修正的范围
            last=min([i-1,len-i-1,wid]);%现在是第i个点，前面有i-1个点，后面有len-i-1个点（因为index_diff只有len-1个数）。设个限制，最多搜索wid个好了。
            for k=10:last%以第i个点为中心搜索，从它两边各10个点开始搜索。如果第k-9次搜索（即第i-k和i+1+k个点）满足：第i-k个点导数大于-0.1且大于第i-k+1个点（绝对值小于0.1，且上升的趋势没有前一次搜索的大），且第i+k个点导数小于0.1且小于第i+k-1个点的导数。
                if index_diff(i-k)>-tol1&&index_diff(i-k)>index_diff(i-k+1) ...%左边：第k次搜索是第i-k个点，如果它的导数绝对值小于0.1（负的，而且大于-0.1），且绝对值比上一次搜索（第i-k+1个点）要小（负的，那个值大）
                        && index_diff(i+k)<tol1&&index_diff(i+k)<index_diff(i+k-1)%右边：第i+k个点，如果导数绝对值小于0.1（正的，且小于0.1），且绝对值比上一次搜索（第i+k-1个点）要小（正的，那个值小）
                    %忽然想起来，一定要index_diff(i-k)>-tol1吗？如果下一个比她更小，岂不是更好？？？
                    forced_area=k;%i-forced_area+1:i+forced_area-1这个范围内，进行鼓包修正。（第i-forced_area个和第i+forced_area个都已经是正常的了啊）
                    break
                end
                forced_area=k;%如果搜索了wid个，还没有满足上述条件的，就弄这wid个了。
            end
            %然后开始修正
            repair_begin=index(i-forced_area);%用第i-forced_area个和第i+forced_area个点的y坐标为基准修正（这两个点都是正常的）。
            repair_end=index(i+forced_area);
            %170522加一段，如果(2*forced_area-1)>3*wid就不修正了。这种情况相当于鼓包太宽（90个以上啊）了。还真有这样的。
            if (2*forced_area-1)<=3*wid%170522
                derta_y=(repair_end-repair_begin)/(2*forced_area);%从第i-forced_area+1个修正到第i+forced_area-1个。那么，有2*forced_area-1个点要插值，但是，每两个点之间的距离，应该是除以2*forced_area而不是2*forced_area-1。
                for m=1:2*forced_area-1
                    index(i-forced_area+m)=index(i-forced_area)+m*derta_y;
                end
            end
        else%这个else修正顶部较平的鼓包。如果顶部较平，就从顶部向两个方向各找10个点，看看有没有梯度绝对值大于0.1的。现在是向上的鼓包，因此就是找左边绝对值小于-0.1，右边绝对值大于0.1
            %flags_front1和rear1不要求大于9了，但是至少也得有那么两个吧，再搜去。人连两个都没有，那也就别搜了吧。。
            if (sum(flags_front1)>=2&&sum(flags_rear1)>=2)%如果前后各有2个斜率绝对值大于0.1，就表示有那么点儿意思，在去搜。。
                flags_front3=[];flags_rear3=[];flags_front4=[];flags_rear4=[];flags_front5=[];flags_rear5=[];
                for l=1:10%往左找10个点
                    if index_diff(i-l)<-tol1%发现了第一个斜率小于-0.1的（向上的鼓包，左侧，斜率绝对值大于0.1）
                        gt0_1_left=i-l;%记录下来这个点
                        flags_front3=(index_diff(max(1,gt0_1_left-9):gt0_1_left)<-tol1);%看看这个点左边的10个点的情况
                        flags_front3_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)<-tol1);%意义见下。
                        flags_front4=(index_diff(max(1,gt0_1_left-9):gt0_1_left)<-tol2);
                        flags_front4_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)<-tol2);
                        flags_front5=(index_diff(max(1,gt0_1_left-9):gt0_1_left)<-2*tol2);%有没有小于-0.4的
                        break
                    end
                end
                for l=1:10
                    if index_diff(i+l)>tol1%发现了第一个斜率大于0.1的（向上的鼓包，右侧，斜率绝对值大于0.1）
                        gt0_1_right=i+l;%记录下来这个点
                        flags_rear3=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>tol1);%看看这个点右边的10个点的情况
                        flags_rear3_long=(index_diff(gt0_1_right:min(gt0_1_right+14,len-1))>tol1);
                        flags_rear4=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>tol2);
                        flags_rear4_long=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>tol2);
                        flags_rear5=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>2*tol2);
                        break
                    end
                end%如果没有找到斜率绝对值大于0.1的点，那么l会等于10，而flags_front3和flags_front4都是空的。
                if ~isempty(flags_front3)&&~isempty(flags_rear3)%如果两侧都找到了斜率绝对值大于0.1的
                    forced_area_left=[];forced_area_right=[];
                    if (sum(flags_front3)>=10&&sum(flags_rear3)>=10)|| sum(flags_front3_long)+sum(flags_rear3_long)>=28 ...%这次要求严一些，就要求这10个点所有的绝对值都大于0.1，而且有6个大于0.2
                            && (sum(flags_front4)>=6&&sum(flags_rear4)>=6)|| sum(flags_front4_long)+sum(flags_rear4_long)>=10 ...
                            ||sum(flags_front5)+sum(flags_rear5)>=2
                        %调试中，经验确定的判断准则。
                        last1=min([gt0_1_left-1,wid]);last2=min([len-gt0_1_right-1,wid]);%开始搜索gt0_1_left和gt0_1_right两边的点，看看哪个点的斜率绝对值开始不大于0.1了。注意是len-gt0_1_right-1，因为index_diff只有len-1个数。
                        for k1=10:last1
                            if index_diff(gt0_1_left-k1)>-tol1&&index_diff(gt0_1_left-k1)>index_diff(gt0_1_left-k1+1) %左边，和原来那个类似，只不过这儿把那个i换成了gt0_1_left，而且不管右边。
                                forced_area_left=k1;
                                break
                            end
                            if k1==last1%如果搜索到last1还没找到的话。1045R
                                index_diff_left=index_diff(gt0_1_left-k1:gt0_1_left);
                                lst0=length(find(index_diff_left<0));%看看里面有多少个<0的
                                if ismember(i,100:len-1-100) && lst0==length(index_diff_left)%如果全<0，即，在一直往上跑。又加了个i要在中间的第101~261列的条件。（为啥？这是因为后面接着往左搜索了，印象中就是因为怕弄到了第一个而造成错误）。
                                    %所以，后面再加一个if，用来处理gt0_1_left-last1==1的情况，即，到了最左边的情况。
                                    lgt0_2=length(find(index_diff_left<-tol2));%看看里面多少个<-0.2（-tol2）的
                                    if lgt0_2>length(index_diff_left)/2%如果有一半的小于-0.2（斜率比较大）
                                        for k1_1=last1+1:2*last1%续命一个last1，接着往左搜索
                                            if index_diff(gt0_1_left-k1_1)>-tol1&&index_diff(gt0_1_left-k1_1)>index_diff(gt0_1_left-k1_1+1)
                                                forced_area_left=k1_1;
                                                break
                                            end
                                        end
                                    end
                                end
                                if gt0_1_left-k1==1%如果搜索到了第一个点（注意前提是搜索到last1还没找到） 69R
                                    index_diff_left=index_diff(gt0_1_left-k1:gt0_1_left);
                                    lst0_1=length(find(index_diff_left<-0.1));%看看里面有多少个<-0.1的
                                    if lst0_1==length(index_diff_left)
                                        lgt0_2=length(find(index_diff_left<-tol2));
                                        if lgt0_2>length(index_diff_left)/2
                                            forced_area_left=k1;%相当于是从1到gt0_1_left都要被修正。
                                        end
                                    end
                                end
                            end
                        end
                        for k2=10:last2
                            if index_diff(gt0_1_right+k2)<tol1&&index_diff(gt0_1_right+k2)<index_diff(gt0_1_right+k2-1)%右边
                                forced_area_right=k2;
                                break
                            end
                            if k2==last2%如果搜索到last1还没找到的话。跟着上面k1的弄的
                                index_diff_right=index_diff(gt0_1_right:gt0_1_right+k2);
                                rst0=length(find(index_diff_right>0));%看看里面有多少个>0的
                                if ismember(i,100:len-1-100) && rst0==length(index_diff_right)%如果全>0，即，在一直往下跑。又加了个i要在中间的第101~261列的条件。
                                    %所以，后面再加一个if，用来处理gt0_1_right+last2==len-1的情况，即，到了最右边的情况。
                                    lgt0_2=length(find(index_diff_right>tol2));%看看里面多少个>0.2（tol2）的
                                    if lgt0_2>length(index_diff_right)/2%如果有一半的大于0.2（斜率比较大）
                                        for k2_1=last2+1:2*last2%续命一个last2，接着往右搜索
                                            if index_diff(gt0_1_right+k2_1)<tol1&&index_diff(gt0_1_right+k2_1)<index_diff(gt0_1_right+k2_1-1)
                                                forced_area_left=k2_1;
                                                break
                                            end
                                        end
                                    end
                                end
                                if gt0_1_right+last2==len-1%如果搜索到了最后一个点（注意前提是搜索到last2还没找到）
                                    index_diff_right=index_diff(gt0_1_right:gt0_1_right+k2);
                                    lst0_1=length(find(index_diff_right>0.1));%看看里面有多少个>0.1的
                                    if lst0_1==length(index_diff_right)
                                        lgt0_2=length(find(index_diff_right>tol2));
                                        if lgt0_2>length(index_diff_right)/2
                                            forced_area_right=k2;%相当于是从gt0_1_right到len-1都要被修正。
                                        end
                                    end
                                end
                            end
                        end
                        if ~isempty(forced_area_left)&&~isempty(forced_area_right)%如果找到了斜率绝对值不再大于0.1的（如果没找到，不修正，说明可能真的是斜的然后平过来）
                            repair_begin=index(gt0_1_left-forced_area_left);
                            repair_end=index(gt0_1_right+forced_area_right);
                            if (forced_area_left+forced_area_right+gt0_1_right-gt0_1_left-1)<=3*wid%170522，太宽（超过3*wid即90个）就不修正了。
                                derta_y=(repair_end-repair_begin)/(forced_area_left+forced_area_right+gt0_1_right-gt0_1_left);
                                for m=1:(forced_area_left+forced_area_right+gt0_1_right-gt0_1_left-1)
                                    index(gt0_1_left-forced_area_left+m)=index(gt0_1_left-forced_area_left)+m*derta_y;
                                end
                            end
                        end
                        f=1;
                    end
                    if f==0%如果不满足那些条件（那个很长的if，如1065R），说明，要么是真的没鼓包，要么是鼓包比较宽。搜索一下大尺度的范围。其实这些似乎可以结合到前面的判断里，现在不想调了。
                        [biggest_count_front1, longest_pos_front1]=find_consecutive_smaller(index_diff(1:gt0_1_left),-0.1);%左边
                        %注意，左边应该满足：longest_pos_front1是第一个连续地<-0.1的点，longest_pos_front1+biggest_count_front1-1是最后一个<-0.1的点（未必是gt0_1_left本身，除非gt0_1_left就是最长的连续的小于-0.1的点的终点）。
                        [biggest_count_rear1, longest_pos_rear1]=find_consecutive_larger(index_diff(gt0_1_right:len-1),0.1);%右边
                        %注意，右边应该满足：gt0_1_right+longest_pos_rear1-1是第一个连续地>0.1的点（未必是gt0_1_right本身，除非gt0_1_right就是最长的连续的大于0.1的点的起点），gt0_1_right+longest_pos1-1+biggest_count_rear1是最后一个>0.1的点。
                        [biggest_count_front2, ~]=find_consecutive_smaller(index_diff(1:gt0_1_left),-0.2);%左边.意义同上。
                        [biggest_count_rear2, ~]=find_consecutive_larger(index_diff(gt0_1_right:len-1),0.2);%右边.意义同上。
                        if biggest_count_front1>=25 && biggest_count_rear1>=25 && biggest_count_front2>=10 && biggest_count_rear2>=10
                            %暂时要求严一些。但是还是担心会有不用修正的被修正掉。。
                            repair_begin=index(longest_pos_front1);
                            repair_end=index(gt0_1_right+longest_pos_rear1-1+biggest_count_rear1);
                            if (gt0_1_right+longest_pos_rear1-1+biggest_count_rear1-longest_pos_front1-1)<=3*wid%170522，太宽（超过3*wid即90个）就不修正了。
                                derta_y=(repair_end-repair_begin)/(gt0_1_right+longest_pos_rear1-1+biggest_count_rear1-longest_pos_front1);
                                for m=1:(gt0_1_right+longest_pos_rear1-1+biggest_count_rear1-longest_pos_front1-1)
                                    index(longest_pos_front1+m)=index(longest_pos_front1)+m*derta_y;
                                end
                                disp('修正了较宽的鼓包')
                            end
                        end
                    end
                end
            end
        end
    end
    
    if (index_diff(i)>0.001 && index_diff(i+1)<-0.001 || (index_diff(i)>0 && index_diff(i+1)<-0  && index_diff(i-1)>0.05 && index_diff(i+1)<-0.05)) %...
            %&& min(index_diff(min(i+21,len-1):min(i+40,len-1)))<0.5
        %修正向下的鼓包，i是最向下的一点。
        flags_front1=(index_diff(max(1,i-9):i)>tol1);%看看它前面的点，有多少个斜率是大于tol1的
        flags_front2=(index_diff(max(1,i-9):i)>tol2);%看看它前面的点，有多少个斜率是大于tol2的
        flags_rear1=(index_diff(min(i+1,len-1):min(i+10,len-1))<-tol1);%看看它后面的点，有多少个斜率是小于-tol1的
        flags_rear2=(index_diff(min(i+1,len-1):min(i+10,len-1))<-tol2);%看看它后面的点，有多少个斜率是小于-tol2的
        if (sum(flags_front1)>=9&&sum(flags_rear1)>=9) || (sum(flags_front2)>=6&&sum(flags_rear2)>=6)%如果前后各有9个斜率绝对值大于0.1或者6个大于0.2的
            %如果满足以上条件，就认为有个鼓包了。
            %先找修正的范围
            last=min([i-1,len-i-1,wid]);%现在是第i个点，前面有i-1个点，后面有len-i-1个点（因为index_diff只有len-1个数）。设个限制，最多搜索wid个好了。
            for k=10:last%以第i个点为中心搜索，从它两边各10个点开始搜索。如果第k-9次搜索（即第i-k和i+1+k个点）满足：第i-k个点导数小于0.1且小于第i-k+1个点，且第i+k个点导数大于-0.1且大于第i+k-1个点
                if index_diff(i-k)<0.1&&index_diff(i-k)<index_diff(i-k+1) ...%左边：第k次搜索是第i-k个点，如果它的导数绝对值小于0.1（正的，而且小于-0.1），且绝对值比上一次搜索（第i-k+1个点）要小（正的，那个值小）
                        && index_diff(i+k)>-0.1&&index_diff(i+k)>index_diff(i+k-1)%右边：第i+k个点，如果导数绝对值小于0.1（负的，且大于0.1），且绝对值比上一次搜索（第i+k-1个点）要小（负的，那个值大）
                    forced_area=k;%i-forced_area:i+forced_area这个范围内，进行鼓包修正。
                    break
                end
                forced_area=k;%如果搜索了wid个，还没有满足上述条件的，就弄这wid个了。
            end
            %然后开始修正
            repair_begin=index(i-forced_area);
            repair_end=index(i+forced_area);
            if (2*forced_area-1)<=3*wid%170522，太宽（超过3*wid即90个）就不修正了。
                derta_y=(repair_end-repair_begin)/(2*forced_area);%从第i-forced_area个修正到第i+forced_area个。那么，有2*forced_area+1个点要插值。
                for m=1:2*forced_area-1
                    index(i-forced_area+m)=index(i-forced_area)+m*derta_y;
                end
            end
        else
            %flags_front1和rear1不要求大于9了，但是至少也得有那么两个吧，再搜。如果连两个都没有，那也就别搜了吧。。
            if (sum(flags_front1)>=2&&sum(flags_rear1)>=2)%如果前后各有2个斜率绝对值大于0.1，就表示有那么点儿意思，在去搜。。
                flags_front3=[];flags_rear3=[];flags_front4=[];flags_rear4=[];flags_front5=[];flags_rear5=[];
                for l=1:10%往左找10个点
                    if index_diff(i-l)>tol1%发现了第一个斜率大于0.1的（向下的鼓包，左侧，斜率绝对值大于0.1）
                        gt0_1_left=i-l;%记录下来这个点
                        flags_front3=(index_diff(max(1,gt0_1_left-9):gt0_1_left)>tol1);%看看这个点左边的10个点的情况
                        flags_front3_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)>tol1);%意义见下。
                        flags_front4=(index_diff(max(1,gt0_1_left-9):gt0_1_left)>tol2);
                        flags_front4_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)>tol2);
                        flags_front5=(index_diff(max(1,gt0_1_left-9):gt0_1_left)>2*tol2);%有没有大于0.4的
                        break
                    end
                end
                for l=1:10
                    if index_diff(i+l)<-tol1%发现了第一个斜率小于-0.1的（向下的鼓包，右侧，斜率绝对值大于0.1）
                        gt0_1_right=i+l;%记录下来这个点
                        flags_rear3=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-tol1);%看看这个点右边的10个点的情况
                        flags_rear3_long=(index_diff(gt0_1_right:min(gt0_1_right+14,len-1))<-tol1);
                        flags_rear4=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-tol2);
                        flags_rear4_long=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-tol2);
                        flags_rear5=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-2*tol2);
                        break
                    end
                end%如果没有找到斜率绝对值大于0.1的点，那么l会等于10，而flags_front3和flags_front4都是空的。
                if ~isempty(flags_front3)&&~isempty(flags_rear3)%如果两侧都找到了斜率绝对值大于0.1的
                    forced_area_left=[];forced_area_right=[];
                    if (sum(flags_front3)>=10&&sum(flags_rear3)>=10) || sum(flags_front3_long)+sum(flags_rear3_long)>=28 ...%这次要求严一些，就要求这10个点所有的绝对值都大于0.1，而且有6个大于0.2
                            && (sum(flags_front4)>=6&&sum(flags_rear4)>=6) || sum(flags_front4_long)+sum(flags_rear4_long)>=10 ...
                            ||sum(flags_front5)+sum(flags_rear5)>=2
                        %按照经验确定的判断条件。
                        last1=min([gt0_1_left-1,wid]);last2=min([len-gt0_1_right-1,wid]);%开始搜索gt0_1_left和gt0_1_right两边的点，看看哪个点的斜率绝对值开始不大于0.1了。注意是len-gt0_1_right-1，因为index_diff只有len-1个数。
                        for k1=10:last1
                            if index_diff(gt0_1_left-k1)<tol1&&index_diff(gt0_1_left-k1)<index_diff(gt0_1_left-k1+1) %左边，和原来那个类似，只不过这儿把那个i换成了gt0_1_left，而且不管右边。
                                forced_area_left=k1;
                                break
                            end
                        end
                        for k2=10:last2
                            if index_diff(gt0_1_right+k2)>-tol1&&index_diff(gt0_1_right+k2)>-index_diff(gt0_1_right+k2-1)%右边
                                forced_area_right=k2;
                                break
                            end
                        end
                        if ~isempty(forced_area_left)&&~isempty(forced_area_right)
                            repair_begin=index(gt0_1_left-forced_area_left);
                            repair_end=index(gt0_1_right+forced_area_right);
                            if ((forced_area_left+forced_area_right+gt0_1_right-gt0_1_left)-1)<=3*wid%170522，太宽（超过3*wid即90个）就不修正了。
                                derta_y=(repair_end-repair_begin)/(forced_area_left+forced_area_right+gt0_1_right-gt0_1_left);
                                for m=1:(forced_area_left+forced_area_right+gt0_1_right-gt0_1_left)-1
                                    index(gt0_1_left-forced_area_left+m)=index(gt0_1_left-forced_area_left)+m*derta_y;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
flags_1_10_gt_p=(index_diff(1:10)>0.2);%看看前面10个点的情况，有多少个>0.2的
if sum(flags_1_10_gt_p)==10 && inclined~=1%如果这10个都大于0.2（前10个往下斜），且整体并不是往下斜的
    for c1=11:30%往后找20列，找斜率不大于0.1的
        if index_diff(c1)<=0.1 && index_diff(c1)>=0
            mark=c1;
            break
        end
        mark=c1;%如果找了30个还没有，那就这样吧。
    end
    aff_area_graymean=grayscale_mean(1:c1);%
    aff_area_graystd=grayscale_std(1:c1);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*c1 && gt15>=0.9*c1 && gt30_std<=1%灰度均值要比较大，标准差不能太大，才修正。
            index(1:c1-1)=index(mark);%全都平过来？
        else
            red_alert(1)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(1)==0%如果red_alert(1)为1的话，就说明，红线是因为不满足灰度条件而没修正。此时，蓝线也不要修正。
            index(1:c1-1)=index(mark);%全都平过来
        end
    end
end
flags_1_10_gt_n=(index_diff(1:10)<-0.2);%看看前面10个点的情况，有多少个<-0.2的
if sum(flags_1_10_gt_n)==10 && inclined~=-1%如果这10个都小于-0.2（前10个往上斜），且整体并不是往上斜的
    for c2=11:30%往后找20列，找斜率不大于0.1的
        if index_diff(c2)>=-0.1 && index_diff(c2)<=0
            mark=c2;
            break
        end
        mark=c2;%如果找了30个还没有，那就这样吧。
    end
    aff_area_graymean=grayscale_mean(1:c2);%
    aff_area_graystd=grayscale_std(1:c2);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*c2 && gt15>=0.9*c2 && gt30_std<=1
            index(1:c2-1)=index(mark);%全都平过来？
        else
            red_alert(2)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(2)==0%如果red_alert(2)为1的话，就说明，红线是因为不满足灰度条件而没修正。此时，蓝线也不要修正。
            index(1:c2-1)=index(mark);%全都平过来
        end
    end
end
flags_len_10_len_gt_p=(index_diff(len-10:len-1)>0.2);%看看后面10个点的情况，有多少个>0.2的
if sum(flags_len_10_len_gt_p)==10 && inclined~=1%如果这10个都大于0.2（后10个往下斜），且整体并不是往下斜的
    for c3=len-11:-1:len-30%往前找20列，找斜率不大于0.1的
        if index_diff(c3)<=0.1 && index_diff(c3)>=0
            mark=c3;
            break
        end
        mark=c3;%如果找了30个还没有，那就这样吧。
    end
    aff_area_graymean=grayscale_mean(c3:len);%
    aff_area_graystd=grayscale_std(c3:len);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*(len-c3+1) && gt15>=0.9*(len-c3+1) || len-c3+1<=20 && gt30_std<=1%len-c3+1<=20是为了1068L耍赖的。。
            index(c3+1:len)=index(mark);%全都平过来？
        else
            red_alert(3)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(3)==0%如果red_alert(3)为1的话，就说明，红线是因为不满足灰度条件而没修正。此时，蓝线也不要修正。
            index(c3+1:len)=index(mark);%全都平过来
        end
    end
end
flags_len_10_len_gt_n=(index_diff(len-10:len-1)<-0.2);%看看后面10个点的情况，有多少个<-0.2的
if sum(flags_len_10_len_gt_n)==10 && inclined~=-1%如果这10个都小于-0.2（后10个往上斜），且整体并不是往上斜的
    for c4=len-11:-1:len-30%往前找20列，找斜率不大于0.1的
        if index_diff(c4)>=-0.1 && index_diff(c4)<=0
            mark=c4;
            break
        end
        mark=c4;%如果找了30个还没有，那就这样吧。
    end
    aff_area_graymean=grayscale_mean(c4:len);%
    aff_area_graystd=grayscale_std(c4:len);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*(len-c4+1) && gt15>=0.9*(len-c4+1) && gt30_std<=1
            index(c4+1:len)=index(mark);%全都平过来？
        else
            red_alert(4)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(4)==0%如果red_alert(4)为1的话，就说明，红线是因为不满足灰度条件而没修正。此时，蓝线也不要修正。
            index(c4+1:len)=index(mark);%全都平过来
        end
    end
end
repaired=index;