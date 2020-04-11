function [repaired,red_alert]=repair_bumpupanddown(I,unrepaired,tol1,tol2,wid,mode,red_alert)
%��������ݻҶȺ�б�ʣ�������������Ĺİ���б����ֵ��������Ĭ��tol1=0.1��tol2=0.2��widΪ������ȡ�
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
count_gt0_2=length(find(index_diff>0.2));%�ж��ٸ�����0.2
count_lt_0_2=length(find(index_diff<-0.2));%�ж��ٸ�С��-0.2
count_all0_2=count_gt0_2+count_lt_0_2;%����֮�ͣ������ٸ�����ֵ����0.2�ġ�
red_alert_default=[0 0 0 0];
if strcmp(mode,'red'), red_alert=[0 0 0 0];end%����Ǻ��ߵĻ����ͳ�ʼ��Ϊ0��
%�����Ȼǰ��10��������б�������岢��������б�ģ�����������߻Ҷ�����������û�������Ļ����ǾͲ��������ġ�
if strcmp(mode,'blue')%��������ߣ�����û�����������룬ҲŪ��0
    if ~exist('red_alert','var')
        red_alert=red_alert_default;
    end
end
inclined=0;
if (count_all0_2>=80 && count_gt0_2>=0.75*count_all0_2) || (count_all0_2>=60 && count_gt0_2>=0.85*count_all0_2)
    inclined=1;%�������ֵ����0.2����100�У����ң�������80%���ϵĶ�����0.2���Ǿ���Ϊȷʵ������б��
end
if (count_all0_2>=80 && count_lt_0_2>=0.75*count_all0_2) || (count_all0_2>=60 && count_lt_0_2>=0.85*count_all0_2)
    inclined=-1;%�������ֵ����0.2����100�У����ң�������80%���ϵĶ�С��-0.2���Ǿ���Ϊȷʵ������б�ġ�69L�ĳ�80�С�41R���ˡ����ߴ���60�У�������85%���϶���һ������ġ�
end
for i=11:len-11%�м���Щ�У�ע��index_diffֻ��i-1��Ԫ�ء�ǰ10�кͺ�10�е������ǣ�����
    if i==257
        pause(0.001)
    end
    if (index_diff(i)<-0.001 && index_diff(i+1)>0.001 || ...
            (index_diff(i)<-0 && index_diff(i+1)>0 && index_diff(i-1)<-0.04 && index_diff(i+2)>0.04)) 
        flags_front1=(index_diff(max(1,i-9):i)<-tol1);%������ǰ��ĵ㣬�ж��ٸ�б����С��-tol1��
        flags_front2=(index_diff(max(1,i-9):i)<-tol2);%������ǰ��ĵ㣬�ж��ٸ�б����С��-tol2��
        flags_rear1=(index_diff(min(i+1,len-1):min(i+10,len-1))>tol1);%����������ĵ㣬�ж��ٸ�б���Ǵ���tol1��
        flags_rear2=(index_diff(min(i+1,len-1):min(i+10,len-1))>tol2);%����������ĵ㣬�ж��ٸ�б���Ǵ���tol2��
        if (sum(flags_front1)>=9&&sum(flags_rear1)>=9) && (sum(flags_front2)>=6&&sum(flags_rear2)>=6)%���ǰ�����9��б�ʾ���ֵ����0.1��������6������0.2��
            %���������������������Ϊ�и��İ��ˡ�
            %���������ķ�Χ
            last=min([i-1,len-i-1,wid]);%�����ǵ�i���㣬ǰ����i-1���㣬������len-i-1���㣨��Ϊindex_diffֻ��len-1��������������ƣ��������wid�����ˡ�
            for k=10:last%�Ե�i����Ϊ�����������������߸�10���㿪ʼ�����������k-9������������i-k��i+1+k���㣩���㣺��i-k���㵼������-0.1�Ҵ��ڵ�i-k+1���㣨����ֵС��0.1��������������û��ǰһ�������Ĵ󣩣��ҵ�i+k���㵼��С��0.1��С�ڵ�i+k-1����ĵ�����
                if index_diff(i-k)>-tol1&&index_diff(i-k)>index_diff(i-k+1) ...%��ߣ���k�������ǵ�i-k���㣬������ĵ�������ֵС��0.1�����ģ����Ҵ���-0.1�����Ҿ���ֵ����һ����������i-k+1���㣩ҪС�����ģ��Ǹ�ֵ��
                        && index_diff(i+k)<tol1&&index_diff(i+k)<index_diff(i+k-1)%�ұߣ���i+k���㣬�����������ֵС��0.1�����ģ���С��0.1�����Ҿ���ֵ����һ����������i+k-1���㣩ҪС�����ģ��Ǹ�ֵС��
                    %��Ȼ��������һ��Ҫindex_diff(i-k)>-tol1�������һ��������С�����Ǹ��ã�����
                    forced_area=k;%i-forced_area+1:i+forced_area-1�����Χ�ڣ����йİ�����������i-forced_area���͵�i+forced_area�����Ѿ����������˰���
                    break
                end
                forced_area=k;%���������wid������û���������������ģ���Ū��wid���ˡ�
            end
            %Ȼ��ʼ����
            repair_begin=index(i-forced_area);%�õ�i-forced_area���͵�i+forced_area�����y����Ϊ��׼�������������㶼�������ģ���
            repair_end=index(i+forced_area);
            %170522��һ�Σ����(2*forced_area-1)>3*wid�Ͳ������ˡ���������൱�ڹİ�̫��90�����ϰ����ˡ������������ġ�
            if (2*forced_area-1)<=3*wid%170522
                derta_y=(repair_end-repair_begin)/(2*forced_area);%�ӵ�i-forced_area+1����������i+forced_area-1������ô����2*forced_area-1����Ҫ��ֵ�����ǣ�ÿ������֮��ľ��룬Ӧ���ǳ���2*forced_area������2*forced_area-1��
                for m=1:2*forced_area-1
                    index(i-forced_area+m)=index(i-forced_area)+m*derta_y;
                end
            end
        else%���else����������ƽ�Ĺİ������������ƽ���ʹӶ����������������10���㣬������û���ݶȾ���ֵ����0.1�ġ����������ϵĹİ�����˾�������߾���ֵС��-0.1���ұ߾���ֵ����0.1
            %flags_front1��rear1��Ҫ�����9�ˣ���������Ҳ������ô�����ɣ�����ȥ������������û�У���Ҳ�ͱ����˰ɡ���
            if (sum(flags_front1)>=2&&sum(flags_rear1)>=2)%���ǰ�����2��б�ʾ���ֵ����0.1���ͱ�ʾ����ô�����˼����ȥ�ѡ���
                flags_front3=[];flags_rear3=[];flags_front4=[];flags_rear4=[];flags_front5=[];flags_rear5=[];
                for l=1:10%������10����
                    if index_diff(i-l)<-tol1%�����˵�һ��б��С��-0.1�ģ����ϵĹİ�����࣬б�ʾ���ֵ����0.1��
                        gt0_1_left=i-l;%��¼���������
                        flags_front3=(index_diff(max(1,gt0_1_left-9):gt0_1_left)<-tol1);%�����������ߵ�10��������
                        flags_front3_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)<-tol1);%������¡�
                        flags_front4=(index_diff(max(1,gt0_1_left-9):gt0_1_left)<-tol2);
                        flags_front4_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)<-tol2);
                        flags_front5=(index_diff(max(1,gt0_1_left-9):gt0_1_left)<-2*tol2);%��û��С��-0.4��
                        break
                    end
                end
                for l=1:10
                    if index_diff(i+l)>tol1%�����˵�һ��б�ʴ���0.1�ģ����ϵĹİ����Ҳ࣬б�ʾ���ֵ����0.1��
                        gt0_1_right=i+l;%��¼���������
                        flags_rear3=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>tol1);%����������ұߵ�10��������
                        flags_rear3_long=(index_diff(gt0_1_right:min(gt0_1_right+14,len-1))>tol1);
                        flags_rear4=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>tol2);
                        flags_rear4_long=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>tol2);
                        flags_rear5=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))>2*tol2);
                        break
                    end
                end%���û���ҵ�б�ʾ���ֵ����0.1�ĵ㣬��ôl�����10����flags_front3��flags_front4���ǿյġ�
                if ~isempty(flags_front3)&&~isempty(flags_rear3)%������඼�ҵ���б�ʾ���ֵ����0.1��
                    forced_area_left=[];forced_area_right=[];
                    if (sum(flags_front3)>=10&&sum(flags_rear3)>=10)|| sum(flags_front3_long)+sum(flags_rear3_long)>=28 ...%���Ҫ����һЩ����Ҫ����10�������еľ���ֵ������0.1��������6������0.2
                            && (sum(flags_front4)>=6&&sum(flags_rear4)>=6)|| sum(flags_front4_long)+sum(flags_rear4_long)>=10 ...
                            ||sum(flags_front5)+sum(flags_rear5)>=2
                        %�����У�����ȷ�����ж�׼��
                        last1=min([gt0_1_left-1,wid]);last2=min([len-gt0_1_right-1,wid]);%��ʼ����gt0_1_left��gt0_1_right���ߵĵ㣬�����ĸ����б�ʾ���ֵ��ʼ������0.1�ˡ�ע����len-gt0_1_right-1����Ϊindex_diffֻ��len-1������
                        for k1=10:last1
                            if index_diff(gt0_1_left-k1)>-tol1&&index_diff(gt0_1_left-k1)>index_diff(gt0_1_left-k1+1) %��ߣ���ԭ���Ǹ����ƣ�ֻ����������Ǹ�i������gt0_1_left�����Ҳ����ұߡ�
                                forced_area_left=k1;
                                break
                            end
                            if k1==last1%���������last1��û�ҵ��Ļ���1045R
                                index_diff_left=index_diff(gt0_1_left-k1:gt0_1_left);
                                lst0=length(find(index_diff_left<0));%���������ж��ٸ�<0��
                                if ismember(i,100:len-1-100) && lst0==length(index_diff_left)%���ȫ<0��������һֱ�����ܡ��ּ��˸�iҪ���м�ĵ�101~261�е���������Ϊɶ��������Ϊ����������������ˣ�ӡ���о�����Ϊ��Ū���˵�һ������ɴ��󣩡�
                                    %���ԣ������ټ�һ��if����������gt0_1_left-last1==1�������������������ߵ������
                                    lgt0_2=length(find(index_diff_left<-tol2));%����������ٸ�<-0.2��-tol2����
                                    if lgt0_2>length(index_diff_left)/2%�����һ���С��-0.2��б�ʱȽϴ�
                                        for k1_1=last1+1:2*last1%����һ��last1��������������
                                            if index_diff(gt0_1_left-k1_1)>-tol1&&index_diff(gt0_1_left-k1_1)>index_diff(gt0_1_left-k1_1+1)
                                                forced_area_left=k1_1;
                                                break
                                            end
                                        end
                                    end
                                end
                                if gt0_1_left-k1==1%����������˵�һ���㣨ע��ǰ����������last1��û�ҵ��� 69R
                                    index_diff_left=index_diff(gt0_1_left-k1:gt0_1_left);
                                    lst0_1=length(find(index_diff_left<-0.1));%���������ж��ٸ�<-0.1��
                                    if lst0_1==length(index_diff_left)
                                        lgt0_2=length(find(index_diff_left<-tol2));
                                        if lgt0_2>length(index_diff_left)/2
                                            forced_area_left=k1;%�൱���Ǵ�1��gt0_1_left��Ҫ��������
                                        end
                                    end
                                end
                            end
                        end
                        for k2=10:last2
                            if index_diff(gt0_1_right+k2)<tol1&&index_diff(gt0_1_right+k2)<index_diff(gt0_1_right+k2-1)%�ұ�
                                forced_area_right=k2;
                                break
                            end
                            if k2==last2%���������last1��û�ҵ��Ļ�����������k1��Ū��
                                index_diff_right=index_diff(gt0_1_right:gt0_1_right+k2);
                                rst0=length(find(index_diff_right>0));%���������ж��ٸ�>0��
                                if ismember(i,100:len-1-100) && rst0==length(index_diff_right)%���ȫ>0��������һֱ�����ܡ��ּ��˸�iҪ���м�ĵ�101~261�е�������
                                    %���ԣ������ټ�һ��if����������gt0_1_right+last2==len-1������������������ұߵ������
                                    lgt0_2=length(find(index_diff_right>tol2));%����������ٸ�>0.2��tol2����
                                    if lgt0_2>length(index_diff_right)/2%�����һ��Ĵ���0.2��б�ʱȽϴ�
                                        for k2_1=last2+1:2*last2%����һ��last2��������������
                                            if index_diff(gt0_1_right+k2_1)<tol1&&index_diff(gt0_1_right+k2_1)<index_diff(gt0_1_right+k2_1-1)
                                                forced_area_left=k2_1;
                                                break
                                            end
                                        end
                                    end
                                end
                                if gt0_1_right+last2==len-1%��������������һ���㣨ע��ǰ����������last2��û�ҵ���
                                    index_diff_right=index_diff(gt0_1_right:gt0_1_right+k2);
                                    lst0_1=length(find(index_diff_right>0.1));%���������ж��ٸ�>0.1��
                                    if lst0_1==length(index_diff_right)
                                        lgt0_2=length(find(index_diff_right>tol2));
                                        if lgt0_2>length(index_diff_right)/2
                                            forced_area_right=k2;%�൱���Ǵ�gt0_1_right��len-1��Ҫ��������
                                        end
                                    end
                                end
                            end
                        end
                        if ~isempty(forced_area_left)&&~isempty(forced_area_right)%����ҵ���б�ʾ���ֵ���ٴ���0.1�ģ����û�ҵ�����������˵�����������б��Ȼ��ƽ������
                            repair_begin=index(gt0_1_left-forced_area_left);
                            repair_end=index(gt0_1_right+forced_area_right);
                            if (forced_area_left+forced_area_right+gt0_1_right-gt0_1_left-1)<=3*wid%170522��̫������3*wid��90�����Ͳ������ˡ�
                                derta_y=(repair_end-repair_begin)/(forced_area_left+forced_area_right+gt0_1_right-gt0_1_left);
                                for m=1:(forced_area_left+forced_area_right+gt0_1_right-gt0_1_left-1)
                                    index(gt0_1_left-forced_area_left+m)=index(gt0_1_left-forced_area_left)+m*derta_y;
                                end
                            end
                        end
                        f=1;
                    end
                    if f==0%�����������Щ�������Ǹ��ܳ���if����1065R����˵����Ҫô�����û�İ���Ҫô�ǹİ��ȽϿ�����һ�´�߶ȵķ�Χ����ʵ��Щ�ƺ����Խ�ϵ�ǰ����ж�����ڲ�����ˡ�
                        [biggest_count_front1, longest_pos_front1]=find_consecutive_smaller(index_diff(1:gt0_1_left),-0.1);%���
                        %ע�⣬���Ӧ�����㣺longest_pos_front1�ǵ�һ��������<-0.1�ĵ㣬longest_pos_front1+biggest_count_front1-1�����һ��<-0.1�ĵ㣨δ����gt0_1_left��������gt0_1_left�������������С��-0.1�ĵ���յ㣩��
                        [biggest_count_rear1, longest_pos_rear1]=find_consecutive_larger(index_diff(gt0_1_right:len-1),0.1);%�ұ�
                        %ע�⣬�ұ�Ӧ�����㣺gt0_1_right+longest_pos_rear1-1�ǵ�һ��������>0.1�ĵ㣨δ����gt0_1_right��������gt0_1_right������������Ĵ���0.1�ĵ����㣩��gt0_1_right+longest_pos1-1+biggest_count_rear1�����һ��>0.1�ĵ㡣
                        [biggest_count_front2, ~]=find_consecutive_smaller(index_diff(1:gt0_1_left),-0.2);%���.����ͬ�ϡ�
                        [biggest_count_rear2, ~]=find_consecutive_larger(index_diff(gt0_1_right:len-1),0.2);%�ұ�.����ͬ�ϡ�
                        if biggest_count_front1>=25 && biggest_count_rear1>=25 && biggest_count_front2>=10 && biggest_count_rear2>=10
                            %��ʱҪ����һЩ�����ǻ��ǵ��Ļ��в��������ı�����������
                            repair_begin=index(longest_pos_front1);
                            repair_end=index(gt0_1_right+longest_pos_rear1-1+biggest_count_rear1);
                            if (gt0_1_right+longest_pos_rear1-1+biggest_count_rear1-longest_pos_front1-1)<=3*wid%170522��̫������3*wid��90�����Ͳ������ˡ�
                                derta_y=(repair_end-repair_begin)/(gt0_1_right+longest_pos_rear1-1+biggest_count_rear1-longest_pos_front1);
                                for m=1:(gt0_1_right+longest_pos_rear1-1+biggest_count_rear1-longest_pos_front1-1)
                                    index(longest_pos_front1+m)=index(longest_pos_front1)+m*derta_y;
                                end
                                disp('�����˽Ͽ�Ĺİ�')
                            end
                        end
                    end
                end
            end
        end
    end
    
    if (index_diff(i)>0.001 && index_diff(i+1)<-0.001 || (index_diff(i)>0 && index_diff(i+1)<-0  && index_diff(i-1)>0.05 && index_diff(i+1)<-0.05)) %...
            %&& min(index_diff(min(i+21,len-1):min(i+40,len-1)))<0.5
        %�������µĹİ���i�������µ�һ�㡣
        flags_front1=(index_diff(max(1,i-9):i)>tol1);%������ǰ��ĵ㣬�ж��ٸ�б���Ǵ���tol1��
        flags_front2=(index_diff(max(1,i-9):i)>tol2);%������ǰ��ĵ㣬�ж��ٸ�б���Ǵ���tol2��
        flags_rear1=(index_diff(min(i+1,len-1):min(i+10,len-1))<-tol1);%����������ĵ㣬�ж��ٸ�б����С��-tol1��
        flags_rear2=(index_diff(min(i+1,len-1):min(i+10,len-1))<-tol2);%����������ĵ㣬�ж��ٸ�б����С��-tol2��
        if (sum(flags_front1)>=9&&sum(flags_rear1)>=9) || (sum(flags_front2)>=6&&sum(flags_rear2)>=6)%���ǰ�����9��б�ʾ���ֵ����0.1����6������0.2��
            %���������������������Ϊ�и��İ��ˡ�
            %���������ķ�Χ
            last=min([i-1,len-i-1,wid]);%�����ǵ�i���㣬ǰ����i-1���㣬������len-i-1���㣨��Ϊindex_diffֻ��len-1��������������ƣ��������wid�����ˡ�
            for k=10:last%�Ե�i����Ϊ�����������������߸�10���㿪ʼ�����������k-9������������i-k��i+1+k���㣩���㣺��i-k���㵼��С��0.1��С�ڵ�i-k+1���㣬�ҵ�i+k���㵼������-0.1�Ҵ��ڵ�i+k-1����
                if index_diff(i-k)<0.1&&index_diff(i-k)<index_diff(i-k+1) ...%��ߣ���k�������ǵ�i-k���㣬������ĵ�������ֵС��0.1�����ģ�����С��-0.1�����Ҿ���ֵ����һ����������i-k+1���㣩ҪС�����ģ��Ǹ�ֵС��
                        && index_diff(i+k)>-0.1&&index_diff(i+k)>index_diff(i+k-1)%�ұߣ���i+k���㣬�����������ֵС��0.1�����ģ��Ҵ���0.1�����Ҿ���ֵ����һ����������i+k-1���㣩ҪС�����ģ��Ǹ�ֵ��
                    forced_area=k;%i-forced_area:i+forced_area�����Χ�ڣ����йİ�������
                    break
                end
                forced_area=k;%���������wid������û���������������ģ���Ū��wid���ˡ�
            end
            %Ȼ��ʼ����
            repair_begin=index(i-forced_area);
            repair_end=index(i+forced_area);
            if (2*forced_area-1)<=3*wid%170522��̫������3*wid��90�����Ͳ������ˡ�
                derta_y=(repair_end-repair_begin)/(2*forced_area);%�ӵ�i-forced_area����������i+forced_area������ô����2*forced_area+1����Ҫ��ֵ��
                for m=1:2*forced_area-1
                    index(i-forced_area+m)=index(i-forced_area)+m*derta_y;
                end
            end
        else
            %flags_front1��rear1��Ҫ�����9�ˣ���������Ҳ������ô�����ɣ����ѡ������������û�У���Ҳ�ͱ����˰ɡ���
            if (sum(flags_front1)>=2&&sum(flags_rear1)>=2)%���ǰ�����2��б�ʾ���ֵ����0.1���ͱ�ʾ����ô�����˼����ȥ�ѡ���
                flags_front3=[];flags_rear3=[];flags_front4=[];flags_rear4=[];flags_front5=[];flags_rear5=[];
                for l=1:10%������10����
                    if index_diff(i-l)>tol1%�����˵�һ��б�ʴ���0.1�ģ����µĹİ�����࣬б�ʾ���ֵ����0.1��
                        gt0_1_left=i-l;%��¼���������
                        flags_front3=(index_diff(max(1,gt0_1_left-9):gt0_1_left)>tol1);%�����������ߵ�10��������
                        flags_front3_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)>tol1);%������¡�
                        flags_front4=(index_diff(max(1,gt0_1_left-9):gt0_1_left)>tol2);
                        flags_front4_long=(index_diff(max(1,gt0_1_left-14):gt0_1_left)>tol2);
                        flags_front5=(index_diff(max(1,gt0_1_left-9):gt0_1_left)>2*tol2);%��û�д���0.4��
                        break
                    end
                end
                for l=1:10
                    if index_diff(i+l)<-tol1%�����˵�һ��б��С��-0.1�ģ����µĹİ����Ҳ࣬б�ʾ���ֵ����0.1��
                        gt0_1_right=i+l;%��¼���������
                        flags_rear3=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-tol1);%����������ұߵ�10��������
                        flags_rear3_long=(index_diff(gt0_1_right:min(gt0_1_right+14,len-1))<-tol1);
                        flags_rear4=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-tol2);
                        flags_rear4_long=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-tol2);
                        flags_rear5=(index_diff(gt0_1_right:min(gt0_1_right+9,len-1))<-2*tol2);
                        break
                    end
                end%���û���ҵ�б�ʾ���ֵ����0.1�ĵ㣬��ôl�����10����flags_front3��flags_front4���ǿյġ�
                if ~isempty(flags_front3)&&~isempty(flags_rear3)%������඼�ҵ���б�ʾ���ֵ����0.1��
                    forced_area_left=[];forced_area_right=[];
                    if (sum(flags_front3)>=10&&sum(flags_rear3)>=10) || sum(flags_front3_long)+sum(flags_rear3_long)>=28 ...%���Ҫ����һЩ����Ҫ����10�������еľ���ֵ������0.1��������6������0.2
                            && (sum(flags_front4)>=6&&sum(flags_rear4)>=6) || sum(flags_front4_long)+sum(flags_rear4_long)>=10 ...
                            ||sum(flags_front5)+sum(flags_rear5)>=2
                        %���վ���ȷ�����ж�������
                        last1=min([gt0_1_left-1,wid]);last2=min([len-gt0_1_right-1,wid]);%��ʼ����gt0_1_left��gt0_1_right���ߵĵ㣬�����ĸ����б�ʾ���ֵ��ʼ������0.1�ˡ�ע����len-gt0_1_right-1����Ϊindex_diffֻ��len-1������
                        for k1=10:last1
                            if index_diff(gt0_1_left-k1)<tol1&&index_diff(gt0_1_left-k1)<index_diff(gt0_1_left-k1+1) %��ߣ���ԭ���Ǹ����ƣ�ֻ����������Ǹ�i������gt0_1_left�����Ҳ����ұߡ�
                                forced_area_left=k1;
                                break
                            end
                        end
                        for k2=10:last2
                            if index_diff(gt0_1_right+k2)>-tol1&&index_diff(gt0_1_right+k2)>-index_diff(gt0_1_right+k2-1)%�ұ�
                                forced_area_right=k2;
                                break
                            end
                        end
                        if ~isempty(forced_area_left)&&~isempty(forced_area_right)
                            repair_begin=index(gt0_1_left-forced_area_left);
                            repair_end=index(gt0_1_right+forced_area_right);
                            if ((forced_area_left+forced_area_right+gt0_1_right-gt0_1_left)-1)<=3*wid%170522��̫������3*wid��90�����Ͳ������ˡ�
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
flags_1_10_gt_p=(index_diff(1:10)>0.2);%����ǰ��10�����������ж��ٸ�>0.2��
if sum(flags_1_10_gt_p)==10 && inclined~=1%�����10��������0.2��ǰ10������б���������岢��������б��
    for c1=11:30%������20�У���б�ʲ�����0.1��
        if index_diff(c1)<=0.1 && index_diff(c1)>=0
            mark=c1;
            break
        end
        mark=c1;%�������30����û�У��Ǿ������ɡ�
    end
    aff_area_graymean=grayscale_mean(1:c1);%
    aff_area_graystd=grayscale_std(1:c1);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*c1 && gt15>=0.9*c1 && gt30_std<=1%�ҶȾ�ֵҪ�Ƚϴ󣬱�׼���̫�󣬲�������
            index(1:c1-1)=index(mark);%ȫ��ƽ������
        else
            red_alert(1)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(1)==0%���red_alert(1)Ϊ1�Ļ�����˵������������Ϊ������Ҷ�������û��������ʱ������Ҳ��Ҫ������
            index(1:c1-1)=index(mark);%ȫ��ƽ����
        end
    end
end
flags_1_10_gt_n=(index_diff(1:10)<-0.2);%����ǰ��10�����������ж��ٸ�<-0.2��
if sum(flags_1_10_gt_n)==10 && inclined~=-1%�����10����С��-0.2��ǰ10������б���������岢��������б��
    for c2=11:30%������20�У���б�ʲ�����0.1��
        if index_diff(c2)>=-0.1 && index_diff(c2)<=0
            mark=c2;
            break
        end
        mark=c2;%�������30����û�У��Ǿ������ɡ�
    end
    aff_area_graymean=grayscale_mean(1:c2);%
    aff_area_graystd=grayscale_std(1:c2);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*c2 && gt15>=0.9*c2 && gt30_std<=1
            index(1:c2-1)=index(mark);%ȫ��ƽ������
        else
            red_alert(2)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(2)==0%���red_alert(2)Ϊ1�Ļ�����˵������������Ϊ������Ҷ�������û��������ʱ������Ҳ��Ҫ������
            index(1:c2-1)=index(mark);%ȫ��ƽ����
        end
    end
end
flags_len_10_len_gt_p=(index_diff(len-10:len-1)>0.2);%��������10�����������ж��ٸ�>0.2��
if sum(flags_len_10_len_gt_p)==10 && inclined~=1%�����10��������0.2����10������б���������岢��������б��
    for c3=len-11:-1:len-30%��ǰ��20�У���б�ʲ�����0.1��
        if index_diff(c3)<=0.1 && index_diff(c3)>=0
            mark=c3;
            break
        end
        mark=c3;%�������30����û�У��Ǿ������ɡ�
    end
    aff_area_graymean=grayscale_mean(c3:len);%
    aff_area_graystd=grayscale_std(c3:len);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*(len-c3+1) && gt15>=0.9*(len-c3+1) || len-c3+1<=20 && gt30_std<=1%len-c3+1<=20��Ϊ��1068Lˣ���ġ���
            index(c3+1:len)=index(mark);%ȫ��ƽ������
        else
            red_alert(3)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(3)==0%���red_alert(3)Ϊ1�Ļ�����˵������������Ϊ������Ҷ�������û��������ʱ������Ҳ��Ҫ������
            index(c3+1:len)=index(mark);%ȫ��ƽ����
        end
    end
end
flags_len_10_len_gt_n=(index_diff(len-10:len-1)<-0.2);%��������10�����������ж��ٸ�<-0.2��
if sum(flags_len_10_len_gt_n)==10 && inclined~=-1%�����10����С��-0.2����10������б���������岢��������б��
    for c4=len-11:-1:len-30%��ǰ��20�У���б�ʲ�����0.1��
        if index_diff(c4)>=-0.1 && index_diff(c4)<=0
            mark=c4;
            break
        end
        mark=c4;%�������30����û�У��Ǿ������ɡ�
    end
    aff_area_graymean=grayscale_mean(c4:len);%
    aff_area_graystd=grayscale_std(c4:len);
    if strcmp(mode,'red')
        gt20=length(find(aff_area_graymean>=20));
        gt15=length(find(aff_area_graymean>=15));
        gt30_std=length(find(aff_area_graystd>=30));
        if gt20>=0.6*(len-c4+1) && gt15>=0.9*(len-c4+1) && gt30_std<=1
            index(c4+1:len)=index(mark);%ȫ��ƽ������
        else
            red_alert(4)=1;
        end
    end
    if strcmp(mode,'blue')
        if red_alert(4)==0%���red_alert(4)Ϊ1�Ļ�����˵������������Ϊ������Ҷ�������û��������ʱ������Ҳ��Ҫ������
            index(c4+1:len)=index(mark);%ȫ��ƽ����
        end
    end
end
repaired=index;