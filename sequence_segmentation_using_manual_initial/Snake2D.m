function [P,Q,J1,fast_r,fast_b]=Snake2D(fr,I,P,Q,Options,neighbor_in,neighbor_out_r,neighbor_out_b,...
    grayscale_rmax_up_mean,grayscale_rmax_down_mean,grayscale_bmax_up_mean,grayscale_bmax_down_mean,...
    inside_IMC_is_bright,inside_IMC_is_dark,severe_noise,medium_noise,exists_noise,...
    severe_fuzzy,medium_fuzzy,exists_fuzzy,IMW,red_dot_ref,blue_dot_ref,fast_r,fast_b,...
    lightflag,lightflag1,noise_flag0,noise_flag1,noise_flag2,fuzzy_flag1)
% ----------------------------------------------����----------------------------------------------
%   I    ��ǰ֡ͼ�����
%   P   ��һʱ�̵ĺ�����ƣ�LI����ɫ�ߣ�����Ϊ��һʱ��LI�ߵĳ�ֵ��
%   Q   ��һʱ�̵ĺ�����ƣ�MA����ɫ�ߣ�����Ϊ��һʱ��MA�ߵĳ�ֵ��
%   Options   �˲�������ѡ�
%   neighbor_in/neighbor_out_r/neighbor_out_b   ��֡�е�������LI������MA������
%   grayscale_rmax_up_mean/grayscale_rmax_down_mean   ��ʼLI��/�·������ڵ�ƽ���Ҷȡ�
%   grayscale_bmax_up_mean/grayscale_bmax_down_mean   ��ʼMA��/�·������ڵ�ƽ���Ҷȡ�
%   inside_IMC_is_bright/inside_IMC_is_dark   ��ʼLI/MA֮�������ڣ��Ҷ��Ƿ���߻���͡�
%   severe_noise/medium_noise/exists_noise    ��ʼLI�Ϸ����Ƿ��������������е����������������
%                                             ����ȷ����ʱ����Grayscale constraint energy��
%   severe_fuzzy,medium_fuzzy,exists_fuzzy    ��ʼMA�·����Ƿ��������������е����������������
%                                             ����ȷ����ʱ����Grayscale constraint energy��
%   IMW   ��ʼ������Ĥ��ȣ�LI-MA��ࣩ
%   red_dot_ref/blue_dot_ref   ��ʼ��LI��MAÿ�����б�ʣ�����Derivative constraint energy��
%   fast_r/fast_b   �Ƿ���LI/MA��ĳЩ֡���˶��ٶȽϿ졣
%   lightflag/lightflag1/noise_flag0/noise_flag1/noise_flag2/fuzzy_flag1   �Ƿ������������������
%                                       ע�⣬�����֮ǰ��һЩ���ã�������ؿ������������������ܵ��¹���ϣ������ұ�����
% ----------------------------------------------���----------------------------------------------
%   P : ��һʱ��LI���죩�ߵĽ��������һʱ�̵Ĺ۲�ֵ
%   Q : ��һʱ��MA�������ߵĽ��������һʱ�̵Ĺ۲�ֵ
%   J : Binary image with the segmented region
%   fast_r/fast_b   �Ƿ���LI/MA��ĳЩ֡���˶��ٶȽϿ졣

% -----------------------0���������롣-----------------------
defaultoptions=struct('Verbose',false,'Wline',0.9,'Wedge',10,'Wterm',2,'Sigma1',1,'Sigma2',1,'Alpha',4,'Beta',2,'Delta',0.1,'Gamma',0.0025,'Kappa',1,'Iterations',200,'GIterations',0,'Mu',0.2,'Sigma3',1);
if(~exist('Options','var'))
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);%
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end%
    if(length(tags)~=length(fieldnames(Options))) 
        warning('snake:unknownoption','unknown options found');
    end
end
% ͼ����double���͡�
I = double(I);
[ylen, xlen]=size(I);
% 1 ��ʼ�����ߵ��������������е�4ʽ�ĵ����������������ͼ��I������Eext�������������������е�Gext��
% ��J��ͼ���y�ĵ���������MATLAB��x��y�ǵߵ��ģ��е��Ի󣩡�
[Eext,J] = ExternalForceImage2D_1(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1,Options.Verbose);
if(Options.Verbose)  
    figure,imshow(Eext,[]),title('Eextern')  
end
if(Options.Verbose)  
    figure,imshow(Eext,[]),title('Eextern with panelty')  
end
% 2 ��ƫ����һ��ͼ�ϴ�ĵط������Ͻϴ�ĳͷ����������Ǻ���Ҫ������ûд��
Estereo=stereo(ylen,xlen,P,Q);
Eext=Eext+0.001*Estereo;
if(nargout>1)
     J1=DrawSegmentedArea2D(P,size(I));
end   
% 3 ��ʼ����Eext���ֱ��LI��MA������
% 3.1 ��ʼ�������ߵ�Eext��Ggc��Gdcֵ��������ֵ��
volcano=200;  % �����е�Ggc����Gdcֵ��
Eext_red=Eext;
Eext_blue=Eext;
eliminate_x_force=0;
yi_red=fix(P(:,2));%LI����
yi_blue=fix(Q(:,2));%MA���ߡ�
red_dot=P(:,2);blue_dot=Q(:,2);%��ʼ������ʵ��Ū��0Ҳ���Եġ�
for i=1:xlen
    % һЩ�����Ĵ���
    % ���ȷֱ�ΪLI��MA����һ�������������ڽ�������������ֵ����Ҫ�ö��ߵ�ͼ���ݶȻ���Ӱ�졣
    Eext_red(yi_blue(i)-fix(IMW(i)/2):ylen,i)=200;%��Eext_red�У�ԭ��LI-MA���м������£�����MA�ߣ��Ĳ��֣���������ֵ����Ҫ��LI���ܵ���������
    Eext_blue(1:yi_red(i)+fix(IMW(i)/2),i)=200;%���ƣ���Eext_blue�У�LI-MA���м����Ϸ��Ĳ��֣���������ֵ��
    % Ȼ�󣬼���LI��MAÿ����ĵ���������Ὣ������Derivative constraint energy������
    if i~=xlen
        red_dot(i)=P(i+1,2)-P(i,2);%���㵼��
        blue_dot(i)=Q(i+1,2)-Q(i,2);
    else
        red_dot(i)=red_dot(i-1);blue_dot(i)=blue_dot(i-1);
    end
end
%��������������ֵ����ʼֵ����Ϊ��һ��ͼ����Ӧֵ��
push_r_up_thr=grayscale_rmax_up_mean;
push_r_down_thr=grayscale_rmax_down_mean;
push_b_up_thr=grayscale_bmax_up_mean;
push_b_down_thr=grayscale_bmax_down_mean;
% �����ĸ�����ʾ������ֵ��push_r_up_thr�����٣�����Ҫ����Ggc�����15 50 40 50ʲô�ģ����ǿ��˺ö�����У����Գ����ġ�
% ������ٿ���������ʱ���������һ��п��ܹ���ϣ���ʱ����ô��ʱ������Ѫ������ʱ�ǲ�������ѧϰ������������
% ����ܰ�����ǿ��ѧϰ����������Զ�ѵ����Щ��ֵ�������Զ��趨�����ֵ��Ӧ�û�����й���ǿ�öࡣ
margin_r_push_up=15*ones(1,xlen);
margin_b_push_up=50*ones(1,xlen);
margin_r_push_down=40*ones(1,xlen);
margin_b_push_down=50*ones(1,xlen);
% ����Iymaxʲô�ģ��������ҵ�LI/MA��������ĻҶ��ݶ����ֵλ�õġ�������ֵλ�ÿ������ڸ����ж��Ƿ�Ҫ����Grayscale constraint energy��
% �����õ�����Iymax_real_r��Iymax_real_b�������Ķ�����ʱ������
Iymax_r=zeros(1,xlen);Iymax_b=zeros(1,xlen);
Iymax_temp_r=zeros(1,xlen);Iymax_temp_b=zeros(1,xlen);
Iymax_real_r=zeros(1,xlen);Iymax_real_b=zeros(1,xlen);
Iymax_real_temp_r=zeros(1,xlen);Iymax_real_temp_b=zeros(1,xlen);
% �����ĸ�����LI/MA�Ҷ��ݶ����ֵλ�õġ����������ڵġ��ҶȾ�ֵ����ʼ��Ϊ0��
mean_up_r=zeros(1,xlen);mean_up_b=zeros(1,xlen);
mean_down_r=zeros(1,xlen);mean_down_b=zeros(1,xlen);
% �����ĸ����������LI/MAλ�õġ����������ڵġ��ҶȾ�ֵ����ʼ��Ϊ0��
mean_red_up=zeros(1,xlen);mean_blue_up=zeros(1,xlen);
mean_red_down=zeros(1,xlen);mean_blue_down=zeros(1,xlen);
% 3.2 ����ÿһ�У��ֱ�ΪLI/MA�ҵ�ԭ��λ�ø�����search_l�����ص��ڣ��ġ����ϻҶ�Ҫ��ġ��Ҷ��ݶ����ֵ���y���ꡣ
% Ŀ���Ǹ��ݻҶ��ݶ����ֵλ�õ���Ϣ�������ж��Ƿ�Ҫ����Grayscale constraint energy��
% Ŀǰ��Щ�������Ǹ��ݾ���ѡȡ�ģ��ô��ǿ��˺ܶ��ͼ��Ӧ�������ıȽ�³�����������ǵ�������̫���
% ����ƪ��������û��д��������ж�׼��
invalid_r_up=zeros(1,xlen);invalid_b_up=zeros(1,xlen);invalid_r_down=zeros(1,xlen);invalid_b_down=zeros(1,xlen);
for i=1:xlen
    search_l=min(8,fix(IMW(i)/2));
    Ji_r=J(yi_red(i)-search_l:yi_red(i)+search_l,i);  % LI
    Ji_b=J(yi_blue(i)-search_l:yi_blue(i)+search_l,i);  % MA
    Ji_rs=sort(Ji_r);%����������Ji_rs�У����һ����Ji_r�е����ֵ�����Ҷ��ݶ����ֵ��
    Ji_bs=sort(Ji_b);
    Jilen=size(Ji_r,1);
    for j=Jilen:-1:1%����LI���ӻҶ��ݶ����Ŀ�ʼ��̽���Ϳ����ǲ�������Ҷ�����
        k=find(Ji_r==Ji_rs(j));%Ji_rs(j)�����ֵ��k�������ֵ��ֵ��λ�á�
        if size(k,1)>1%�����ֹһ�����ֵ����ô�������������������
            %��Ϊ���������һ����������LI�����ܵúܿ죬�ܵ���ɫ������ȥ�ˣ����Ծ���Iymax_r(i)>search_l������ֱ�����ơ�
            Iymax_r(i)=search_l+1;
            break;
        end
        Iymax_temp_r(i)=k;
        Iymax_real_temp_r(i)=Iymax_temp_r(i)+yi_red(i)-search_l;
        mean_up_r_temp=mean2(I(Iymax_real_temp_r(i)-neighbor_out_r(i)-1:Iymax_real_temp_r(i)-1,max(1,i-2):min(i+2,xlen)));
        mean_down_r_temp=mean2(I(Iymax_real_temp_r(i)+1:Iymax_real_temp_r(i)+neighbor_in(i)+1,max(1,i-2):min(i+2,xlen)));
        if abs(mean_up_r_temp-grayscale_rmax_up_mean(i))<abs(mean_up_r_temp-grayscale_bmax_up_mean(i)) ...
                && abs(mean_down_r_temp-grayscale_rmax_down_mean(i))<abs(mean_down_r_temp-grayscale_bmax_down_mean(i))...
                && abs(Q(i,2)-Iymax_real_temp_r(i))>abs(P(i,2)-Iymax_real_temp_r(i))
            % ���if���ж�����Ҷ��ݶ����ֵ�㣬�Ƿ�������������һЩ�Ҷ��������Ϳ��ܲ�����
            Iymax_r(i)=k;%�������Ҷ�������˵��Ji_r�еĵ�k����ȷʵ����LI�ĻҶ��ݶ���󴦡�
            Iymax_real_r(i)=Iymax_r(i)+yi_red(i)-search_l;
            mean_up_r(i)=mean_up_r_temp;
            mean_down_r(i)=mean_down_r_temp;
            break;%�������Ҷ������Ļ�����ô������for j=Jilen:-1:1���ѭ����
        end
        %���������Ҷ���������ôj�ͼ���1���ҵڶ�����������¿��Ƿ�����Ҷ�������   
        if j==1
            invalid_r_up(i)=1;invalid_r_down(i)=1;
            Iymax_r(i)=search_l+1;%�������������м���Ǹ��㣬�൱�ڻָ���ʼλ���ˡ�
            Iymax_real_r(i)=P(i,2);%Ҳ�ָ���ʼλ�ã�����Ӧ����yi_red(i)+1�ģ������ǻ��������ָ���P���ˡ�
        end%����ҵ����һ����û�ҵ����������Ǹ�if�ĵ㣬�Ǿ�˵���Ҷ��ݶ����ֵ�ķ�����Ч�������Ͳ��������ж��ˡ�
    end
    for j=Jilen:-1:1%�������棬���ڴ���MA��
        k=find(Ji_b==Ji_bs(j));%Ji_bs(j)�����ֵ����仰��˵����Ji_b���ҵ�����Ji_bs(j)��Ҳ�������ֵ��ֵ��λ�á�
        Iymax_temp_b(i)=k;
        Iymax_real_temp_b(i)=Iymax_temp_b(i)+yi_blue(i)-search_l;
        mean_up_b_temp=mean2(I(Iymax_real_temp_b(i)-neighbor_in(i)-1:Iymax_real_temp_b(i)-1,max(1,i-2):min(i+2,xlen)));
        mean_down_b_temp=mean2(I(Iymax_real_temp_b(i)+1:Iymax_real_temp_b(i)+neighbor_out_b(i)+1,max(1,i-2):min(i+2,xlen)));
        if abs(mean_up_b_temp-grayscale_bmax_up_mean(i))<abs(mean_up_b_temp-grayscale_rmax_up_mean(i)) ...
                && abs(mean_down_b_temp-grayscale_bmax_down_mean(i))<abs(mean_down_b_temp-grayscale_rmax_down_mean(i))...
                && abs(P(i,2)-Iymax_real_temp_b(i))>abs(Q(i,2)-Iymax_real_temp_b(i)) %20170502�������������Ҫ�������߸�Զ��
            Iymax_b(i)=k;%����ҶȺ����ߵĸ��ӽ���˵��Ji_r�еĵ�k����ȷʵ�������ߵĻҶ��ݶ���󴦡�
            Iymax_real_b(i)=Iymax_b(i)+yi_blue(i)-search_l;
            mean_up_b(i)=mean_up_b_temp;
            mean_down_b(i)=mean_down_b_temp;
            break;%�������Ҷ������Ļ�����ô������for j=Jilen:-1:1���ѭ����
        end
        %���������Ҷ���������ôj�ͼ���1���ҵڶ�����������¿��Ƿ�����Ҷ������� 
        if j==1
            invalid_b_up(i)=1;invalid_b_down(i)=1;
            Iymax_b(i)=search_l+1;
            Iymax_real_b(i)=Q(i,2);
        end%����ҵ����һ����û�ҵ����������Ǹ�if�ĵ㣬�ǾͲ����ˣ���������������ԭ��λ�ðɡ�
    end
end
% ���£������Ҷ��ݶ����ֵλ���еġ��İ�����
Iymax_real_r=smooth(Iymax_real_r,39,'sgolay',2);%ͼ���˲������ִ��������˶�����������Ե������������˶������
Iymax_real_b=smooth(Iymax_real_b,39,'sgolay',2);
Iymax_r=smooth(Iymax_r,39,'sgolay',2);
Iymax_b=smooth(Iymax_b,39,'sgolay',2);
if sum(abs(red_dot_ref)>0.1)<=0.5*xlen && sum(abs(red_dot_ref)>0.2)<=0.25*xlen ...
        && ~(sum(red_dot_ref>0)>=0.8*xlen||sum(red_dot_ref<0)>=0.8*xlen)
    [Iymax_real_r,red_alert]=repair_bumpupanddown(I,Iymax_real_r,0.1,0.2,30,'red',[]);
else
    red_alert=[0 0 0 0];
end
if sum(abs(blue_dot_ref)>0.1)<=0.5*xlen && sum(abs(blue_dot_ref)>0.2)<=0.25*xlen
    [Iymax_real_b,~]=repair_bumpupanddown(I,Iymax_real_b,0.1,0.2,30,'blue',red_alert);
end
% 3.3 �۲쵱ǰLI/MA���ߵ�б�ʣ��Լ��ƶ��������¼������
red_dot_front_is_neg=red_dot(1:7)<0;%�жϺ���ǰ7����б���Ƿ��Ǹ�����������red_dot_is_neg���7��������1��б���Ǹ�����ֵ�ϴ�ĸ���������˵��ǰ7�������ֵһֱ���½�����matlab��ͼ����λ����������
blue_dot_front_is_neg=blue_dot(1:7)<0;%�ж�����ǰ7����б���Ƿ��Ǹ�����
red_dot_front_is_pos=red_dot(1:7)>0;%�жϺ���ǰ7����б���Ƿ�������������ǣ�˵������ǰ7�������ֵһֱ����������matlab��ͼ����λ���½�����˵���п��ܺ���ǰ�����£���matlab��ͼ�������ϣ����ˡ�
blue_dot_front_is_pos=blue_dot(1:7)>0;%�ж�����ǰ7����б���Ƿ���������
red_dot_tail_is_neg=red_dot(xlen-7:xlen)<0;%�жϺ��ߺ�7����б���Ƿ��Ǹ���������ǣ�˵�����ߺ�7�������ֵһֱ���½�����matlab��ͼ����λ�����������п��������£���matlab��ͼ�������ϣ��̡�
blue_dot_tail_is_neg=blue_dot(xlen-7:xlen)<0;%�ж����ߺ�7����б���Ƿ��Ǹ�����
red_dot_tail_is_pos=red_dot(xlen-7:xlen)>0;%�жϺ��ߺ�7����б���Ƿ�������������ǣ�˵�����ߺ�7�������ֵһֱ����������matlab��ͼ����λ���½������п��������ϣ���matlab��ͼ�������£��̡�
blue_dot_tail_is_pos=blue_dot(xlen-7:xlen)>0;%�ж����ߺ�7����б���Ƿ���������
movement_r=Iymax_real_r-yi_red-1;%�ղ��Ǹ���forѭ��������Ǹ�if���ƣ����ƺ����߲�Ҫ��̫Զ��������������ԣ����Ƕ���1040R�����Ķ����Ƚϴ�ľͲ��У���������������䡣
movement_b=Iymax_real_b-yi_blue-1;%
movement_r_count_down=sum(movement_r>5);%���߻Ҷ��ݶ���� �� ԭ������λ�� �·� �ĵ�ĸ���
movement_r_count_down1=sum(movement_r>0);
movement_r_count_down2=sum(movement_r>7);
movement_b_count_down=sum(movement_b>5);
movement_b_count_down1=sum(movement_b>0);
movement_b_count_down2=sum(movement_b>7);
% 3.4 ��������ĻҶȡ�б�ʡ��˶����ȵ���ֵ��LI/MA�ߵľ��룬ȷ��Ҫ��Ggc��Gdc���ںδ����Ե���Eext��
blue_suspect=0;
if (~eliminate_x_force)%ֻ�Ժ�����ֻ��y�����˶���ʱ�����Ч
    for i=1:xlen
        % ���ȣ����������LI/MAλ�õġ����������ڵġ��ҶȾ�ֵ��
        mean_red_up(i)=mean2(I(yi_red(i)-neighbor_out_r(i):yi_red(i),max(1,i-2):min(i+2,xlen)));% ������yi_red(i)����-1����Ϊȡ��fix��
        mean_red_down(i)=mean2(I(yi_red(i)+1:yi_red(i)+1+neighbor_in(i),max(1,i-2):min(i+2,xlen)));% ������yi_red(i)Ҫ+1����Ϊfix�õ���������LI���Ϸ�����Ϣ��
        mean_blue_up(i)=mean2(I(yi_blue(i)-neighbor_in(i):yi_blue(i),max(1,i-2):min(i+2,xlen)));
        mean_blue_down(i)=mean2(I(yi_blue(i)+1:yi_blue(i)+1+neighbor_out_b(i),max(1,i-2):min(i+2,xlen)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%˼·�����ߺ����߶��������Ĵ���ʽ����һ�������ȣ��������û�����������߱Ƚ�������%%%%%%%%
        %%%%%%�Ͷ����ñȽϴ��margin_r����margin_b��������Ҫ���GgcȥӰ�����ǡ�Ȼ�����LI��%%%%%%%%%
        %%%%%%����������ô������֪��Ӧ���������������˶��������������˶������ԣ�����ֱ�Ӹı���ֵ%%%%%%%
        %%%%%%�������������˶��������������˶���Ȼ����������߲������Ļ������ǲ���ȷ������Ӧ������%%%%%%
        %%%%%%�����˶��������£�����Ҫ��С������ֵ������������ֵ���������������ϣ����������£���%%%%%%%%
        %%%%%%�������ϻ����£����Ҷ��ݶ����ֵ��λ�ã�Iymax_b(i)��search_l+1�ıȽϣ����ж�����%%%%%%%%
        %%%%%%�����ϻ������¡�Ҳ����˵���Ҷ�ֵ�����Ӱ�챻���������Ҷ��ݶ�λ�þ����˶��ķ��� %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3.4.1����Ggc��LI�����˶��������ǣ���������Ϸ��Ҷ�ֵ̫�������ǲ���ֵ������20%�����ͼ�Ggc��
        % ���������У��������������Ҷ��ݶȵȵĸ����ж�������
        if ~ismember(i,exists_noise)
            % ͨ���ж�����������ı����Ggc����ֵ��
            if noise_flag0==1%��������������Ļ�����ô�����������û������Ҳ��margin_r_push_up(i)����һЩ�����������ơ�
                margin_r_push_up(i)=20;
            end%
            if mean_red_up(i)>20%��������Ϸ�ƽ��ֵ���ˣ���һ��ʼ���Ǳ��ж�Ϊû������������ô�ͺܿ��ɣ��п�����һ��ʼû�������������ˡ�
                margin_r_push_up(i)=30;%��ʵ���Ǳ��������ˡ��������mean_red_up(i)�ܴ���30����
            end
            %�����Ҫ�����Լ�һ�£����������������û����������ô�졣
        end
        if(ismember(i,exists_noise) && ~ismember(i,medium_noise) && ~ismember(i,severe_noise))
            push_r_up_thr(i)=grayscale_rmax_down_mean(i);%���ڵ�i���㣬������������down�ĻҶ����жϺ����Ƿ������ơ�
            %margin_r_push_up(i)������ʼֵ10;
            if noise_flag0==1%��������������Ļ�����ô��������margin_r_push_up(i)����һЩ�����������ơ�
                margin_r_push_up(i)=20;
            end
        end
        if(ismember(i,medium_noise))
            push_r_up_thr(i)=grayscale_rmax_down_mean(i);
            %margin_r_push_up(i)������ʼֵ10;
            if noise_flag0==1%��������������Ļ�����ô��������margin_r_push_up(i)����һЩ�����������ơ�
                margin_r_push_up(i)=20;
            end
        end
        if(ismember(i,severe_noise))
            push_r_up_thr(i)=grayscale_rmax_down_mean(i);
            %margin_r_push_up(i)������ʼֵ10;
            if noise_flag0==1%��������������Ļ�����ô��������margin_r_push_up(i)����һЩ�����������ơ�
                margin_r_push_up(i)=20;
            end
        end
        if mean_red_up(i)>=push_r_up_thr(i)+margin_r_push_up(i) && (invalid_r_up(i)==1 || Iymax_r(i)<=search_l && mean_down_r(i)>mean_up_r(i))...  %|| (~ismember(i,exists_noise) && sum(I(yi_red(i),max(1,i-5):min(i+5,xlen))>=15)>=10 && I(yi_red(i),i)>=15 && I(yi_red(i)-1,i)>=15 && I(yi_red(i)-2,i)>=15) ...
                || (~ismember(i,exists_noise) && movement_r(i)<-2 && I(yi_red(i),i)>=15 && I(yi_red(i)-1,i)>=15 && I(yi_red(i)-2,i)>=15) ...
                || grayscale_rmax_up_mean(i)<=20 && min(I(yi_red(i),max(1,i-10):min(i+10,xlen)))>=30 && min(I(yi_red(i)-2:yi_red(i),i))>=30 ...
                || movement_r(i)<-3.5 && I(yi_red(i),i)>=50 && I(yi_red(i)-1,i)>=50 && I(yi_red(i)-2,i)>=50 ...
                || (mean_red_up(i)>=10&&mean_red_down(i)<10&&inside_IMC_is_dark(i)==1)
            % ��һ��������grayscale_rmax_up_mean(i)+margin_r_push_up(i)ȡ��Խ��
            % ��Խ��ʹ���������������������Ҳ����Խ��ͨ��Ggc���������˶���
            Eext_red(yi_red(i):yi_red(i)+7,i)=Eext_red(yi_red(i):yi_red(i)+7,i)+volcano;
            Eext_red(yi_red(i):yi_red(i)+7,max(1,i-1):min(i+1,xlen))=Eext_red(yi_red(i):yi_red(i)+7,max(1,i-1):min(i+1,xlen))+volcano/2;
            Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+1.5*volcano;
        end
        % 3.4.2����Ggc��LI�����˶��������ǣ���������·��Ҷ�ֵ̫С�������ǲ���ֵ������20%���������ơ�
        % ���������У�������ǰ�棬Ҳ�Ǽ������������Ҷ��ݶȵȵĸ����ж�������
        if(~ismember(i,exists_noise))
            if noise_flag0==0%�������û����������ô��������ֵ��һ����С��������Ҫ���ơ�
                push_r_down_thr(i)=5;
                margin_r_push_down(i)=0;%���û�����Ļ�������mean_red_down(i)<push_r_down_thr(i)-margin_r_push_down(i)=10�����ơ�
            end
            if noise_flag0==1%�����������������ֵ�����
                if ~ismember(i,[1:20 xlen-19:xlen])  % �Ƿ�Ϊ��ǰ������󼸸��㣨Ӧ���ǵ����۲��⼸����ĻҶ����ڲ�ͬ��
                    push_r_down_thr(i)=15;
                else
                    push_r_down_thr(i)=10;
                end
                margin_r_push_down(i)=0;%���û�����Ļ�������mean_red_down(i)<push_r_down_thr(i)-margin_r_push_down(i)=10�����ơ�
            end
        end
        if(ismember(i,exists_noise) && ~ismember(i,medium_noise) && ~ismember(i,severe_noise))%
            if grayscale_rmax_down_mean(i)>=30
                push_r_down_thr(i)=grayscale_rmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_rmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_rmax_down_mean
                %����Ӧ�ڣ�����Ϸ�̫������������������������ơ�Ȼ���أ���������Ʋ�һ�����ǣ�����Ϸ���������Ҫ�������������ƣ�����margin_r_push_down���ڱ�ġ���
            end
            margin_r_push_down(i)=15;%�·�̫��������30���ģ�����ԭ���Ϸ��ļ�15���������ԭ���·��ģ�С��30�ˣ�����15�������15���¡���Ϊ���������ֵ������15�����������ǣ�����·�С��15�����ˣ��Ǿ����ưɡ�
            %���������invalid_r(i)Ԥ��Ϊ0��
            %����Ԥ��һ�£����noise_flag0,noise_flag1==1������ô��.
        end
        if(ismember(i,medium_noise))
            if grayscale_rmax_down_mean(i)>=30
                push_r_down_thr(i)=grayscale_rmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_rmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_rmax_down_mean
            end
            margin_r_push_down(i)=10;%������е���������ô��������ֵ�ĳ�10�������·�̫���ģ������Ϸ��ļ�10���������·��ļ�10�����20���£���������һЩ��
            %���������invalid_r(i)Ԥ��Ϊ0��
            %����Ԥ��һ�£����noise_flag0,noise_flag1==1������ô��.
        end
        if(ismember(i,severe_noise))
            if grayscale_rmax_down_mean(i)>=30
                push_r_down_thr(i)=grayscale_rmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_rmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_rmax_down_mean
            end
            margin_r_push_down(i)=5;%��������
            %���������invalid_r(i)Ԥ��Ϊ0��
            %����Ԥ��һ�£����noise_flag0,noise_flag1==1������ô��.
        end
        if mean_red_down(i)<=max(5, push_r_down_thr(i)-margin_r_push_down(i)) &&~(mean_red_down(i)<=15&&mean_red_up(i)<=15&&max(J(yi_red(i):yi_blue(i)-fix(IMW(i)/2),i))<3)... %
                || mean_red_down(i)<=grayscale_rmax_down_mean(i) && movement_r_count_down1>xlen*(300/361) && movement_r(i)>2 && mean(movement_r)>1.5 ...%170525��300�ĳ���xlen*(300/361) %&& ~(min(I(yi_red(i),max(1,i-20):min(i+20,xlen)))>20)
                && (invalid_r_down(i)==1 || Iymax_r(i)>=search_l && mean_down_r(i)>mean_up_r(i))...
                || (I(yi_red(i),i)<=5 && I(yi_red(i)+1,i)<=5 && I(yi_red(i)+2,i)<=5) ...
                && inside_IMC_is_dark(i)==0
        %���������˶��ĸ��������ܽᡣ
            Eext_red(yi_red(i)-7:yi_red(i),i)=Eext_red(yi_red(i)-7:yi_red(i),i)+volcano;
            Eext_red(yi_red(i)-7:yi_red(i),max(1,i-1):min(i+1,xlen))=Eext_red(yi_red(i)-7:yi_red(i),max(1,i-1):min(i+1,xlen))+volcano/2; 
            Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+1.5*volcano;
        end
        %3.4.3����Ggc��MA�����˶���ÿ���㶼ȷ���Ƿ�������Խ��������Խû��Ҫ������ӵ�����
        if(~ismember(i,exists_fuzzy))
            push_b_up_thr(i)=max(210,grayscale_bmax_up_mean(i));%ԭ��˵��Ҳ�ǣ�ֻҪ�����������Ͳ�Ҫ�������˾͵��ˡ�
            %margin_b_push_up(i)������ʼֵ����50;
        end
        if(ismember(i,exists_fuzzy) && ~ismember(i,medium_fuzzy) &&~ismember(i,severe_fuzzy))
            %�������Ĳ��֣�������ֵ������С�����������������˶���
            margin_b_push_up(i)=40;
        end
        if(ismember(i,medium_fuzzy))
            %push_b_up_thr(i)������ʼֵ����grayscale_bmax_up_mean(i)
            margin_b_push_up(i)=35;
        end
        if(ismember(i,severe_fuzzy))
            %push_b_up_thr(i)������ʼֵ����grayscale_bmax_up_mean(i)
            margin_b_push_up(i)=30;
        end
        if mean_blue_up(i)>=push_b_up_thr(i)+margin_b_push_up(i) && (invalid_b_up(i)==1 || Iymax_b(i)<=search_l && mean_down_b(i)>mean_up_b(i))
        %��������Ϸ��������ص�ƽ��ֵ�е���ڵİף�����push_b_up_thr(i)+margin_b��������Ϊ�����ܵ��˰�ɫ����֯���档
        %push_b_up_thr(i)+margin_b_push_up(i)ȡ��Խ�󣬾�Խ�������ơ�
        %Iymax_b(i)<search_l+1��ʾ���Ҷ��ݶ����ֵ��Ҫ��ԭ����������λ���Ϸ���<search_l+1�����������ơ�mean_down_b(i)>mean_up_b(i)�����͡�
            Eext_blue(yi_blue(i):yi_blue(i)+7,i)=Eext_blue(yi_blue(i):yi_blue(i)+7,i)+volcano;
            Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))=Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))+volcano/2; %�ݼ���Χ
            Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+1.5*volcano;
        end
        if mean_blue_up(i)>mean_blue_down(i) && ~(max(I(yi_blue(i)+1:yi_blue(i)+10,i))>=max(I(yi_blue(i)-neighbor_in(i):yi_blue(i),i))) ...%20170502
                && ~(grayscale_bmax_up_mean(i)>grayscale_bmax_down_mean(i)) ...
                && ~(sum(grayscale_bmax_up_mean(max(1,i-50):min(xlen,i+50))>150)>10) ...
                && inside_IMC_is_bright(i)==0 ...
                && max(abs(blue_dot(max(1,i-20):min(i+20,xlen))))>=0.05 ... 
                &&~(sum(blue_dot(max(1,i-20):i)<0)>10&&sum(blue_dot(i:min(i+20,xlen))>0)>10)
            Eext_blue(yi_blue(i):yi_blue(i)+7,i)=Eext_blue(yi_blue(i):yi_blue(i)+7,i)+2*volcano;
            Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))=Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-1):min(i+1,xlen))+1.5*volcano; %�ݼ���Χ
            Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+4*volcano;
        end
        %3.4.4����Ggc��MA�����˶������������棬ͨ���Ҷ�б�ʵ��������ж���ֵ�ı仯��
        if(~ismember(i,exists_fuzzy))
            if grayscale_bmax_down_mean(i)>=180%Ӧ����һ��Ѵ���180�ģ���������·�Ҳ��һ��Ѵ���30�ġ�
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_bmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_bmax_down_mean
                %�ͺ������ƺ���
                margin_b_push_down(i)=10;%��������Ļ�����margin_b_push_down(i)�ĳ�10���Ͼ��Ѿ�������ĻҶ�ȥ���ˡ�����
            end
            if mean_blue_down(i)<=110 && std2(I(yi_blue(i)-neighbor_in(i):yi_blue(i),max(1,i-2):min(i+2,xlen)))<=20%�������һ��ʼ�������ģ��������ڲ������ˣ��ǾͱȽϿ��ɡ������ں������ƣ�ֻ��������˸��Ҷ�ֵ�Ƚ�С��������
                blue_suspect=1;
            end%������߿��ɣ�����ͼ�����߲�������2����
        end
        if(ismember(i,exists_fuzzy) && ~ismember(i,medium_fuzzy) &&~ismember(i,severe_fuzzy))%�������Ĳ��֣�������ֵ������С��
            if grayscale_bmax_down_mean(i)>=180
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_bmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_bmax_down_mean
                %�ͺ������ƺ���
                margin_b_push_down(i)=10;%��������Ļ�����margin_b_push_down(i)�ĳ�10���Ͼ��Ѿ�������ĻҶ�ȥ���ˡ�����
            end
        end
        if(ismember(i,medium_fuzzy))
            if grayscale_bmax_down_mean(i)>=180
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_bmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_bmax_down_mean
                %�ͺ������ƺ���
                margin_b_push_down(i)=10;%��������Ļ�����margin_b_push_down(i)�ĳ�10���Ͼ��Ѿ�������ĻҶ�ȥ���ˡ�������ʵ10��δ���У�Ū���õ�-10������
            end
        end
        if(ismember(i,severe_fuzzy))
            if grayscale_bmax_down_mean(i)>=180
                push_b_down_thr(i)=grayscale_bmax_up_mean(i);
                %����·�̫�����Ǿ���grayscale_bmax_up_meanȥ�ƣ�������Ĭ�ϵ�grayscale_bmax_down_mean
                %�ͺ������ƺ���
                margin_b_push_down(i)=10;%��������Ļ�����margin_b_push_down(i)�ĳ�10���Ͼ��Ѿ�������ĻҶ�ȥ���ˡ�����
            end
        end
        if (mean_blue_down(i)<=push_b_down_thr(i)-margin_b_push_down(i) || mean_blue_down(i)<=grayscale_bmax_down_mean(i) && movement_b_count_down1>xlen*(300/361) && movement_b(i)>2 && mean(movement_b)>1.5)...
                && (invalid_b_down(i)==1 || Iymax_b(i)>=search_l && mean_down_b(i)>mean_up_b(i))
            Eext_blue(yi_blue(i)-7:yi_blue(i),i)=Eext_blue(yi_blue(i)-7:yi_blue(i),i)+volcano;
            Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-1):min(i+1,xlen))=Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-1):min(i+1,xlen))+volcano/2;
            Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+1.5*volcano;
        end
        % ���ϣ��ûҶ��ж��Ƿ����Ggc�Ĳ��ֽ��������棬�����������������������
        % 3.4.5���LI/MAֻ����2���ص����ڣ��Ϳ��������⡣���ݾֲ��Ҷȸ���һ�£���������⣬�ͼ���Ggc������λ�á�
        tol=min(IMW(i),6);
        if (Q(i,2)-P(i,2)<=tol)
            yi_mid=fix((yi_red(i)+yi_blue(i))/2);
            I_local=I(yi_mid-2:yi_mid+4,max(P(i,1)-2,1):min(P(i,1)+2,xlen));%����Ϊ�˷�ֹ��������0���ߴ���xlen�������
            mean_local=mean(mean(I_local));%�ֲ�ƽ���Ҷ�
            volcano1=200*(tol/(Q(i,2)-P(i,2)));
            if abs(mean_local-grayscale_rmax_down_mean(i))<10 || abs(mean_local-(grayscale_bmax_up_mean(i)/2+grayscale_bmax_down_mean(i)/2))>10 && mean_blue_down(i)<=grayscale_bmax_down_mean(i)
                %��if�����������ϵ����
                if (ismember(i,1:7) && min(abs(blue_dot_front_is_pos))==1) || ...
                        (ismember(i,xlen-7:xlen) && min(abs(blue_dot_tail_is_neg))==1)
                    volcano1=volcano1*2;%�����������߻����ұ������̣��Ӵ�Ggc���������ߡ�
                end
                if Q(i,2)-P(i,2)<=tol*0.8%���̫С�ˣ��Ǿ�ֱ������������������
                    Q(i,2)=Q(i,2)+tol/2;
                    disp('blue pushed up manually because it is getting too close to red')
                    disp(fr)
                end
                Eext_blue(yi_blue(i)-7:yi_blue(i),i)=Eext_blue(yi_blue(i)-7:yi_blue(i),i)+2*volcano1;
                Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+4*volcano1;%��β�Ҫ�ݼ���Χ�ˣ���Ϊ��һ���㻹��ѭ����������ġ���Ҫ�ǻݼ���Χ�ˣ��൱�ڼ��˺ܶ�λ�ɽ��
            end
            if abs(mean_local-grayscale_bmax_up_mean(i))<10 || abs(mean_local-(grayscale_rmax_up_mean(i)/2+grayscale_rmax_down_mean(i)/2))>10 && mean_red_up(i)>max(10,grayscale_rmax_up_mean(i))
                %��if����������µ����
                if (ismember(i,1:7) && min(abs(red_dot_front_is_neg))==1) || ...
                        (ismember(i,xlen-7:xlen) && min(abs(red_dot_tail_is_pos))==1)
                    volcano1=volcano1*2;%����Ǻ�����߻����ұ�������������ģ��Ӵ��ɽ�����ơ�����tmd�β������ˡ�
                end
                if Q(i,2)-P(i,2)<=tol*0.8%���̫С�ˣ��Ǿ�ֱ�����Ѻ���������.����������һ���ټӻ�ɽ�ɣ������̫���ˣ��³��¶���
                    P(i,2)=P(i,2)-tol/2;
                    disp('red pushed up manually because it is getting too close to blue')
                    disp(fr)
                end
                if P(i,2)-yi_red(i)>0.8, trick=1;else trick=0;end%�������������***.9�Ļ�����ô�������������һ�����ٿ�ʼ�ӻ�ɽ��������ܷ������������ˡ�1060R��88�š�
                %���߲��ü������������Ϊ�������ƣ�������ӻ�ɽ�������Ǹ�yi_blue�Ͱ�С����Ĩ���ˡ�
                Eext_red(yi_red(i)+trick:yi_red(i)+trick+7,i)=Eext_red(yi_red(i)+trick:yi_red(i)+trick+7,i)+2*volcano1;
                Eext_red(yi_red(i)+trick,i)=Eext_red(yi_red(i)+trick,i)+4*volcano1;
            end
        end
        % 3.4.6����İ�/�����������
        % �ú����ͼ�͵�һ��ͼÿһ���б��ȥ�Աȣ������̫���ˣ��Ǿͼ�Gdc���ϡ�ʵ������ʱ�򣬵�һ��ͼ��б�ʲ���0����1�������õڶ���ͼ����Ϊ�ж���׼��
        if fr~=1
            vol_times=10;%��Ҫ���Ǹ�б��ƽ��ֵ̫С�ˣ�һ����0.������0.0������������Ӵ�һЩGdc�ı�����
            mean5_abs_red_dot_ref=mean(abs(red_dot_ref(max(1,i-5):min(i+5,xlen))));%��i�㸽���ġ���2��ͼ�ĺ�����б�ʾ���ֵ�ľ�ֵ��
            mean5_abs_blue_dot_ref=mean(abs(blue_dot_ref(max(1,i-5):min(i+5,xlen))));
            mean10_abs_red_dot_ref=mean(abs(red_dot_ref(max(1,i-10):min(i+10,xlen))));%��i�㸽���ġ���2��ͼ�ĺ�����б�ʾ���ֵ�ľ�ֵ��
            mean10_abs_blue_dot_ref=mean(abs(blue_dot_ref(max(1,i-10):min(i+10,xlen))));
            if (i>7 && i<xlen-7)%ǰ7����ͺ�7�������Ҳ��ܣ����浥�����ǡ�
                if (red_dot(i)<-0.01 && red_dot(i+1)>0.01)%LI�߳��ּ�Сֵ��ͼ�б�ʾ���ϵĹİ�����ע�⣬��ʱ����ֱ�Ӿ͸�����Gdc�������ж�б�ʾ���ֵ������һ��if��
                    %������0����Ϊ�˷�ֹ��Щ1*10^-4������ֵ������������������0ȫ���ĳ�0.01�ˡ�
                    mean5_abs_red_dot=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_red_dot=mean(abs(red_dot(max(1,i-10):min(i+10,xlen))));
                    if (mean5_abs_red_dot>mean5_abs_red_dot_ref+0.05 && sum(Iymax_real_r(i-5:i+5)>P(i-5:i+5,2))>=8) ...
                            ||(mean10_abs_red_dot>mean10_abs_red_dot_ref+0.05 && sum(Iymax_real_r(max(1,i-10):min(i+10,xlen))>P(max(1,i-10):min(i+10,xlen),2))>=16) ...
                            || sum(Iymax_real_r(i-5:i+5)>P(i-5:i+5,2)+7)>=8
                        % �ùİ�������б�ʣ��������ֲ�ͬ�����ڵ�б�ʣ����Ҷ��ݶ����ֵ�㣬��Щ���ص��ۺϿ��ǣ��ж��Ƿ�����Gdc��
                        % ��Щ�ֹ����ԣ���Ҫ�����ܷ����Զ��������档
                        for try_loca=1:30%���forѭ�������ҹİ��ĳ���
                            if mean(abs(red_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(red_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;%i-forced_area:i+forced_area�����Χ�ڣ���Gdc��
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))-red_dot_ref(max(1,i-5):min(i+5,xlen))));
                        volcano1=volcano*vol_times*scalarco;%��Gdc�����ƽ���ҶȾ���ֵԽ��Gdc��Խ��
                        Eext_red(yi_red(i)-7:yi_red(i),max(1,i-forced_area):min(i+forced_area,xlen))=Eext_red(yi_red(i)-7:yi_red(i),max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
                if (red_dot(i)>0.01 && red_dot(i+1)<-0.01)%LI�߳��ּ���ֵ��ͼ�б�ʾ���µĹİ���
                    mean5_abs_red_dot=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_red_dot=mean(abs(red_dot(max(1,i-10):min(i+10,xlen))));
                    if (mean5_abs_red_dot>mean5_abs_red_dot_ref+0.05 && sum(Iymax_real_r(i-5:i+5)<P(i-5:i+5,2))>=8) ...
                            ||(mean10_abs_red_dot>mean10_abs_red_dot_ref+0.05 && sum(Iymax_real_r(max(1,i-10):min(i+10,xlen))<P(max(1,i-10):min(i+10,xlen),2))>=8)...
                            || sum(Iymax_real_r(i-5:i+5)<P(i-5:i+5,2)-7)>=8 && lightflag==0
                        %�߼�������һ����
                        for try_loca=1:30
                            if mean(abs(red_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(red_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;%i-forced_area:i+forced_area�����Χ�ڣ���Gdc��
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(red_dot(max(1,i-5):min(i+5,xlen))-red_dot_ref(max(1,i-5):min(i+5,xlen))));
                        volcano1=volcano*vol_times*scalarco;
                        Eext_red(yi_red(i):yi_red(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))=Eext_red(yi_red(i):yi_red(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
                if (blue_dot(i)<-0.01 && blue_dot(i+1)>0.01)%MA�߳��ּ�Сֵ��ͼ�б�ʾ���ϵĹİ�����
                    mean5_abs_blue_dot=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_blue_dot=mean(abs(blue_dot(max(1,i-10):min(i+10,xlen))));
                    if lightflag==1, mean5_abs_blue_dot=mean5_abs_blue_dot+0.1;mean10_abs_blue_dot=mean10_abs_blue_dot+0.1;end%�Ҿ�������1055L1�����ϵĹİ����͸�������ȥ��������������ã���Ϊ�е��ƶ��ˡ��Ȳ����ˣ�������
                    if (mean5_abs_blue_dot>mean5_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(i-5:i+5)>Q(i-5:i+5,2))>=8) ...
                            || sum(Iymax_real_b(i-5:i+5)>Q(i-5:i+5,2)+7)>=8 ...
                            || (mean10_abs_blue_dot>mean10_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(max(1,i-10):min(i+10,xlen))>Q(max(1,i-10):min(i+10,xlen),2))>=16)
                        %ͬ��,�ùİ�������б�ʣ��������ֲ�ͬ�����ڵ�б�ʣ����Ҷ��ݶ����ֵ�㣬��Щ���ص��ۺϿ��ǣ��ж��Ƿ�����Gdc��
                        for try_loca=1:30%Ӧ�ò����г���30�����ص�Ĺİ��ˡ�����
                            if mean(abs(blue_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(blue_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;%i-forced_area:i+forced_area�����Χ�ڣ���Gdc��
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))-blue_dot_ref(max(1,i-5):min(i+5,xlen))));
                        if lightflag==1, vol_times=vol_times*4;end
                        volcano1=volcano*vol_times*scalarco;
                        Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-forced_area):min(i+forced_area,xlen))=Eext_blue(yi_blue(i)-7:yi_blue(i),max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
                if (blue_dot(i)>0.01 && blue_dot(i+1)<-0.01)%MA�߳��ּ���ֵ��ͼ�б�ʾ���µĹİ�����
                    mean5_abs_blue_dot=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))));
                    mean10_abs_blue_dot=mean(abs(blue_dot(max(1,i-10):min(i+10,xlen))));
                    if (mean5_abs_blue_dot>mean5_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(i-5:i+5)<Q(i-5:i+5,2))>=8)  ...
                            ||(mean10_abs_blue_dot>mean10_abs_blue_dot_ref+0.05 && sum(Iymax_real_b(max(1,i-10):min(i+10,xlen))<Q(max(1,i-10):min(i+10,xlen),2))>=16)...
                            || sum(Iymax_real_b(i-5:i+5)<Q(i-5:i+5,2)-7)>=8 && lightflag==0
                    % �ùİ�������б�ʣ��������ֲ�ͬ�����ڵ�б�ʣ����Ҷ��ݶ����ֵ�㣬��Щ���ص��ۺϿ��ǣ��ж��Ƿ�����Gdc��
                        for try_loca=1:30
                            if mean(abs(blue_dot(max(1,i-try_loca):min(i+try_loca,xlen))))<mean(abs(blue_dot(max(1,i-(try_loca-1)):min(xlen,i+(try_loca-1)))))
                                forced_area=try_loca;
                                break;
                            end
                            forced_area=try_loca;
                        end
                        scalarco=mean(abs(blue_dot(max(1,i-5):min(i+5,xlen))-blue_dot_ref(max(1,i-5):min(i+5,xlen))));
                        volcano1=volcano*vol_times*scalarco;
                        Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))=Eext_blue(yi_blue(i):yi_blue(i)+7,max(1,i-forced_area):min(i+forced_area,xlen))+volcano1;
                    end
                end
            end
            v1=5;
            % ���������ƵĴ���İ��ĳ���ֻ�����ǿ���ǰ��7��������7����ġ�
            if ismember(i,1:7)
                m_rd_f=mean(red_dot(1:7));m_rd_f_r=mean(red_dot_ref(1:7));
                m_bd_f=mean(blue_dot(1:7));m_bd_f_r=mean(blue_dot_ref(1:7));
                if min(abs(red_dot_front_is_neg))==1 && abs(m_rd_f)>0.1 && ~(m_rd_f_r<-0.1) || (m_rd_f_r<-0.2 && abs(m_rd_f)>=abs(m_rd_f_r)*1.5)
                    %�Ǹ�if�����ķ��������ǰ7�������б�ʶ��Ǹ��ģ��Ҿ���ֵƽ��ֵ����0.1����˵������ǰ��7����������ϣ���matlab��ͼ�������£��̣���Ҫ��Gdc������
                    %���ǣ�����ڶ���ͼǰ7����б��ƽ��ֵҲС��-0.1������б�ʶ��Ǹ��ģ��Ҿ���ֵ����0.1������ô������Ϊ���������ˡ�
                    volcano1=volcano*v1*max(0,mean(red_dot_ref(1:7)-red_dot(1:7)));%�����һ��ͼ�Ƚ�ƽ����б��û����ô�ࣩ��
                    if sum(Iymax_real_r(1:7)<P(1:7,2)-5)==8,volcano1=volcano1*2;end%������֣�ǰ8������ÿһ����ĺ����ݶ����ֵ������ԭ�������Ϸ�5���ϣ���ô��Gdc�ӱ���
                    if movement_r_count_down1<=xlen*(180/361) %�������һ��ĵ��λ��С��0���������ˣ���
                        %���if��Ϊ�˷�ֹ�е�ʱ��LI/MA�����������ƣ�����Ҳ���㡰min(abs(red_dot_front_is_neg))==1...�����if������
                        %�����¼��������Gdc�����������Gdc����������ʹ�����ڳ��������������¿�ʼ��ʼ�Ҷ�����
                        %movement_r_count1�ǡ�LI���ݶ����ֵ������LI���·��ĵ���������������С��һ�룬��֤��������ͼ�У�
                        %LI������û�з����������µ��ƶ����ŷ��ĵؼ�Gdc��
                        Eext_red(yi_red(i):yi_red(i)+4,i)=Eext_red(yi_red(i):yi_red(i)+4,i)+volcano1;%�·���Gdc��
                        Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                    end
                end
                if min(abs(blue_dot_front_is_neg))==1 && abs(m_bd_f)>0.1 && ~(m_bd_f_r<-0.1) || (m_bd_f_r<-0.2 && abs(m_bd_f)>=abs(m_bd_f_r)*1.5)%�����������˵������ǰ��7�������ϣ���matlab��ͼ�������£�����
                    volcano1=volcano*v1*max(0,mean(blue_dot_ref(1:7)-blue_dot(1:7)));%
                    if sum(Iymax_real_b(1:7)<Q(1:7,2)-5)==8,volcano1=volcano1*2;end%
                    if movement_b_count_down1<=xlen*(180/361) 
                        %���������ƣ��������һ��ĵ��λ��С��0���������ˣ�����MA��û�����Ե������ƶ����ŷ��ĵؼ�Gdc��
                        Eext_blue(yi_blue(i):yi_blue(i)+4,i)=Eext_blue(yi_blue(i):yi_blue(i)+4,i)+volcano1;
                        Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                    end
                end
                if min(abs(red_dot_front_is_pos))==1 && abs(m_rd_f)>0.1 && ~(m_rd_f_r>0.1) || (m_rd_f_r>0.2 && abs(m_rd_f)>=abs(m_rd_f_r)*1.5)%�����������˵������ǰ��7�������£���matlab��ͼ�������ϣ�����
                    volcano1=volcano*v1*max(0,mean(red_dot(1:7)-red_dot_ref(1:7)));%���ʱ��Ӧ�÷���������
                    if sum(Iymax_real_r(1:7)>P(1:7,2)+5)==8,volcano1=volcano1*2;end%
                    %������ټ�if movement_r_count1<=180�����������ˣ���Ҫ��û����LI/MA�����Կ��ٵ������Ƶġ�
                    Eext_red(yi_red(i)-4:yi_red(i),i)=Eext_red(yi_red(i)-4:yi_red(i),i)+volcano1;
                    Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                end
                if min(abs(blue_dot_front_is_pos))==1 && abs(m_bd_f)>0.1 && ~(m_bd_f_r>0.1) || (m_bd_f_r>0.2 && abs(m_bd_f)>=abs(m_bd_f_r)*1.5)%�����������˵������ǰ��7�������£���matlab��ͼ�������ϣ�����
                    volcano1=volcano*v1*max(0,mean(blue_dot(1:7)-blue_dot_ref(1:7)));%���ʱ��Ӧ�÷���������
                    if sum(Iymax_real_b(1:7)>Q(1:7,2)+5)==8,volcano1=volcano1*2;end%
                    %������ټ�if movement_b_count1<=180�����������ˣ���Ҫ��û����LI/MA�����Կ��ٵ������Ƶġ�
                    Eext_blue(yi_blue(i)-4:yi_blue(i),i)=Eext_blue(yi_blue(i)-4:yi_blue(i),i)+volcano1;
                    Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                end
            end
            if ismember(i,xlen-7:xlen)  % ��������7��������������ǰ7��������һ���ġ�
                m_rd_l=mean(red_dot(xlen-7:xlen));m_rd_l_r=mean(red_dot_ref(xlen-7:xlen));
                m_bd_l=mean(blue_dot(xlen-7:xlen));m_bd_l_r=mean(blue_dot_ref(xlen-7:xlen));
                if min(abs(red_dot_tail_is_neg))==1 && abs(m_rd_l)>0.1 && ~(m_rd_l_r<-0.1) || (m_rd_l_r<-0.2 && abs(m_rd_l)>=abs(m_rd_l_r)*1.5)%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %������㣺�����7����б�ʶ�С��0���Ҿ���ֵƽ��ֵ����0.1���ڵڶ���ͼ��7����б��ƽ��ֵ����С��-0.1��
                    %���ߣ���һ��ͼ��7����б��ƽ��ֵ�����ڵڶ���ͼ��7����б��ƽ��ֵ��1.5�����ҵڶ���ͼ��7����б��ƽ��ֵС��-0.2���Ƚ�б��
                    %��˵��LI�ߺ���7������matlab��ͼ�������̣�����ƫС����
                    volcano1=volcano*v1*max(0,mean(red_dot_ref(xlen-7:xlen)-red_dot(xlen-7:xlen)));%����Gdcֵ
                    if sum(Iymax_real_r(xlen-7:xlen)>P(xlen-7:xlen,2)+5)==8,volcano1=volcano1*2;end
                    %������֣����8������ÿһ����ĺ����ݶ����ֵ������ԭ�������·�5���£���ô��Gdc�ӱ���
                    Eext_red(yi_red(i)-4:yi_red(i),i)=Eext_red(yi_red(i)-4:yi_red(i),i)+volcano1;
                    Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                end
                if min(abs(blue_dot_tail_is_neg))==1 && abs(m_bd_l)>0.1 && ~(m_bd_l_r<-0.1) || (m_bd_l_r<-0.2 && abs(m_bd_l)>=abs(m_bd_l_r)*1.5)
                    %�����������˵�����ߺ���7������matlab��ͼ�������̣�����ƫС����
                    volcano1=volcano*v1*max(0,mean(blue_dot_ref(xlen-7:xlen)-blue_dot(xlen-7:xlen)));
                    if sum(Iymax_real_b(xlen-7:xlen)>Q(xlen-7:xlen,2)+5)==8,volcano1=volcano1*2;end
                    Eext_blue(yi_blue(i)-4:yi_blue(i),i)=Eext_blue(yi_blue(i)-4:yi_blue(i),i)+volcano1;
                    Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                end
                if min(abs(red_dot_tail_is_pos))==1 && abs(m_rd_l)>0.1 && ~(m_rd_l_r>0.1) || (m_rd_l_r>0.2 && abs(m_rd_l)>=abs(m_rd_l_r)*1.5)
                    volcano1=volcano*v1*max(0,mean(red_dot(xlen-7:xlen)-red_dot_ref(xlen-7:xlen)));
                    if sum(Iymax_real_r(xlen-7:xlen)<P(xlen-7:xlen,2)-5)==8,volcano1=volcano1*2;end%
                    if movement_r_count_down1<=xlen*(180/361) %ȷ��û����LI�����������Ƶ����
                        Eext_red(yi_red(i):yi_red(i)+4,i)=Eext_red(yi_red(i):yi_red(i)+4,i)+volcano1;
                        Eext_red(yi_red(i),i)=Eext_red(yi_red(i),i)+2*volcano1;
                    end
                end
                if min(abs(blue_dot_tail_is_pos))==1 && abs(m_bd_l)>0.1 && ~(m_bd_l_r>0.1) || m_bd_l_r>0.2 && abs(m_bd_l)>=abs(m_bd_l_r)*1.5
                    %�����������˵�����ߺ���7�������ϣ���matlab��ͼ�������£�����
                    volcano1=volcano*v1*max(0,mean(blue_dot(xlen-7:xlen)-blue_dot_ref(xlen-7:xlen)));%����Gdcֵ
                    if sum(Iymax_real_b(xlen-7:xlen)<Q(xlen-7:xlen,2)-5)==8,volcano1=volcano1*2;end%
                    if movement_b_count_down1<=xlen*(180/361) %ȷ��û����MA�����������Ƶ����
                        Eext_blue(yi_blue(i):yi_blue(i)+4,i)=Eext_blue(yi_blue(i):yi_blue(i)+4,i)+volcano1;
                        Eext_blue(yi_blue(i),i)=Eext_blue(yi_blue(i),i)+2*volcano1;
                    end
                end
            end
        end
    end
    % 3.4.7 �ܵ��������������������˶�����
    fl_r=0;fl_b=0;r=0;%ÿһ��ͼ�����㣬�����һ��ͼ�����󣬲���һ��
    % ���Ǵ���LI�����ƶ��������
    if movement_r_count_down>xlen*(180/361) && movement_r_count_down2>xlen*(100/361)  
        % �����180��100Ҳ���Գ����ġ����������180���������ƶ��ˣ�5�����ص㣩������100�������ƶ��ϴ�7�����ص㣩���ż���Ggc��
        P(:,2)=P(:,2)+abs(mean(movement_r));%����һ��˼·��ֱ����������
        for i=1:xlen%%
            Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)=Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)+volcano;
            Eext_red(fix(P(i,2)),i)=Eext_red(fix(P(i,2)),i)+volcano;
        end
        disp('red moves down fast'),disp(fr)
        fl_r=1;fast_r=1;
    else
        if movement_r_count_down1>xlen*(330/361) && fast_r==0 && abs(mean(movement_r))>3 && ...
                size(severe_noise,2)<xlen*(80/361) || ...
                (noise_flag1==1 && movement_b_count_down1>xlen*(330/361) &&...
                fast_b==0 && abs(mean(movement_b))>3 && abs(mean(movement_r))>1)
            %���û���������¶�����ô���ң��������¶��ĵ����࣬Ҳ��Ggc���������ߡ�
            for i=1:xlen  % ������Ը���ϸ����ͬ�ĵ㣬�����˶�״̬����ͬ��ȥ��Ggc����������ûʱ�����ˡ�
                if (~ismember(i,exists_noise) && I(yi_red(i),i)<10) || (ismember(i,exists_noise))
                    %���û��������LI�����Ҷ�С���ż�Ggc����Ҫ�Ǿ���LI�������������ߵ��ǲ����ڡ�
                    if movement_r(i)>5, P(i,2)=P(i,2)+3;disp('red pushed down manually'),disp(fr),end
                    Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)=Eext_red(fix(P(i,2))-4:fix(P(i,2)),i)+volcano*max(0,abs(mean(movement_r))-1)/8;
                    Eext_red(fix(P(i,2)),i)=Eext_red(fix(P(i,2)),i)+volcano*max(0,abs(mean(movement_r))-1)/4;
                end
            end
            disp('red moves down'),disp(abs(mean(movement_r))),disp(fr)
            r=1;
        end
    end
    % ����������˼·������MA����ֵ�᲻̫һ����
    if movement_b_count_down>xlen*(180/361)  && movement_b_count_down2>xlen*(50/361)%
        Q(:,2)=Q(:,2)+abs(mean(movement_b));%
        for i=1:xlen
            Eext_blue(max(1,fix(Q(i,2))-4):fix(Q(i,2)),i)=Eext_blue(max(1,fix(Q(i,2))-4):fix(Q(i,2)),i)+volcano;
            Eext_blue(fix(Q(i,2)),i)=Eext_blue(fix(Q(i,2)),i)+volcano;
        end
        disp('blue moves down fast'),disp(fr)
        fl_b=1;fast_b=1;
    else
        if movement_b_count_down1>xlen*(330/361) && fast_b==0 && abs(mean(movement_b))>3 || (r==1 && abs(mean(movement_b))>1)
            %MAû��LI��ô�������½��ܶ࣬��ʵ������û�������ǿ������������ԣ������LI���ж���׼��̫һ������������֮һ��
            %��ƽ���½���3������ƽ���½���1��������LI�����½���red moves down����r==1
            for i=1:xlen
                if movement_b(i)>6 && grayscale_bmax_down_mean(i)>200, Q(i,2)=Q(i,2)+3;disp('blue pushed down manually'),disp(fr),end%
                %���Ƶأ��ҵ���Ggc��������
                Eext_blue(fix(Q(i,2))-4:fix(Q(i,2)),i)=Eext_blue(fix(Q(i,2))-4:fix(Q(i,2)),i)+volcano*max(0,abs(mean(movement_b))-1)/8;
                Eext_blue(fix(Q(i,2)),i)=Eext_blue(fix(Q(i,2)),i)+volcano*max(0,abs(mean(movement_b))-1)/4;
            end
            disp('blue moves down'),disp(abs(mean(movement_b))),disp(fr)
        end
    end
    % 3.4.8 ���������˶����ƣ��趨���㷨���н�������
    % ����LI
    devi_r=abs(mean(movement_r));
    if fl_r==1%
        Gamma_red=Options.Gamma*min(10,max(2,devi_r));%��
    else
        Gamma_red=Options.Gamma;
    end
    if (noise_flag1==1 && fl_r==0) || (noise_flag2==1 && fl_r==0)%noise_flag1�Ƕ��������ж��Եġ�
        Gamma_red=Options.Gamma/2;
    end
    % MA�ߡ�
    devi_b=abs(mean(movement_b));
    if fl_b==1%
        Gamma_blue=Options.Gamma*min(10,max(2,devi_b));
    else
        Gamma_blue=Options.Gamma;
    end
    if lightflag==1, Gamma_blue=Options.Gamma/4;end%
    if (blue_suspect==1 && lightflag==0 && fl_b==0),Gamma_blue=Options.Gamma/3;end%
    if (fuzzy_flag1==1 || lightflag1==1 && blue_suspect==0 && lightflag==0 && fl_r==0)%
        Gamma_blue=Gamma_blue/2;
    end
    % 4 �������������������Eext����ʼ�ݻ��ߵ�������%
    % 4.1 ����ͼ�������������Լ�LI��MA�����ߵ��ⲿ����������4ʽ��Gext�ĵ�������
    f0=1;%���������Ǹ��������ӡ�
    Fx_red=ImageDerivatives2D(Eext_red,Options.Sigma2,'x');  % LI�ߵ�ͼ�����������е�5ʽ������������ͼ��������һ�����ӡ�
    Fext_red(:,:,1)=-f0*Fx_red*2*Options.Sigma2^2;
    Fext_red(:,:,2)=0;
    Fx_blue=ImageDerivatives2D(Eext_blue,Options.Sigma2,'x');  % MA�ߵ�ͼ����
    Fext_blue(:,:,1)=-f0*Fx_blue*2*Options.Sigma2^2;
    Fext_blue(:,:,2)=0;
    if(Options.Verbose)
        [x,y]=ndgrid(1:10:size(Fext_red,1),1:10:size(Fext_red,2));
        figure,imshow(Eext_red,[]),title('Eextern after volcano (red)')
        figure,imshow(I,[]); hold on, quiver(y,x,Fext_red(1:10:end,1:10:end,2),Fext_red(1:10:end,1:10:end,1));        
        title('The external force field after volcano (red)')
        [x,y]=ndgrid(1:10:size(Fext_blue,1),1:10:size(Fext_blue,2));
        figure,imshow(Eext_blue,[]),title('Eextern after volcano (blue)')
        figure,imshow(I,[]); hold on, quiver(y,x,Fext_blue(1:10:end,1:10:end,2),Fext_blue(1:10:end,1:10:end,1));        
        title('The external force field after volcano (blue)')
    end
    % 4.2 �����ߵ��ڲ�����������4ʽ��Gint�ĵ�����
    S=SnakeInternalForceMatrix2D(xlen,Options.Alpha,Options.Beta);  % ���ܵ�������
    for i=1:Options.Iterations
    % 4.3 �����������ƶ��ߵ��˶���
        P=SnakeMoveIteration2D(S,P,Fext_red,Gamma_red,Options.Kappa,Options.Delta,eliminate_x_force);
        Q=SnakeMoveIteration2D(S,Q,Fext_blue,Gamma_blue,Options.Kappa,Options.Delta,eliminate_x_force);
    end
end