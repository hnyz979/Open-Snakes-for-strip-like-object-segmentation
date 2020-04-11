function J=DrawSegmentedArea2D(P,Isize)
% 画出闭合轮廓，轮廓内部的为1，外部的为0。
% 输入：
%  P        2*N的轮廓点，不必须首尾相接
%  Isize    输出的图像大小
% 输出：
%  J        二值图像，轮廓内部的为1，外部的为0。
% 例如：
%   y=[182 233 251 205 169];
%   x=[163 166 207 248 210];
%   P=[x(:) y(:)];
%   J=DrawSegmentedArea2D(P,[400 400]);
%   figure, imshow(J); 
%   hold on; plot([P(:,2);P(1,2)],[P(:,1);P(1,1)]);
% 忘记这个函数是从哪里借用的了，但是效果不错。。。

J=false(Isize+2);
% Loop through all line coordinates
x=round([P(:,1);P(1,1)]); x=min(max(x,1),Isize(1));
y=round([P(:,2);P(1,2)]); y=min(max(y,1),Isize(2));
for i=1:(length(x)-1)
   % Calculate the pixels needed to construct a line of 1 pixel thickness
   % between two coordinates.
   xp=[x(i) x(i+1)];  yp=[y(i) y(i+1)]; 
   dx=abs(xp(2)-xp(1)); dy=abs(yp(2)-yp(1));
   if(dx==dy)
     if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
     if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
   elseif(dx>dy)
     if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
     yline=linspace(yp(1),yp(2),length(xline));
   else
     if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
     xline=linspace(xp(1),xp(2),length(yline));   
   end
   % Insert all pixels in the fill image
   J(round(xline+1)+(round(yline+1)-1)*size(J,1))=1;
end
J=bwfill(J,1,1); J=~J(2:end-1,2:end-1);