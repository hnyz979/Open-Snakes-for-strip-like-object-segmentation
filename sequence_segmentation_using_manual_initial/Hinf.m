function [xnew,K] = Hinf(xold,y,u,F,G,H,L,R,Q,S,I,theta)
Sdot = L' * S * L;
E = (H' * R^(-1) * H  - theta * Sdot)^(-1); 
[P,~,~] = dare(F',I,Q,E);
v = P^(-1) - theta*Sdot + H'*R^(-1)*H;
eigValue = eig(v);
if ~isempty(find(eigValue<=0))  % 判断矩阵v是否正定
    error('func_Hinf_filter():矩阵P^(-1) - theta*Sdot + H^T*R^(-1)*H非正定，建议重新选择theta');
end
K = P * (I - theta * Sdot * P + H' * R^(-1) * H * P)^(-1) * H' * R^(-1);
% 进入H-infinity计算
xnew = F * xold + G * u + F * K * (y-H*xold);
end