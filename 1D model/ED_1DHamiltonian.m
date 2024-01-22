clear;
%%%----------------------——————————————————————————————-------------%%%%%%
% The following demonstration pertains to the one-dimensional 
% XXZ model (Heisenberg Model), Ising model, 
% and their representations in terms of transverse and longitudinal fields. 
% Users are encouraged to selectively comment out other sections as needed when employing a specific model. 
% It is important to note that the modules for the Hamiltonian matrix 
% cannot be utilized concurrently, but transverse and longitudinal fields 
% can be simultaneously implemented.
%%%----------------------——————————————————————————————-------------%%%%%%


% Preceding specifications for the following variables:
% hb, sh, szh represent variables associated with the longitudinal field;
% Gamma, sgamma, sxg denote variables associated with the transverse field;
% It is advisable to refrain from modifying variables related to the Hamiltonian matrix.

% Note that J has not been introduced explicitly in the Heisenberg model, 
% and by default, J is set to 1. If its utilization is desired, please exercise caution.
%%%----------------------——————————————————————————————-------------%%%%%%


% 格点数 numbers of spin
L = 3;

J = 1; % coupling

% XXZ 参数: XXZ perematers
Delta = 1.0;

% 纵场强度 longitudinal fields
hb = 0;

% 横场强度 transverse fields
gamma = 0;

% 准备算符矩阵 operators:
sz = 0.5*[1.0,0;0,-1.0];
sp = [0,1;0,0];
id = eye(2,2); % 2x2 I matrix

% 初始化 initial:
H = zeros(2,2);% 初始化为0矩阵
Sp0 = sp; %S+
Sz0 = sz;
Id = eye(2,2); %单位矩阵 I matrix

Sp1 = Sp0;
Sz1 = Sz0;
% Sp2 = sp';%求转置的操作，果然变成了S-



%%%%-------------------------------------------------------------------------%
%% 1D Heisenberg Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------------------------------------------------------------------------%
% for site=1:L-1 %L-1
%     h = 0.5*kron(Sp1, sp'); %kronecker 张量积,此时的Sp1是site-1的
%     H = kron(H, id) + h + h' + kron(Sz1,sz);%完成一个哈密顿量 此时的Sz1是site-1的
% 
%     if L>2
%         %对最后一个格点的处理
%         if site == L-1
%             h = 0.5*kron(Sp0, sp');
%             H = H + h + h' + kron(Sz0,sz);
%         end
% 
%         %
%         if site < L-1
%             Sp1 = kron(Id, sp); %此时的Sp1是site的
%             Sz1 = kron(Id, sz); %此时的Sz1是site的
% 
%             Sp0 = kron(Sp0, id);%此时的Sp0是site的
%             Sz0 = kron(Sz0, id);%此时的Sz0是site的
%             Id = eye(size(Sp0)); 
%         end
%     end
% end
%%%%-------------------------------------------------------------------------%
%% 1D Heisenberg Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------------------------------------------------------------------------%



% %-------------------------------------------------------------------------%
% %%% 1D XXZ Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %-------------------------------------------------------------------------%
% for site=1:L-1 %L-1
%     h = 0.5*kron(Sp1, sp'); %kronecker 张量积,此时的Sp1是site-1的
% 
%     H = kron(H, id)+ h + h' + Delta * kron(Sz1,sz);
%     %完成一个哈密顿量 此时的Sz1是site-1的
% 
% 
%     if L>2
%         %对最后一个格点的处理
%         if site == L-1
%             h = 0.5*kron(Sp0, sp');
%             H = H + h + h' + Delta * kron(Sz0,sz);
%         end
% 
% 
%         if site < L-1
%             Sp1 = kron(Id, sp); %此时的Sp1是site的
%             Sz1 = kron(Id, sz); %此时的Sz1是site的
% 
%             Sp0 = kron(Sp0, id);%此时的Sp0是site的
%             Sz0 = kron(Sz0, id);%此时的Sz0是site的
%             Id = eye(size(Sp0)); 
%         end
%     end
% end
% %-------------------------------------------------------------------------%
% %%% 1D XXZ Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %-------------------------------------------------------------------------%


% %-------------------------------------------------------------------------%
% %%% 纵场 longitudinal fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %-------------------------------------------------------------------------%
sh = zeros(2,2);
idnew = id;
szh = sz;
for site=1:L-1 %L-1

    %site-1的：
    sh = kron(sh,id) + hb * kron(szh,id);%初始化

    if L>2
        %对最后一个格点的处理
        if site == L-1
            sh = sh + hb*kron(idnew,sz);
        end
    
    
        if site < L-1
            szh = kron(idnew, sz); %此时的Sz1是site的
            idnew = eye(size(szh)); 
        end
    else %2-spin的情况
    % 纵场对最后一个格点(site)的处理(2-spin)
    sh = sh + hb * kron(id,sz);
    end

end
% %-------------------------------------------------------------------------%
% %%% 纵场 longitudinal fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%%%% 横场 transverse fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
sx = [0,1.0;1.0,0];%注意，这里用的是sigam矩阵，不是S矩阵 sigma^x
iidnew = id;
sxg = sx;
sgamma = zeros(2,2);
for site=1:L-1 %L-1

    %site-1的：
    sgamma = kron(sgamma,id) + gamma * kron(sxg,id);%初始化

    if L>2
        %对最后一个格点的处理
        if site == L-1
            sgamma = sgamma + gamma * kron(iidnew,sx);
        end


        if site < L-1
            sxg = kron(iidnew, sx); %此时的Sz1是site的
            iidnew = eye(size(sxg)); 
        end
    else %2-spin的情况
    % 纵场对最后一个格点(site)的处理(2-spin)
    sgamma = sgamma + gamma * kron(id,sx);
    end

end

%%%若处理的是S=1/2模型，要用S矩阵，因此对矩阵整体有1/2系数：
% sgamma = 1/2* sgamma;

%-------------------------------------------------------------------------%
%%% 横场 transverse fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%





%%%%-------------------------------------------------------------------------%
%% 1D Ising Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------------------------------------------------------------------------%

for site=1:L-1 %L-1
    H = kron(H, id) + J*kron(Sz1,sz);
    %完成一个哈密顿量 此时的Sz1是site-1的


    if L>2
        %对最后一个格点的处理
        if site == L-1
            H = H + J*kron(Sz0,sz);
        end
    
    
        if site < L-1
            Sz1 = kron(Id, sz); %此时的Sz1是site的
            Sz0 = kron(Sz0, id);%此时的Sz0是site的
            Id = eye(size(Sz0)); 
        end
    end
end

H = 4*H; %消去S矩阵表示1/4系数的影响

%%%%-------------------------------------------------------------------------%
%% 1D Ising Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------------------------------------------------------------------------%


%加上纵场：H + longitudinal fields
H = H + sh;

%加上横场: H + transverse fields
H = H + sgamma; 


[V,D] = eig(H);
% det(H)
D = diag(D)/L;
GS = V(:,1); % gound state
