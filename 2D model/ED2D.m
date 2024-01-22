%2D Heisenberg model

clear;

% 格点数
L = 3;

J = 1;

% 准备算符矩阵
sz = 0.5*[1.0,0;0,-1.0];
sp = [0,1;0,0];
id = eye(2,2);

% 初始化
H = zeros(2,2);% 初始化为0矩阵
Sp0 = sp; %S+
Sz0 = sz;
Id = eye(2,2); %单位矩阵

Sp1 = Sp0;
Sz1 = Sz0;
% Sp2 = sp';%求转置的操作，果然变成了S-



row=0;
for site=1:L-1 %L-1
    h = 0.5*kron(Sp1, sp'); %kronecker 张量积,此时的Sp1是site-1的
    H = kron(H, id) + h + h' + kron(Sz1,sz);%完成一个哈密顿量 此时的Sz1是site-1的

    if L>2
        %对最后一个格点的处理
        if site == L-1
            h = 0.5*kron(Sp0, sp');
            H = H + h + h' + kron(Sz0,sz);
        end

        %
        if site < L-1
            Sp1 = kron(Id, sp); %此时的Sp1是site的
            Sz1 = kron(Id, sz); %此时的Sz1是site的

            Sp0 = kron(Sp0, id);%此时的Sp0是site的
            Sz0 = kron(Sz0, id);%此时的Sz0是site的
            Id = eye(size(Sp0)); 
        end
    end
end

for i=1:L
    Sp{i}=1;
    Sz{i}=1;
end
for i=1:L
for j=1:L
    if (j==i)
        Sp{j}=kron(Sp{j},sp);
        Sz{j}=kron(Sz{j},sz);
    else
        Sp{j}=kron(Sp{j},id);
        Sz{j}=kron(Sz{j},id);
    end
end
end

ssp=Sp;
ssz=Sz;
    %准备好了row1行的算符

for row=1:L-1
    II=eye(2^(row*L),2^(row*L));
    III=eye(2^(row*L-L),2^(row*L-L));
    Sp0 = kron(II,sp);%row行的算符
    Sz0 = kron(II,sz);
    Sp1 = kron(II,sp);
    Sz1 = kron(II,sz);
    Id = eye(size(Sp0));

    for i=1:L
        Sp{i}=kron(III,ssp{i});%row-1行的算符
        Sz{i}=kron(III,ssz{i});
    end
    
    for site=0:L-1
        h = 0.5*kron(kron(Sp{site+1},eye(2^site,2^site)), sp'); %kronecker 张量积
        H = kron(H, id) + h + h' + kron(kron(Sz{site+1},eye(2^site,2^site)),sz);%完成一个纵向哈密顿量
    
        if site>0
            h = 0.5*kron(Sp1, sp'); %kronecker 张量积
            H = H + h + h' + kron(Sz1,sz);%完成一个横向哈密顿量
    
    
            if site == L-1
                h = 0.5*kron(Sp0, sp');
                H = H + h + h' + kron(Sz0,sz);
            end
    
            if site < L-1
                Sp1 = kron(Id, sp);
                Sz1 = kron(Id, sz);
    
                Sp0 = kron(Sp0, id);
                Sz0 = kron(Sz0, id);
                Id = eye(size(Sp0));        
            end
    
        end
    end
end

row=L;
III=eye(2^(row*L-2*L),2^(row*L-2*L));
for i=1:L
    Sp{i}=kron(III,ssp{i});
    Sz{i}=kron(III,ssz{i});
end

for i=1:L
    h = 0.5*kron(ssp{i},Sp{i}); %kronecker 张量积
    H = H + h + h' + kron(ssz{i},Sz{i});%完成一个纵向哈密顿量
end




[V,D] = eig(H);
% det(H)
D = diag(D)/L;
GS=V(:,1) % gound state