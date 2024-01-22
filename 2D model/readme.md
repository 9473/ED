# 2D Ising Model

We have already presented the one-dimensional representation of Hamiltonians with off-diagonal terms. Readers may notice that the construction of off-diagonal terms is consistent with the $S^zS^z$ interaction. Next, we will primarily focus on the construction of two-dimensional Hamiltonians, with a main emphasis on the operator construction involving $S^zS^z$. The construction of off-diagonal terms follows a similar approach.

Its overall construction approach is as follows:

```matlab
clear;
L = 2;
J = 1;

sz = 0.5*[1.0,0;0,-1.0];
sp = [0,1;0,0];
id = eye(2,2);

H = zeros(2,2);
Sp0 = sp; 
Sz0 = sz;
Id = eye(2,2); 

Sp1 = Sp0;
Sz1 = Sz0;

row=0;
for site=1:L-1 %L-1
    h = 0.5*kron(Sp1, sp'); 
    H = kron(H, id) + h + h' + kron(Sz1,sz);

    if L>2
        if site == L-1
            h = 0.5*kron(Sp0, sp');
            H = H + h + h' + kron(Sz0,sz);
        end


        if site < L-1
            Sp1 = kron(Id, sp); 
            Sz1 = kron(Id, sz); 

            Sp0 = kron(Sp0, id);
            Sz0 = kron(Sz0, id);
            Id = eye(size(Sz0)); 
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


 for row=1:L-1
        II=eye(2^(row*L),2^(row*L));
        III=eye(2^(row*L-L),2^(row*L-L));  
        Sp0 = kron(II,sp);
        Sz0 = kron(II,sz); 
        Sp1 = kron(II,sp);
        Sz1 = kron(II,sz); 
        Id = eye(size(Sp0));

        for i=1:L
            Sp{i}=kron(III,ssp{i});
            Sz{i}=kron(III,ssz{i});
        end
    for site=0:L-1
        h = 0.5*kron(kron(Sp{site+1},eye(2^site,2^site)), sp');
        H = kron(H, id) + h + h' + kron(kron(Sz{site+1},eye(2^site,2^site)),sz);

        if site>0
            % h = 0.5*kron(Sp1, sp');
            H = H + kron(Sz1,sz);


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
        h = 0.5*kron(ssp{i},Sp{i});
        H = H + h + h' + kron(ssz{i},Sz{i});
    end
```

For a 2x2-spin system, the Hilbert space is given by $D = 2^{4} = 16$, and consequently, the size of the final Hamiltonian matrix should be 16x16. We mainly introduce 3x3-spin systems.
##### 3x3-spin system

Hamiltonian
$$H_1 = \color{blue}S_1 S_2\color{black} I_3I_4I_5I_6I_7I_8I_9$$
$$H_2 = \color{blue}S_1 \color{black} I_2  I_3  \color{blue}S_4 \color{black}I_5I_6I_7I_8I_9$$
$$H_3 = I_1 \color{blue}S_2 S_3  \color{black}I_4 I_5I_6I_7I_8I_9$$
$$H_4 = I_1\color{blue}S_2\color{black}I_3I_4\color{blue}S_5\color{black}I_6I_7I_8I_9$$
$$H_5 = \color{blue}S_1\color{black}I_2\color{blue}S_3\color{black}I_4I_5I_6I_7I_8I_9$$
$$H_6 = I_1I_2\color{blue}S_3\color{black}I_4I_5\color{blue}S_6\color{black}I_7I_8I_9$$
$$H_7 = I_1I_2I_3\color{blue}S_4S_5\color{black}I_6I_7I_8I_9$$
$$H_8 = I_1I_2I_3\color{blue}S_4\color{black}I_5I_6\color{blue}S_7\color{black}I_8I_9$$
$$H_9 = I_1I_2I_3I_4\color{blue}S_5S_6\color{black}I_7I_8I_9$$
$$H_{10} = I_1I_2I_3I_4\color{blue}S_5\color{black}I_6I_7\color{blue}S_8\color{black}I_9$$
$$H_{11} = I_1I_2I_3\color{blue}S_4\color{black}I_5\color{blue}S_6\color{black}I_7I_8I_9$$
$$H_{12} = I_1I_2I_3I_4I_5\color{blue}S_6\color{black}I_7I_8\color{blue}S_9$$
$$H_{13} = I_1I_2I_3I_4I_5I_6\color{blue}S_7S_8\color{black}I_9$$
$$H_{14} = \color{blue}S_1\color{black} I_2I_3I_4I_5I_6\color{blue}S_7\color{black}I_8I_9$$
$$H_{15} = I_1I_2I_3I_4I_5I_6I_7\color{blue}S_8S_9$$
$$H_{16} = I_1\color{blue}S_2\color{black}I_3I_4I_5I_6I_7\color{blue}S_8\color{black}I_9$$
$$H_{17} = I_1I_2I_3I_4I_5I_6\color{blue}S_7\color{black}I_8\color{blue}S_9$$
$$H_{18} = I_1I_2\color{blue}S_3\color{black}I_4I_5I_6I_7I_8\color{blue}S_9$$

#### Construction of the Transverse Hamiltonian

```matlab
row=0;
for site=1:L-1 
    H = kron(H, id) + kron(Sz1,sz);

    if L>2
        if site == L-1
            H = H + kron(Sz0,sz);
        end

        if site < L-1

            Sz1 = kron(Id, sz); 

            Sz0 = kron(Sz0, id);
            Id = eye(size(Sz0)); 
        end
    end
end
```

At `row=0`, by adopting the one-dimensional Hamiltonian construction method, we can successfully construct:
$$H =[\frac{1}{2}(S^+_1 \otimes S^-_2\otimes I_3 + S^-_1 \otimes S^+_2\otimes I_3) + S^z_1\otimes  S^z_2\otimes I_3$$
$$+\frac{1}{2}( I_1\otimes S^+_2 \otimes S^-_3 + I_1\otimes S^-_2 \otimes S^+_3) +I_1 \otimes  S^z_2\otimes  S^z_3$$
$$+\frac{1}{2}( S^+_1 \otimes I_2 \otimes S^-_3 +  S^-_1 \otimes I_2 \otimes S^+_3) +  S^z_1 \otimes I_2\otimes  S^z_3]$$
Automatic completion will be performed during the subsequent `kron(H, id)` operation.

```matlab
for row=1:L-1
				II=eye(2^(row*L),2^(row*L)); 
      
        Sz0 = kron(II,sz); 
        Sz1 = kron(II,sz); 
        Id = eye(size(Sz0));
	for site=0:L-1
        H = kron(H, id) + Longitudinal Hamiltonian;

        if site>0
            H = H + kron(Sz1,sz);

            if site == L-1
                H = H + kron(Sz0,sz);
            end

            if site < L-1
            		Sz1 = kron(Id, sz);
                Sz0 = kron(Sz0, id);
                Id = eye(size(Sz0));        
            end

        end
  end
end
```

At `row=1`, $II = I_{2^3\times 2^3}$.(L=3)

```matlab
% Sz1 = kron(II,sz);
H = ... + kron(Sz1,sz);
```

So:
$$I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \color{orange} \otimes S^z_{2\times 2}$$
1. For `0<site<L-1=3-1=2`：

```matlab
Sz1 = kron(Id, sz);
Sz0 = kron(Sz0, id);
Id = eye(size(Sz0)); 
```

$$\text{Sz1: }I_{16\times 16}\color{orange} \otimes S^z_{2\times 2},\ \ \\ \text{Sz0: }I_{8\times 8} \otimes S^z_{2\times 2} \color{orange} \otimes I_{2\times 2}$$

2. For `site=2=L-1`

```matlab
H = ... + kron(Sz1,sz);
```

$$I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2} \otimes  S^z_{2\times 2}\color{orange} \otimes S^z_{2\times 2} \to 64\times 64 \text{ matrix}$$

And,

```matlab
if site == L-1
    h = 0.5*kron(Sp0, sp');
    H = H + h + h' + kron(Sz0,sz);
end
```

$$I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}  \otimes I_{2\times 2}\color{orange} \otimes S^z_{2\times 2} \to 64\times 64 \text{ matrix}$$

At this point, the realization includes $H_7, H_9, H_{11}$.

However, for a 3x3-spin system, there is an additional loop when `row=L-1=2`. It is easy to anticipate that this loop can implement $H_{13},H_{15},H_{17}$.

#### Construction of the Longitudinal Hamiltonian:

Operator Preparation Part I:

```matlab
for i=1:L
    Sz{i}=1;
end
for i=1:L
for j=1:L
    if (j==i)
        Sz{j}=kron(Sz{j},sz);
    else
        Sz{j}=kron(Sz{j},id);
    end
end
end

ssz=Sz;
```

For a **3x3-spin** system, the output of Sz{i} is a 1x3 array, with each element being an 8x8 matrix:


```bash
Sz =

  1×3 cell 数组

  {8×8 double}    {8×8 double}    {8×8 double}
```

```bash
Sz{1} =   
    0.5000         0         0         0         0         0         0         0
         0    0.5000         0         0         0         0         0         0
         0         0    0.5000         0         0         0         0         0
         0         0         0    0.5000         0         0         0         0
         0         0         0         0   -0.5000         0         0         0
         0         0         0         0         0   -0.5000         0         0
         0         0         0         0         0         0   -0.5000         0
         0         0         0         0         0         0         0   -0.5000
```

It's:
$$\text{Sz(1)} :S^z_{2\times 2} \otimes I_{4\times 4}=S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}$$


```bash
Sz{2} = 
		0.5000         0         0         0         0         0         0         0
         0    0.5000         0         0         0         0         0         0
         0         0   -0.5000         0         0         0         0         0
         0         0         0   -0.5000         0         0         0         0
         0         0         0         0    0.5000         0         0         0
         0         0         0         0         0    0.5000         0         0
         0         0         0         0         0         0   -0.5000         0
         0         0         0         0         0         0         0   -0.5000
```

$$\text{Sz(2)}: I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}$$

```bash
Sz{3} = 
    0.5000         0         0         0         0         0         0         0
         0   -0.5000         0         0         0         0         0         0
         0         0    0.5000         0         0         0         0         0
         0         0         0   -0.5000         0         0         0         0
         0         0         0         0    0.5000         0         0         0
         0         0         0         0         0   -0.5000         0         0
         0         0         0         0         0         0    0.5000         0
         0         0         0         0         0         0         0   -0.5000
```

$$\text{Sz(3)}: I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}$$




During the construction of the Longitudinal Hamiltonian:

```matlab
for row=1:L-1
        III=eye(2^(row*L-L),2^(row*L-L)); 
        for i=1:L
            Sz{i}=kron(III,ssz{i});
        end
    for site=0:L-1
        H = kron(H, id) + kron(kron(Sz{site+1},eye(2^site,2^site)),sz);

        if site>0
            Transverse Hamiltonian
        end
    end
 end
 
 row=L;
 ...
```

There,

```matlab
for i=1:L
    Sz{i}=kron(III,ssz{i});
end
% III=eye(2^(row*L-L),2^(row*L-L))
```
##### At `row=1`,

$$III = I_{2^0,2^0} = I_{1\times 1}.$$

1. For `site=0`:

```matlab
kron(kron(Sz{site+1},eye(2^(site),2^(site)),sz)
```

$$\text{Sz(1)} \otimes I_{1\times 1}$$ 
$$\to S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \color{orange}\otimes I_{1\times1}\\ \to S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \color{orange}\otimes I_{1\times1}\color{red}\otimes S^z_{2\times 2} \to 16\times 16 \text{ matrix}$$

2. For `site=1`

```
H = kron(H,id) +...+ kron(kron(Sz{site+1},eye(2^site,2^site)),sz)
```

$$S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes S^z_{2\times 2}\color{orange}\otimes I_{2\times2}$$
$$+I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\color{orange}\otimes I_{2\times 2}\otimes S^z_{2\times 2}$$

gives 32x32 matrix.

3. For `site=2=L-1`

```
H = kron(H,id) +...+ kron(kron(Sz{site+1},eye(2^site,2^site)),sz)
```

$$S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes S^z_{2\times 2}\otimes I_{2\times2}\color{orange}\otimes I_{2\times 2}$$
$$+I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\color{orange}\otimes I_{2\times 2}$$
$$+I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}\color{orange}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}$$


##### At `row=2=L-1`

$$III = I_{2^3\times 2^3} = I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}$$

1. For `site=0`:

```matlab
kron(kron(Sz{site+1},eye(2^(site),2^(site)),sz)
```

$$I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes \color{orange} I_{1\times1}$$
$$\to I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes \color{red} S^z_{2\times 2} \to 2^7\times 2^7 \text{ matrix}$$



2. For `site=1`

```
H = kron(H,id) +...+ kron(kron(Sz{site+1},eye(2^site,2^site)),sz)
```

$$I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes \color{red} S^z_{2\times 2}\color{orange}\otimes I_{2\times2}$$
$$+I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\otimes \color{orange}I_{2\times 2}\otimes S^z_{2\times 2}$$

gives $2^8 \times 2^8$ matrix.

3. For `site=2=L-1`

```
H = kron(H,id) +...+ kron(kron(Sz{site+1},eye(2^site,2^site)),sz)
```

$$I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes \color{red} S^z_{2\times 2}\color{orange}\otimes I_{2\times2}\otimes I_{2\times 2}$$
$$+I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\otimes \color{orange}I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}$$
$$+I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}\otimes \color{orange}I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}$$

Gives $2^9 \times 2^9$ matrix.

Therefore, we summarize the overall Longitudinal Hamiltonian when `row = 1, 2`:


$$=S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes S^z_{2\times 2}\otimes I_{2\times2}\color{orange}\otimes I_{2\times 2}\color{green}\otimes I_{8\times 8}$$
$$+I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\color{orange}\otimes I_{2\times 2}\color{green}\otimes I_{8\times 8}$$
$$+I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}\color{orange}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\color{green}\otimes I_{8\times 8}$$
$$+I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes \color{red} S^z_{2\times 2}\color{orange}\otimes I_{2\times2}\color{orange}\otimes I_{2\times 2}$$
$$+I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\otimes \color{orange}I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}$$
$$+I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}\otimes \color{orange}I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}$$


The green part is the completion caused by the `kron(H,id)` in the loop. Thus, we observe that this section of the code accomplishes the implementation of $H_2, H_4, H_6, H_{8}, H_{10}, H_{12}$.

So, there is another part to implement in `row=L`:

##### At `row=L`

```fortran
row = L;
III=eye(2^(row*L-2*L),2^(row*L-2*L));
for i=1:L
    Sp{i}=kron(III,ssp{i});
    Sz{i}=kron(III,ssz{i});
end
for i=1:L
		H = H + kron(ssz{i},Sz{i})
end
```

$$III : (row-2)*L = 3 \to I_{2^3 \times 2^3}$$

$$\text{kron(ssz(i)},\color{red}\text{Sz(i})\color{black})$$  
$$S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}  \color{red} \otimes I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}$$
$$+ I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}\color{red }\otimes I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2}\otimes I_{2\times 2}\otimes S^z_{2\times 2}\otimes I_{2\times 2}$$
$$+ I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}\color{red } \otimes I_{2\times 2} \otimes I_{2\times 2}\otimes I_{2\times 2} \otimes I_{2\times 2} \otimes I_{2\times 2} \otimes S^z_{2\times 2}$$

Therefore, `row=L` completes the construction of the Longitudinal Hamiltonian.

It is evident that the construction of the Longitudinal Hamiltonian involves building the first row at `row=1`, the (L-1)th row at `row=L-1`, and the last row separately at `row=L`. The Hamiltonian of the last row is periodically "connected" to the elements of the first row.
