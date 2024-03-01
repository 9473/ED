# ED(Exact diagonalization)

------
Exact diagonalization method applied to quantum spin systems, dealing with some simple models: Ising Model + transverse field, $S=1/2$ Heisenberg Model, XXZ Model + longitudinal field  

------

### Ising Model
the Hamiltonian of Ising Model with transverse field gives:  

$$
H_{TFIM} = -J\sum_{\langle i,j \rangle}\sigma^z_i\sigma^z_j - \sum_{i} h_{x}\sigma^x_i
$$

If $h_x$ is set to $0$ , the model returns to the simplest standard Ising model. Here, $h_x$ is a constant and we use a uniform transverse field.

------
### XXZ Model
Now we give the Hamiltonian of XXZ model with longitudinal field,  

$$
H =\sum_{\langle i,j\rangle} \left[ (S^x_i  S^x_j + S^y_i S^y_j) + \Delta S^z_i S^z_j \right]- \sum_i h^z S^z_i
$$

1. where $\Delta =\frac{J_z}{J_{x,y}}$ . If $h^z=0$ gives the XXZ model: $H =\sum_{\langle i,j\rangle} \left[ (S^x_i  S^x_j + S^y_i S^y_j) + \Delta S^z_i S^z_j \right]$ .

2. $\Delta = 1.0$ in XXZ model gives the standard Heisenberg Model: $H = J\sum_{<i,j>}\mathbf{S}_i\cdot \mathbf{S}_j$ $\to$ .  

$$
H = J\sum_{\langle i,j \rangle } (S^x_iS^x_j+S^y_iS^y_j+S^z_iS^z_j)
$$

3. Our more common form is to use $S^+, S^-$,

$$
H=\sum_{\langle i,j\rangle} \left[ \frac{1}{2}(S^+_i S^-_j + S^-_i S^+_j) + S^z_i S^z_j \right] - \sum_i h^z S^z_i
$$

