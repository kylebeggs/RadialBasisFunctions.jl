# Radial Basis Functions Theory

Radial Basis Functions (RBF) use only a distance (typically Euclidean) when constructing the basis. For example, if we wish to build an interpolator we get the following linear combination of RBFs

```math
f(\mathbf{x})=\sum_{i=1}^{N} \alpha_{i} \phi(\lvert \mathbf{x}-\mathbf{x}_{i} \rvert)
```

where ``\mid \cdot \mid`` is a norm (we will use Euclidean from here on) and so ``\lvert \mathbf{x}-\mathbf{x}_{i} \rvert = r`` is the Euclidean distance (although it can be any) and ``N`` is the number of data points.

There are several types of RBFs to choose from, some with a tunable shape parameter, ``\varepsilon``. Here are some popular ones:

| Type                 | Function                               |
| -------------------- | -------------------------------------- |
| Polyharmonic Spline  | ``\phi(r) = r^n`` where ``n=1,3,5,7,\dots``                      |
| Multiquadric             |    ``\phi(r)=\sqrt{ (r \varepsilon)^{2}+ 1 }``                                    |
| Inverse Multiquadric | ``\phi(r) = 1 / \sqrt{(r \varepsilon)^2+1}}`` |
| Gaussian             | ``\phi(r) = e^{-(r \varepsilon)^2}``               |

## Augmenting with Monomials

The interpolant may be augmented with a polynomial as

```math
f(\mathbf{x})=\sum_{i=1}^{N} \alpha_{i} \phi(\lvert \mathbf{x}-\mathbf{x}_{i} \rvert) + \sum_{i=1}^{N_{p}} \gamma_{i} p_{i}(\mathbf{x})
```

where ``N_{p}=\begin{pmatrix} m+d \\ m \end{pmatrix}`` is the number of monomials (``m`` is the monomial order and ``d`` is the dimension of ``\mathbf{x}``) and ``p_{i}(\mathbf{x})`` is the monomial term, or:

```math
p_{i}(\mathbf{x})=q_{i}(\lvert \mathbf{x}-\mathbf{x}_{i} \rvert)
```

where ``q_{i}`` is the ``i``-th monomial in ``\mathbf{q}=\begin{bmatrix} 1, x, y, x^2, xy, y^2 \end{bmatrix}`` in 2D, for example. By collocation the expansion of the augmented interpolant at all the nodes ``\mathbf{x}_{i}`` where ``i=1\dots N``, there results a linear system for the interpolant weights as:

```math
\begin{bmatrix}
\mathbf{A} & \mathbf{P} \\
\mathbf{P}^\mathrm{T} & 0
\end{bmatrix}
\begin{bmatrix}
\boldsymbol{\alpha} \\
\boldsymbol{\gamma}
\end{bmatrix}=
\begin{bmatrix}
\mathbf{f} \\
0
\end{bmatrix}
```

where

```math
\mathbf{A}=
\begin{bmatrix}
\phi(\lvert \mathbf{x}_{1}-\mathbf{x}_{1} \rvert) & \dots & \phi(\lvert \mathbf{x}_{1}-\mathbf{x}_{N} \rvert) \\
\vdots & & \vdots \\
\phi(\lvert \mathbf{x}_{N}-\mathbf{x}_{1} \rvert) & \dots & \phi(\lvert \mathbf{x}_{N}-\mathbf{x}_{N} \rvert)
\end{bmatrix}
\hspace{2em}
\mathbf{p}=
\begin{bmatrix}
p_{1}(\mathbf{x}_{1}) & \dots & p_{N}(\mathbf{x}_{1}) \\
\vdots & & \vdots \\
p_{1}(\mathbf{x}_{N}) & \dots & p_{N}(\mathbf{x}_{N})
\end{bmatrix}
```

and ``\mathbf{f}`` is the vector of dependent data points

```math
\mathbf{f}=
\begin{bmatrix}
f(\mathbf{x}_{1}) \\
\vdots \\
f(\mathbf{x}_{N})
\end{bmatrix}
```

and ``\boldsymbol{\alpha}`` and ``\boldsymbol{\gamma}`` are the interpolation coefficients. Note that the equations relating to ``\mathbf{P}^\mathrm{T}`` are included to ensure optimal interpolation and unique solvability given that conditionally positive radial functions are used and the nodes in the subdomain form a unisolvent set. See (Fasshauer, et al. - Meshfree Approximation Methods with Matlab) and (Wendland, et al. - Scattered Data Approximation).

This augmentation of the system is highly encouraged for a couple main reasons:

1. Increased accuracy especially for flat fields and near boundaries
2. To ensure the linear system has a unique solution

See (Flyer, et al. - On the role of polynomials in RBF-FD approximations: I. Interpolation and accuracy) for more information on this.

## Local Collocation

The original RBF method employing the Kansa approach which connects all the nodes in the domain and, as such, is a _global_ method. Due to ill-conditioning and computational cost, this approach scales poorly; therefore, a _local_ approach is used instead. In the _local_ approach, each node is influenced only by its ``k`` nearest neighbors which helps solve the issues related to _global_ collocation.

## Constructing an Operator

In the Radial Basis Function - Finite Difference method (RBF-FD), a stencil is built to approximate derivatives using the same neighborhoods/subdomains of $N$ points. This is used in the [[MeshlessMultiphysics.jl]] package. For example, if ``\mathcal{L}`` represents a linear differential operator, one can express the differentiation of the field variable ``u`` at the center of the subdomain ``\mathbf{x}_{c}`` in terms of some weights ``\mathbf{w}`` and the field variable values on all the nodes within the subdomain as

```math
\mathcal{L}u(\mathbf{x}_{c}) = \sum_{i=1}^{N}w_{i}u(\mathbf{x}_{i})
```

We can find $\mathbf{w}$ by satisfying

```math
\sum_{i=1}^{N}w_{i}\phi_{j}(\mathbf{x}_{i}) = \mathcal{L}\phi_{j}(\mathbf{x}_{c})
```

for each $\phi_{j}$ where $j=1,\dots, N$ and if you wish to augment with monomials, we also must satisfy

```math
\sum_{i=1}^{N_{p}}\lambda_{i}p_{j}(\mathbf{x}_{i}) = \mathcal{L}p_{j}(\mathbf{x}_{c})
```

which leads to an overdetermined problem

```math
\mathrm{min} \left( \frac{1}{2} \mathbf{w}\mathbf{A}^{\intercal}\mathbf{w} - \mathbf{w}^{\intercal} \mathcal{L}\phi \right), \text{ subject to } \mathbf{P}^{\intercal}\mathbf{w}=\mathcal{L}\mathbf{p}
```

which is practically solved as a linear system for the weights $\mathbf{w}$ as

```math
\begin{bmatrix}\mathbf{A} & \mathbf{P} \\
\mathbf{P}^\mathrm{T} & 0
\end{bmatrix}
\begin{bmatrix}
\mathbf{w} \\
\boldsymbol{\lambda}
\end{bmatrix}=
\begin{bmatrix}
\mathcal{L}\boldsymbol{\phi} \\
\mathcal{L}\mathbf{p}
\end{bmatrix}
```

where ``\boldsymbol{\lambda}`` are treated as Lagrange multipliers and are discarded after solving the linear system and

```math
\mathcal{L}\boldsymbol{\phi}=
\begin{bmatrix}
\mathcal{L}\boldsymbol{\phi}(\lvert \mathbf{x}_{1}-\mathbf{x}_{c} \rvert) \\
\vdots \\
\mathcal{L}\boldsymbol{\phi}(\lvert \mathbf{x}_{N}-\mathbf{x}_{c} \rvert)
\end{bmatrix}
\hspace{2em}
\mathcal{L}\mathbf{p}=
\begin{bmatrix}
\mathcal{L}p_{1}(\mathbf{x}_{c}) \\
\vdots \\
\mathcal{L}p_{N_{p}}(\mathbf{x}_{c})
\end{bmatrix}
```
