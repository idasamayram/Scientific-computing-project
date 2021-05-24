### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 524d5cef-debd-47a3-8ff5-2f25341d57b3
begin 
	ENV["LANG"]="C"
	using Printf
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add(["Plots", "PyPlot","PlutoUI","ExtendableGrids", "GridVisualize","Documenter","VoronoiFVM","Triangulate","SimplexGridFactory","DataFrames"])
	using PlutoUI,PyPlot,ExtendableGrids,VoronoiFVM,Triangulate,DataFrames,SimplexGridFactory,GridVisualize,Documenter
TableOfContents()
end

# ╔═╡ 832fd72b-0d51-44f5-97f0-7ad261d65e13
begin

	Pkg.add("DifferentialEquations")
	using DifferentialEquations
	Pkg.add("LinearAlgebra")
	using LinearAlgebra
	
end

# ╔═╡ 16dec684-3cbb-45e6-8761-16734643a2cb
using Markdown

# ╔═╡ 1a3f82ef-9ac3-48ae-abea-9b150b59a955
using InteractiveUtils

# ╔═╡ 584c650e-a819-4c98-809d-1101ec1e83f8
md"
>##### Gas Transport in a Porus Medium
>###### Scientific Computing, Winter Semester 2020/21 
>###### By Prof. Jürgen Fuhrmann  
>###### Group 08, Students: Marina Matthaiou, Maryam Asadi, Ria Pachal"

# ╔═╡ 1909662f-efa8-4f60-a7d3-33ccbff69897
md"
# Introduction
Gas transport in a porus medium equation is an eliptic parabolic equation, the problem is transient, so the solution is time dependant in every point of the domain. In this report it has been tried to elaborate the problem, introduce different possible solution in discritized mode and solve the problem with a finite volume package in Julia, named VoronoiFVM compiled by Prof. Fruhmann and co.
At first, the problem has been solved with Neumann boundary condition with two different initial values, first the penalty method, then the support of the Barenblatt solution in different time has been examined, and the initial values has been set with Barenblatt solution, moreover the Dirichlet boundary condition with initial value of the Barenblatt solution has been observed, and at the end the different output from two different solution method, namely VoronoiFVM and Differential Equation, has been compared.
All of the steps has been done in 1D and 2D space domain.
## Physics of the Problem
Gas transport in a porus medium has two part, adective and diffusive. The advection is generally analyzed by Darcy's law(Darcy, 1856), it states that gas Darcy's velocity $u_g$ is proportional to gass-phase presuure gradient $\nabla{P_g}$, and the gas phas permeability $K_g$, or we can write:

$\bar{u_g}=-\frac{k_g}{\mu_g}(\nabla{P_g}-\rho_g\bar{g})$
Where $\mu_g$ is the gas-phase viscosity and $g$ is the gravitational constant, which is usually being ignored due to its very low impact in the gas state, for the sake of simplicity we neglect the coefficient and rewrite the Darcy's law as:

$\bar{u_g}=-\nabla{P_g}$
## Porus equation for Gas Transport
$\partial_t u + \nabla(u^m)  = 0         \qquad{(1)}$
The entities describing the discrete system can be subdivided into two categories:
- geometrical data: $|\omega_k|, \gamma_k, \sigma_{kl}, h_{kl}$ together with the connectivity information of the triangles
- physical data: the number $m$ and the function $u,g$ describing the particular problem, where $g$ is a flux function approximating $\vec j$.
The solution of the nonlinear systems will be approximated with Voroinoi Finite Volume Method in 1D and 2D. At first let's see what is the flux function of this equation.

"

# ╔═╡ 98deb536-382e-42d9-8809-d2bcda17420a
md"
# Finite Volume Discretization approach 
To solve partial differential equations like we have in (Eq. 1), we form a matrices which evolve from the discretization of these PDEs.
We can approximate continuous functions by piecewise linear functions defined by the values $f_i= f(x_i)$. Using more points yields a better approximation. in order to do so we have to subdivide our domain into finite number of closed subsets. 

Fix a number of discretization points N
- Let $h=\frac{l}{N-1}$ in which $l$ is the ,length of the domain.
- Let $x_i=l(i-1)h\;$ $i=1\dots N$ be discretization points, 

this yeilds to finite difference approximation of derivitives:

$u'(x_{i+\frac12})\approx \frac{u_{i+1}-u_{i}}{h}$ and Same approach for second derivative:

$u''(x_i)=\frac{u'(x_{i+\frac12})-u'(x_{i-\frac12})}{h}$ 
This derivitives can be generalized into definition of the divergence, when the function has more than one variable, like $u(x_1, x_2, .., x_d)$ in $R^d$ domain.
 
This procedure can be implemented in the boundaries by the help of mirror values $u_0$ and $u_{N+1}$. This forms an $n\times n$ matrix.

For 2D We use 2D regular discretization $n\times n$ grid with grid points $x_{ij}=(l_x(i-1)h, l_y(j-1)h)$. This forms an $n^2\times n^2$ matrix.
In order to solve these matrices, there are several iterative methods plus preconditioning like: Incomplete LU decomposition, Gauss-Seidel, Jacobi, Generalized minimal residual method(GMRES), which are implemented in the Julia packages like:  [IterativeSolvers.jl](https://juliamath.github.io/IterativeSolvers.jl), [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl), [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl). We have to be careful that this standard PDE calculus is valid only if the subset is [Lipschitz] (https://en.wikipedia.org/wiki/Lipschitz_domain).

"


# ╔═╡ a382d48f-835e-470f-8cf4-a8ac35a9d45e
md" ## Discretization in Space
### Triangulation
Above mentioned grid forms a rectangular mesh within a rectangular domain, but it is always not the case and we might need to use other mesh generation methods, like triangles or quadrilaterals for 2D domain, the boundary of a polygonal domain $\Omega$  which is $\partial\Omega=\Gamma$ consists of finite subset of hyper planes in $\mathbb R^n$ (line segments for 2D domain). One of the useful meshes is simplex grid which has been used in this project. A condition for generating meshes is that they should be admissible.
Voronoi diagram is a method for generating admissible grids in a way that the distance of the point inside domain is being minimized at the Voronoi point sets, which means that it subdivide the region into the set of nearest neighbors which are convex set (dividing into half planes).

Delaunay triangulation connect the points that share an edge in Voronoi diagram, in a way that the triangle containing 3 point set will be on a circle, and no 4 point can be on its circumference. So the boundary of the triangle is containd inside the circle and the edge points(vertices) are placed on the circumference of this circle.

In [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl) package we can creat the the Voronoi diagram, Delaunay triangulation and the boundary of a given point set.
Below, 5 random points have been created with specifying the the [commond line switches](https://www.cs.cmu.edu/~quake/triangle.switch.html), we see the Voronoi diagrams, Delaunay triangulation as well as the boundary of the point set. 

we can see that Voronoi diagrams are prependicular to the boundary edges and every three point is placed on a circle containing the Delauny triangles. The minimum triagnle we can get is three in this case (it depends on the edges that Voronoi diagrams share, so on the distance of the points) and with different points this number will change (the Delauny triangulation of given point set is unique).

"

# ╔═╡ ed6f83f7-3db4-4a0e-b867-0cace9aa444e
md"In triangulate, we can also specify the constrained edge and the maximum area that the triangles can cover, in this case the triangles are not Delauny anymore, we can also specify the region that we want to refine the meshes, like close to the era that we need more accurate result due to the high error or complexity of the domain, it is also possible to divide the domain into different region with different mesh size."

# ╔═╡ 0c59788a-6e6b-4a13-b305-3ceee53fd8f0
md"
## Constructing Discrete System of Equations
After mesh generation, we can write the governing equation of the functions and partial differential equations at REVs: $w  \subset \Omega$ (Representative Elementary  Volume).

 $u(\vec x,t)$is a time dependant amount of a certain species on subset domain $\Omega$ within time interval $(o,T]$, $f(\vec x,t)$ is the function for species source in same domain and time interval and $\vec j(\vec x,t)$ is the vector of species flux.

We can calculate the flux $j(t)$going through $\partial w$ (the boundary of REV) by integrating flux function over its boundary, the amount of species and source function $u(t)$ and $f(t)$ respectively at $w$ also can be calculated by integration over $w$.
Conservation of the species at $w$ leads to form a continuity equation as below:

$u(t_1)-u(t_0)+ \int_{t_0}^{t_1} j(t).dt=\int_{t_0}^{t_1}  f(t).dt$
using Gauss theoreom we can rewrite this equation in the form of partial differential equation as below:

$\partial_t u(\vec x,t) + \nabla.\vec j(\vec x,t)=f(\vec x , t) \qquad{(2)}$
Flux $\vec j(\vec x,t)$ is proportional to the change of $u(\vec x,t)$ but in negative direction, which means we can write it as $\vec j(\vec x,t)=-\delta\nabla.\vec u(\vec x,t)$. Here $\delta$ can be constant, or space depndant, or even dependant on $u$, in the case of non-constant $\delta$ the equation becomes nonlinear.
If we rewrite the equation 2 with this definition, it leads to a second order PDE, namely parabolic PDE for transient problems:

$\partial_t u(\vec x,t) - \nabla.(\delta\vec\nabla u(\vec x,t))=f(\vec x , t) \qquad{(3)}$
Note that if there is a convective term, the flux can be written as:  $\quad\vec j=-\delta\vec\nabla u+ u\vec v$ where $\vec v$ is convective velocity.
In order to solve such equations, we need initial value: $u(x,t_0)$ as well as boundary conditions: governing equations on $\partial \Omega$.
"

# ╔═╡ fc300c0e-12db-4ce3-a26e-0b27591f23a1
md"
### Boundary Condition
Note that the equation (1) is satisfied in the interior of its defined domain, so the values at the boundary have to be determined by certian charecteristics, there are three different boundary condition:

- Dirichlet $\quad u(x) = g_i(x) \qquad\forall \vec x \in \Gamma_i$  (fixed solution at boundary)
- Neumann: $\quad\delta\nabla u(\vec{x},t).\vec{n}=g_i(\vec{x},t)\qquad \forall \vec x \in \Gamma_i$ (fixed boundary flux)
- Robin: $\quad\delta\nabla u(\vec{x},t).\vec{n}+\alpha_i(\vec{x},t)u(\vec{x},t)=g_i(\vec{x},t) \qquad \forall \vec x \in \Gamma_i$ (Boundary flux proportional to solution)
Where $\Gamma_i$s are finite non-itersecting subsets $\partial\Omega$.
"

# ╔═╡ 744983d7-eb8d-4cfe-8d71-06f631e4b13e
md"
## Constructing Contorl Volumes
Now we have everything we need to implement the FVM to our domain $\Omega$ which is subdivided into finite REVs: $w_i$, we have to assign a value of $u$ at each $w$: $u_i$, this value can be calculated at collocation points $x_i \in w$,
find PDE approximation as explained above, assume $\Gamma_i$ are planer that subdivide the boundary of polygonal domain $\Omega$ such that: $\vec n_i \perp \Gamma_i$

Picture below from [Juliahub](The Voronoi finite volume method) shows a schematic view of the finite volume approach.
![](https://juliahub.com/docs/VoronoiFVM/XOFEJ/0.8.5/vor.png)
Note that the domains are not intersecting, the collocation points of the neighbor REVs build a  line which direction is a vector: $\vec \nu_{kl}$ normal to the boundary: $\sigma_{kl}$ (admissibility of grids). This allows us to aproximate the $\nabla u. \vec \nu_{kl} \approx \frac{u_l-u_k}{|x_l-x_k|}$.
Whenever the REV $w_k$ has mutual edge with $\partial \Omega$ at $\Gamma_m$ its called $\gamma_km$, note that in this case collocation point $x_m \in \partial \Omega$. This allows us to assign boundary values at this collocation points. The set of $w_k$'s neighbor and the boundary part of $w_k$ are called $\mathcal N_l$ and $\mathcal G_m$ respectively, note that the union of this two sets build the boundary of $w_k$.
"

# ╔═╡ 04cf9a64-3c8d-426e-bae3-93867508d978
md"
In 1D we can construct $w_k$ as following:
$w_k=\left\{ \begin{array}{rcl}
(x_1,\frac{x_1+x_2}{2})
   \qquad \qquad \qquad \quad k=1 \\(\frac{x_{k-1}+x_k}{2}, \frac{x_k+x_{k+1}}{2}) \qquad \quad 1<k<n \\
(\frac{x_{n-1}+x_n}{2}, x_n) \qquad \qquad \quad \quad k=n
\end{array}\right.$
In 2D we construct $w_k^x$ and $w_k^y$ as above, then collocation points will be: $\vec x_{kl}=(x_k,yl)$ and $w_{kl}=w_k^x \times w_k^y$, but again this gives rectangular meshes, if we construct the delauny triangulation on collocation points $\vec x_k$ and restricted Voroni cells $w_k$ such that $x_k \in w_k$, the corner of Voronoi cells are the centre of circle connecting midpoint of the boundary edges, so the admissibility condition is fulfilled naturally (boundary edge prependicular to line connecting collocation points). The trangulation edges are line connecting $x_k$ and $x_l$ and the trianhulation nodes are collocation points: $x_k$.
"

# ╔═╡ 13f173f3-15ea-4110-861a-4aa663a506ca
md"
To solve equation 3 we need to integrate over control volume $w_k$ and approximate the flux between neighbouring control volumes, moreover the boundary condition should be specified, in this project the Dirichlet boundary condition can describe the problem very well, because it has a fixed function at the boundaries, also the Barenblatt equation which is a exact solution of the problem can help to define this function at boundary. Initial values should be introduced into the FVM method, to start the iteration scheme for time dependant value $\vec u(\vec x,t)$.

$\nabla u. \vec \nu_{kl} \approx \frac{u_l-u_k}{|h_{kl}|} \qquad\text{:approximation of normal derivitive.}$ 
$\int_{\sigma_{kl}} \vec j .\vec \nu_{kl}.ds \approx \frac{|\sigma_{kl}|}{|h_{kl}|}\delta (u_k-u_l) \colon= \frac{|\sigma_{kl}|}{|h_{kl}|}g(u_k,u_l) \quad\text{:approximation of flux between REV's}$ 
$\int_{\gamma{km}} \vec j .\vec n_m.ds \approx |\gamma{km}|g(x_m) \qquad \text{:approximation of flux at boundary for Neumann BC.}$ 
$f_k=\frac1{w_k}\int_{w_k}f(\vec x)dw \qquad \text{:approximation of source function(RHS)}$
where g(.,.) is flux function and $h_{kl}=|x_k-x_l|$
The above equations can form a discrete system of equation with $N= |\mathcal N|$ equation at each REV, and $N$ unknowns for each collocation points.
### Matrix Assembely Algorithm:
the discrete system can be written based on triangulation an matrix will be assembled in two loops, first loop over triangles and calculate its contribution to the matrix entries, then loop over the boundaries and calculate their contribution, then for the RHS the loop over the triangles will add up their contribution to the source function. The solution of this matrix gives a piecewise solution at each REV: $w_k$ and $x_k$ are the pieacewise linear function on triangles.
"

# ╔═╡ e1777730-babd-42f6-8022-36ff9e6e7d6e
md"## Nonlinear Problems
Finite volume for nonlinear problems follows the same rules, but as the flux is proportional to $-\vec \nabla u$ the approximation of flux function: $\quad  g(u_k,u_l) \approx D(u)(u_k-u_l)$ can be done in two different way: 
- integrating: $\mathcal D(u)=\int_{u_k}^{u_l}D(\xi)d \xi$
- averaging: $\quad g(u_k,u_l)=D(\frac{u_k+u_l}2)(u_k-u_l)$

Equation 2 is a general form of equation 1 in which the RHS is zero, Let $\vec j= -D(u)\vec\nabla u$  so it gives: $D(u)=mu^{m-1}$ 
and by integration the flux will be:

$\mathcal D(u)= \int_0^u D(\xi) d\xi = u^m$
We can define flux function between $u_k$ and $u_l$ as: $g(u_k, u_l)=u_k^m - u_l^m$.
This discrete system of the nonlinear equation gives an $N \times N$ matrix with N equation and N unknown, note that the Jacobi matrix of this system should be written as well to find the residula at each step as well as the derivitive of matrix assemply (we can use Julia dual number to do so).
"

# ╔═╡ af2facd0-ba67-4821-b83b-2f1e09157770
md"
# Time Discretization in Transient Problems
Now that we have the method to construct the FVM in nonlinear problems, we need to define a way to descretize time for the transient problem, so we need the boundary condition, flux function and initial value at boundaries.
There two different way to do so:
- first discretize in time then space (Rothe method)
- first discretize in space (get huge ODE) then decretize time (method of lines)
We have to choose the end time, and refinement of our time discretization: $t_0=t^0<t^1<...< t^N=T \quad$ and let: $\quad \tau^n=t^n-t^{n-1}$ for $i=1, .., N$ solve equation 3 as:

$\frac{u^n-u^{n-1}}{\tau^n}- \nabla.D\vec \nabla u^{\theta}=f \quad \text{in}\quad  \Omega \times [t_0,T]$
$D\vec \nabla u^{\theta}.\vec n= g \quad \text{on}\quad \partial \Omega \times [t_0,T]$ 
Where $u^{\theta}=\theta u^n +(1-\theta)u^{n-1}$
- $ \theta=1\quad $ Implicit Euler method (backward). Solve PDE problem in each timestep. First order accuracy in time
- $ \theta=0\quad $ Explicit Euler method (forward). First order accuracy in time, does not contain PDE 
- $ \theta=\frac1{2}$ Crank-Nicolson scheme. Solve PDE problem in each timestep. Second order accuracy in time
The problem for explicit euler method is that it approximate the solution to an ODE by using first derivitive which is a poor approximation and can add up the error in itteration steps, the only solution to this is refining the discretization time until the stability of the results, which is not an easy task for all of the problems.

Search function for iteration should be defined such that: $u(x,o)=u_0(x)\quad\text{in }\Omega$, integratation of equation over space-time control volume $w_k\times (t^{n-1},t^n)$ and dividing by $\tau^n$ will result to the matrix assembly.

The package used here to construct FVM and implement the grid generation and implementing the discrete system assembly, [VoronoiFVM.jl](https://j-fu.github.io/VoronoiFVM.jl/stable/) uses implicit Euler method plus damped Newton method to solve the problem.
"

# ╔═╡ 408a53a3-beaa-497c-aa18-c1c2601eda4f
md"
# Solution Methods for Discrete Systems
For problems like :$A(u)=f(u)$ where $A$ is anonlinear operator, there is not a Guassian elimination schem, so the problem has to be solved itteratively, there are several schem to do so, some of them are explained below:
- **Fixpoint iteration scheme**, $A(u)=M(u)u$, where $M(u)$ is linear at each u,choose an initial value $u_0$ and at each iteration step, solve: 
$M(u^i)u^{i+1}=f$ terminate if the residual is small enough to reach the desired accuracy. This method has large domain of convergence that may be slow. 
- **Newton iteration scheme**, $A'(u)=(a_{kl})$ where $A'(u)$ is the Jacobi matrix of first partial derivatives of $A$ at point $u$, and $a_{kl}= \frac{\partial}{\partial u_l}A_k(u_1\dots u_n)$

In the $i$-th iteration step:

$u_{i+1}=u_i - (A'(u_i))^{-1}(A(u_i) -f)$

Calculate residual: $r_i=A(u_i)-f$, solve linear system for update: $A'(u_i)h_i = r_i$, update solution: $u_{i+1}=u_i-h_i$

This method has potenially small domain of convergence (good initial value is necessary), initial convergence is possibly slow, but quadratic convergence are close to the solution, has linear and quadratic convergence, uses automatic differentiation for Newton's method, damped Newton iteration speeds the solution, parameter embedding improves it.
The dual numbers in Julia can be used  to find the derivitives used in PDE's, this is implemented in forwardiff package.
- **Damped Newton iteration**, do not use the full update, but damp it by some factor which we increase during the iteration process, linesearch also can be used (automatic detection of a damping factor).
- **Parameter Embedding** is another way to solve discrete system in parameter dependant problems. If we use parameter embedding along with damping and update based convergence control, we can solve even more difficult nonlinear problems.
VoronoiFVM package solve stationary problem via parameter embedding, and uses implicit Euler method + damped Newton's method to solve time dependent problem. Time step control is performed according to the data in control.
"


# ╔═╡ e3e56550-825e-4905-bcb8-3c8086671c28
md"
# Simulation results

Grids 1D domain $\Omega=(-1,1)$  consisting of N=$(@bind N Scrubbable(10:10:100,default=20)) points at x axis 

Grid 2D domain $\Omega=(-1,1)×(-1,1)$ at x, y axses with N point.
"

# ╔═╡ 1ce4a24c-2be3-4b9b-9045-72b90f329df1
md"
## Neumann Boundary Condition Solution
The VoronoiFVM has been used in order to solve the diffusion problem with neumann boundary condition, in order to do this, a function has been introduced to generate the initial values at boundary. The function is defined as bellow for 1D and 2D domain respectively:

$F_{peak}(x)=exp(-100×x^2)$

$F_{peak}(x,y)=exp(-100(x^2+y^2))$
### 1D Domain $\Omega=(-1,1)$"


# ╔═╡ 14730986-22e5-4c63-a076-2a84686f3edf
# Define function for initial value $u_0$ with two methods - for 1D and 2D problems

begin

    fpeak(x)=exp(-100*x^2)

    fpeak(x,y)=exp(-100*(x^2+y^2))
end;

# ╔═╡ eb3b35da-199f-4e5a-b294-a3f428f53535
md" In order to estimate the accuracy of the result, the exact Barenbllatt solution has been used to find the error. It is clear that this solution lead to a high error, this might be a result of the bad initial value, or the time interval in which the itteration took place. Also with change of time, the error will not change significantly.
### 2D Domain $\Omega=(-1,1)×(-1,1)$
Same approach is applied here, we can see that the results are also very far from the exact solution in 2D as well.
"

# ╔═╡ 3d56289e-c1ad-4751-be54-14e35f88dd5a
md""" 
## Barenblatt solution support
we saw that initializaing the solution with a small value function $F_{peak}$ did not lead to an accurate solution. Now in this section we want to use a better initial values for Neumann BC, which is  Barenblatt solution, but first we need to know the time interval in which the solution supports the domain.

For space dimension $d$ in the domain $R_d × (0, ∞)$ the equation has a radially symmetrict exact solution, the so-called
Barenblatt solution:


$u(x,t)= max(0,t^{-\alpha}(1-\frac{\alpha(m-1)r^{2}}{2dmt^{\frac{2\alpha}{d}}})^{\frac{1}{m-1}})$ 

Where:
$r=|x|$ and $\alpha=\frac{1}{m-1+\frac{2}{d}}$. 

This solution has finite support and spreads a finite amount of mass over the space domain
In the this section, it has been tried to find the supporting of the solution in different times at boundaries of the domain
value of m, m=$(@bind m Scrubbable(1:1:10,default=2)) points.

In order to find the support of the solution,  we widened our domain from $(-1,1)$ to $(-2,2)$ and in 2D from $(-1,1)×(-1,1)$  to $(-2,2)×(-2,2)$ to get an impression at which time the Barenblatt will leave the domain at 2D . The values on the tables shows whenevr the value of u is greater than zero at marginal (boundary) points of the domains. this can be seen in the plots as well, by increasing time, we can find the $t_1$ so that Barenblatt supports the solution at $(t_0,t_1)$
###### 1D Domain $\Omega=(-1,1)$.
"""


# ╔═╡ 5b923416-0c28-4a10-8d0b-b4cd62e5ffba
function barenblatt1d(x,t,m,dim=1)
	
	alpha=(1/(m-1+(2/dim)))
	K=(alpha*(m-1)*x^2)/(2*dim*m*(t^(2*alpha/dim)))
	
	B=(t^(-alpha))*((1-K)^(1/(m-1)));
	
    if B<0.0
        B=0.0
    end
	
    return B
end;

# ╔═╡ 6e1eb8c9-9e6b-4f5a-b01b-0cca7b10f80e
x=collect(range(-2,2,length=2N+1));

# ╔═╡ 2408a0e9-532c-4be5-ab38-61fb5437aba1
begin 
	t_s=[10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),1]
	u_L=zeros(length(t_s))
	u_R=zeros(length(t_s))
		for i=1:length(t_s)
		u_L[i]=barenblatt1d(-1,t_s[i],m)
		u_R[i]=barenblatt1d(1,t_s[i],m)	
		u_L, u_R
		end
end;

# ╔═╡ 8bc22a4f-0def-4d51-b9bf-310a40ea8482
DataFrame(time= t_s[:], LeftValue = u_L[1:length(t_s)], RightValue =u_R[1:length(t_s)])

# ╔═╡ ba9f286c-b6ab-4986-976b-761b6c53d26a
md"""
time=$(@bind t Slider(10^(-5):0.001:1,default=0.01,show_value=true))
"""

# ╔═╡ 0f2f3926-34ad-4660-8c13-35dee4ec6474
scalarplot(x,map(x->barenblatt1d(x,t,m),x),Plotter=PyPlot,resolution=(1200,600))

# ╔═╡ bf1a6b3a-f87f-470e-9074-47f475d760f8
md""" 
###### 2D Domain $\Omega=(-1,1)×(-1,1)$
"""

# ╔═╡ 998fb1e8-ff11-46a7-bb41-98bbb3eb76f5
function barenblatt2d(x,y,t,m,dim=2)
	
	alpha=(1/(m-1+(2/dim)))
	r=sqrt(x^2+y^2)
	K=(alpha*(m-1)*r^2)/(2*dim*m*(t^(2*alpha/dim)))
	
	B=(t^(-alpha))*((1-K)^(1/(m-1)));
	
    if B<0.0
        B=0.0
    end
	
    return B
end;

# ╔═╡ 63a41d85-79ad-49b1-aefd-c7a4fd250cea
y=collect(range(-2,2,length=2N+1));

# ╔═╡ e6f3a74c-9161-4b83-bd5d-61028b4e81dd
grid2d=simplexgrid(x,y);

# ╔═╡ 041a1492-0aaf-4187-b7ed-cd49b092b727
begin 
	t_s2d=[10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),1]
	u_L2d=zeros(length(t_s2d))
	u_R2d=zeros(length(t_s2d))
	u_U2d=zeros(length(t_s2d))
	u_D2d=zeros(length(t_s2d))
		for i=1:length(t_s2d)
		u_L2d[i]=barenblatt2d(-1,0,t_s2d[i],m)
		u_R2d[i]=barenblatt2d(1,0,t_s2d[i],m)
		u_U2d[i]=barenblatt2d(0,-1,t_s2d[i],m)
		u_D2d[i]=barenblatt2d(0,1,t_s2d[i],m)
		
		end
end;

# ╔═╡ 988f0e80-bf87-4c78-bba2-a04d53db519c
DataFrame(time= t_s2d[:], LeftValue = u_L2d[1:length(t_s2d)], RightValue=u_R2d[1:length(t_s2d)], UpperValue=u_U2d[1:length(t_s2d)],BottomValue=u_D2d[1:length(t_s2d)])

# ╔═╡ d9300a5a-0682-4f72-862f-ca90215b553f
scalarplot(grid2d, map((x,y)->barenblatt2d(x,y,t,m), grid2d),Plotter=PyPlot,resolution=(1200,600))

# ╔═╡ 36f4fd5c-d093-4eb9-8956-4dc4af0d34f7
md" We can see that Barenblatt supports the domain till around 0.1 $\sec$ for 1D and around 0.01 $\sec$ for 2D,
now we can set the function to solve the PME with Neumann BC, within this domain and initial value specified by Barenblatt solution." 

# ╔═╡ 98f441ac-f7fa-4979-b0b6-c4ced36f84b6
md"""
## Initial value u(x, t0) = b(x, t0)
### 1D Domain $\Omega=(-1,1)$
"""

# ╔═╡ b037aa3b-c4d3-44c5-9b75-a075657fdd81
md"""
### 2D Domain $\Omega=(-1,1)×(-1,1)$
"""

# ╔═╡ 24a8924f-b9ef-4cce-bf11-d0803f254b95
md"Here a 3D plot of the solution in different time steps is shown so that the impression of the solution is more clear."

# ╔═╡ 276ebe5b-7799-4a41-901c-931972a65e95
md"2D plots are also shown to compare with other methods."

# ╔═╡ a03c19d7-2b7b-456b-8dbe-db2809221f2d
md"As we can see the solution has significantly improved by fixing the Neumann BC initial values adopted from Barenblatt solution within desired time intervals"

# ╔═╡ 944b1ade-abd3-49f8-91d9-b01e8254911e
md"""

## Dirichlet boundary Condition
We can still improve the solution if we use the dirichlet boundary condition, because the equation (1) shows that the value of the species at boundary can be defined by a function, rather than having a fixed flux function at boundaries, in this way it is not necessary to find the initial value randomly like we did at first for Neumann because if we can find the function at boundary, we can easily set the time to $t_0$ and find the initial value. Below, the solution has been shown for the Dirichlet boundary condition, in 1D the initial values are taken from the Barenblatt at left and right of the domain, in 2D the initial values have been taken from the midpoint of each edge of the boundary.
### 1D Domain $\Omega=(-1,1)$
"""

# ╔═╡ 5ab09487-01bf-437d-8e45-53f97fbc9ff0
md"""
##### 2D Domain $\Omega=(-1,1)×(-1,1)$
"""

# ╔═╡ f7fd554e-3885-4c58-afa8-26b00dba8b96
md"
The solution has been improved impressively, which means that initial values that are following the Dirichlet BC are a better choice for this problem.
There are still several method to increase the accuracy and speed of the solution, some of them includes refining mesh when we are moving from edges to the centre of domain, because we get higher error at these points, using [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) will give a better accuracy as well as faster solution, we can see the result of this comparison as following:
![](https://github.com/idasamayram/Scientific-computing-project/blob/main/comparison.JPG)
"

# ╔═╡ 35a21196-41c4-4ed7-8eae-8d559bd5500e


# ╔═╡ 84ad06f5-72c9-456a-aefe-f4fff86ce60d


# ╔═╡ c46db835-f4f0-4c09-9e16-b093d73c8312
import Base:push!

# ╔═╡ c14c98b0-0fb0-4572-a721-a9cffba63e3e
function example_convex_hull_voronoi_delaunay(;n=10,circumcircles=false)
    triin=Triangulate.TriangulateIO()
    triin.pointlist=rand(Cdouble,2,n)
    (triout, vorout)=triangulate("vcDQ", triin)
    plot_in_out(PyPlot,triin,triout,voronoi=vorout,circumcircles=circumcircles)
end;

# ╔═╡ f13a6d3b-546b-48c2-8190-97644bc17837
example_convex_hull_voronoi_delaunay(;n=5,circumcircles=true);gcf()

# ╔═╡ 1e0c8cbb-58b2-4019-9e30-1589f0ebf92f
# Create discretization grid in 1D or 2D with approximately n nodes
function create_grid(n,dim)
	nx=n
	if dim==2
		nx=ceil(sqrt(n))
	end
	X=collect(-1:1.0/nx:1)
	if dim==1
      grid=simplexgrid(X)
	else
      grid=simplexgrid(X,X)
	end
end;

# ╔═╡ d2c895f0-4b44-4260-a2c2-49b06b43f063
gridplot(create_grid(N, 1),Plotter=PyPlot,resolution=(600,200))

# ╔═╡ 9c9dcd8a-e79a-45a4-be26-7a5f40348f16
gridplot(create_grid(N, 2),Plotter=PyPlot,resolution=(600,200))

# ╔═╡ 7f4b9a30-ad19-4e21-8bf1-a443f503ca05

   function diffusion(;n=100,dim=1,tstep=1.0e-4,tend=0.1, dtgrowth=1.1)
	grid=create_grid(n,dim)
	
	## Diffusion flux between neigboring control volumes
	function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
	end
	

    ## Storage term (under time derivative)
	function storage!(f,u,node)
		f[1]=u[1]
	end
	
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
    
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])
    
    ## Create a solution array
    inival=unknowns(sys)
	t0=0.001

	## Broadcast the initial value
    inival[1,:].=map(fpeak,grid)

	control=VoronoiFVM.NewtonControl()
	control.Δt_min=tstep
	control.Δt=tstep
	control.Δt_max=tend
	control.Δu_opt=1
	control.Δt_grow=dtgrowth
	
	tsol=solve(inival,sys,[t0,tend];control=control)
	return grid,tsol
end;

# ╔═╡ 9a384177-ab71-4b9e-9a04-263fddb68718
grid_diffusion,tsol_diffusion=diffusion(dim=1,n=N);

# ╔═╡ 7ee60df1-a27f-476b-a40e-9cd1a942d0ba
md"""
Solution at different time steps: Controlling time step number:$(@bind t_diffusion Slider(1:length(tsol_diffusion),default=10,show_value=true))
"""

# ╔═╡ d0f8c2ce-52dd-4451-986a-cea9fd3bd35e
begin
	np=GridVisualizer(Plotter=PyPlot,layout=(3,1),fast=true,resolution=(600,600))
    	
		
        	scalarplot!(np[1,1],grid_diffusion,tsol_diffusion[1,:,t_diffusion],title=@sprintf("Numerical Solution, Neumann BC, t=%.3g",tsol_diffusion.t[t_diffusion]))
		
        	scalarplot!(np[2,1],grid_diffusion,map(x->barenblatt1d(x,tsol_diffusion.t[t_diffusion],m),grid_diffusion),title=@sprintf("Exact Barenblatt Solution, t=%.3g",tsol_diffusion.t[t_diffusion]))
		
		scalarplot!(np[3,1],grid_diffusion,map(x->barenblatt1d(x,tsol_diffusion.t[t_diffusion],m),grid_diffusion)-tsol_diffusion[1,:,t_diffusion],title=@sprintf("Error of the Numerical Solution, t=%.3g",tsol_diffusion.t[t_diffusion]))
	
        reveal(np)
    	
	gcf()
end

# ╔═╡ 6b373f2e-61e9-4219-aeb4-6f86cc7a2189
function diffusion2D(;n=100,dim=2,tstep=1.0e-4,tend=1, dtgrowth=1.1)
	grid=create_grid(n,dim)
	
	## Diffusion flux between neigboring control volumes
	function flux!(f,u0,edge)
		
		u=unknowns(edge,u0)
		
        
        f[1]=((u[1,1])^m-(u[1,2])^m)
	end
	

    ## Storage term (under time derivative)
	function storage!(f,u,node)
		f[1]=u[1]
	end
	
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
    
    sys2D=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys2D,1,[1])
    
    ## Create a solution array
    inival=unknowns(sys2D)
	t0=0.001

	
	
    inival[1,:].=map(fpeak,grid)
	
	

	control=VoronoiFVM.NewtonControl()
	control.Δt_min=0.01*tstep
	control.Δt=tstep
	control.Δt_max=0.1*tend
	control.Δu_opt=0.1
	control.Δt_grow=dtgrowth
	
	tsol2D=solve(inival,sys2D,[t0,tend];control=control)
	return grid,tsol2D
	
end;

# ╔═╡ ae1adee2-cb6f-471c-a9b7-8c67ffd2f844
grid_diffusion2D,tsol_diffusion2D=diffusion2D(dim=2,n=N);

# ╔═╡ 6207ed96-c2d6-4fc3-8a71-29f6f0101b61
md"""
Solution at different time steps, Controlling time step number:$(@bind t_diffusion2D Slider(1:length(tsol_diffusion2D),default=10,show_value=true))
"""

# ╔═╡ 8ea3bb12-d0f2-487d-a8ad-3d5030a4ffa7
function porous_diffusion(;n=100,dim=1,tstep=1.0e-4,tend=0.1, dtgrowth=1.1)
	grid=create_grid(n,dim)

		
	#flux 
	function flux!(f,u0,edge)  
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
	end
		
	#storage
	function storage!(f,u,node)  
		f[1]=u[1]
	end
		
	
	
		
	#create physics	
	physics=VoronoiFVM.Physics(
	flux=flux!,
	storage=storage!)   
		
	#create system	
	sys=VoronoiFVM.DenseSystem(grid,physics)   
	enable_species!(sys,1,[1])
	
	#initial value	
	inival=unknowns(sys)
	t0=0.001
		
	# Broadcast the initial value
	inival[1,:].=map(x->barenblatt1d(x,t0,m),grid)
		
	# Create solver control info for constant time step size	
	control=VoronoiFVM.NewtonControl()
	control.Δt_min=0.01*tstep
	control.Δt=tstep
	control.Δt_max=0.1*tend
	control.Δu_opt=1
	control.Δt_grow=dtgrowth
	
	
	tsol=solve(inival,sys,[t0,tend],control=control)
	
	return [grid,tsol]
end;

# ╔═╡ 46b92be8-c262-4c80-83c7-ae3b12d01c71
grid_b,tsol_b=porous_diffusion(dim=1,n=N);

# ╔═╡ 0a00fa12-0e55-44f3-9685-e5d1cd7adb52
md"""
Solution at different time steps, time: $(@bind t_b Slider(1:length(tsol_b),default=10,show_value=true))
"""

# ╔═╡ 2dd023d3-f1b6-404f-b3f0-d8997f121ee2
begin
	
	rp=GridVisualizer(Plotter=PyPlot,layout=(3,1),fast=true,resolution=(600,600))
    	
		
        	scalarplot!(rp[1,1],grid_b,tsol_b[1,:,t_b],title=@sprintf("Numerical solution, t=%.3g",tsol_b.t[t_b]))
		
        	scalarplot!(rp[2,1],grid_b,map(x->barenblatt1d(x,tsol_b.t[t_b],m),grid_b),title=@sprintf("Exact barenblatt, t=%.3g",tsol_b.t[t_b]))
		
		scalarplot!(rp[3,1],grid_b,map(x->barenblatt1d(x,tsol_b.t[t_b],m),grid_b)-tsol_b[1,:,t_b],title=@sprintf("error 1D, t=%.3g",tsol_b.t[t_b]))
	
        reveal(rp)
    	
	gcf()
end

# ╔═╡ 034beddc-3a6b-4ea5-bcf7-73718adf80b8
function porous_diffusion2D(;n=100,dim=2,tstep=1.0e-4,tend=0.01, dtgrowth=1.1)
	grid=create_grid(n,dim)

		
	#flux 
	function flux!(f,u0,edge)  
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
	end
		
	#storage
	function storage!(f,u,node)  
		f[1]=u[1]
	end
		
	
	
		
	#create physics	
	physics=VoronoiFVM.Physics(
	flux=flux!,
	storage=storage!)   
		
	#create system	
	sys2D=VoronoiFVM.DenseSystem(grid,physics)   
	enable_species!(sys2D,1,[1])
	
	#initial value	
	inival=unknowns(sys2D)
	t0=0.001
	
   
	# Broadcast the initial value
	inival[1,:].=map((x,y)->barenblatt2d(x,y,t0,m),grid)
		
	# Create solver control info for constant time step size	
	control=VoronoiFVM.NewtonControl()
	control.Δt_min=0.01*tstep
	control.Δt=tstep
	control.Δt_max=0.1*tend
	control.Δu_opt=1
	control.Δt_grow=dtgrowth
	
	
	tsol=solve(inival,sys2D,[t0,tend],control=control)
	
	return [grid,tsol[1:end]]
end;

# ╔═╡ a62b283f-c407-4c8f-911d-486a99b07bec
grid_b2d,tsol_b2d=porous_diffusion2D(dim=2,n=N);

# ╔═╡ 190d69d8-8980-4ccc-b978-73256cd56de5
begin
	
pn2=GridVisualizer(Plotter=PyPlot,layout=(3,1),fast=true,resolution=(600,600))
    	
		
scalarplot!(pn2[1,1],grid_diffusion2D,tsol_diffusion2D[1,:,t_diffusion2D],title=@sprintf("Numerical Solution, Neumann BC, t=%.3g",tsol_diffusion2D.t[t_diffusion2D]))
		
scalarplot!(pn2[2,1],grid_diffusion2D,map((x,y)->barenblatt2d(x,y,tsol_diffusion2D.t[t_diffusion2D],m),grid_diffusion2D),title=@sprintf("Exact Barenblatt Solution, t=%.3g",tsol_diffusion2D.t[t_diffusion2D]))
		
		scalarplot!(pn2[3,1], grid_b2d,map((x,y)->barenblatt2d(x,y,tsol_diffusion2D.t[t_diffusion2D],m),grid_diffusion2D)-tsol_diffusion2D[1,:,t_diffusion2D],title=@sprintf("Error of the Numerical Solution, t=%.3g",tsol_diffusion2D.t[t_diffusion2D]))
	
        reveal(pn2)
    	
	gcf()
end

# ╔═╡ 978aca9c-7e19-485d-a943-45e72d4b83a6
md"""
Time step number:$(@bind t_b2d Slider(1:length(tsol_b2d),default=10,show_value=true))
"""

# ╔═╡ 0de74633-25d1-4c06-a38d-f2deb982cd2a
let
	clf()
	
	suptitle("Numerical Solution ")
	
	
	surf(grid_b2d[Coordinates][1,:],grid_b2d[Coordinates][2,:], tsol_b2d[1,:,t_b2d],cmap=:summer) # 3D surface plot
	ax=gca(projection="3d")  # Obtain 3D plot axes
	
surf(grid_b2d[Coordinates][1,:],grid_b2d[Coordinates][2,:], map((x,y)->barenblatt2d(x,y,tsol_b2d.t[t_b2d],m),grid_b2d),cmap=:autumn) # 3D surface plot
	ax=gca(projection="3d")
	
	xlabel("x")
	ylabel("y")
	gcf()
end

# ╔═╡ 7c1eb0a7-ee28-4fb5-b560-8040793202c8
let
	clf()
	
	suptitle("Exact Barenblatt Solution")
	surf(grid_b2d[Coordinates][1,:],grid_b2d[Coordinates][2,:], map((x,y)->barenblatt2d(x,y,tsol_b2d.t[t_b2d],m),grid_b2d),cmap=:autumn) # 3D surface plot
	ax=gca(projection="3d")
	xlabel("x")
	ylabel("y")
	gcf()
end

# ╔═╡ 86da8092-f771-4ef0-a6c0-4bcb45ec4173
let
	clf()
	
	suptitle("Error of the Numerical Solution")
	surf(grid_b2d[Coordinates][1,:],grid_b2d[Coordinates][2,:], map((x,y)->barenblatt2d(x,y,tsol_b2d.t[t_b2d],m),grid_b2d)-tsol_b2d[1,:,t_b2d],cmap=:autumn) # 3D surface plot
	ax=gca(projection="3d")
	xlabel("x")
	ylabel("y")
	gcf()
end

# ╔═╡ 084cb624-fdeb-42e6-9c1b-c5f03beaa32f
begin
	
	pd=GridVisualizer(Plotter=PyPlot,layout=(3,1),fast=true,resolution=(600,600))
    	
		
scalarplot!(pd[1,1],grid_b2d,tsol_b2d[1,:,t_b2d],title=@sprintf("Numerical Solution with Neumann Boundary Condition, t=%.3g",tsol_b2d.t[t_b2d]))
		
scalarplot!(pd[2,1],grid_b2d,map((x,y)->barenblatt2d(x,y,tsol_b2d.t[t_b2d],m),grid_b2d),title=@sprintf("Exact Barenblatt Solution, t=%.3g",tsol_b2d.t[t_b2d]))
		
		scalarplot!(pd[3,1], grid_b2d,map((x,y)->barenblatt2d(x,y,tsol_b2d.t[t_b2d],m),grid_b2d)-tsol_b2d[1,:,t_b2d],title=@sprintf("Error of the Numerical Solution in 2D, t=%.3g",tsol_b2d.t[t_b2d]))
	
        reveal(pd)
    	
	gcf()
end

# ╔═╡ 39afbd01-44bd-47b1-9ab1-aadd019bd08d
   function Dirichlet(;n=100,dim=1,tstep=1.0e-4,tend=0.01, dtgrowth=1.1,DiffEq=nothing)
	grid=create_grid(n,dim)
	
	## Diffusion flux between neigboring control volumes
	function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
	end
	

    ## Storage term (under time derivative)
	function storage!(f,u,node)
		f[1]=u[1]
	end
	
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
	t0=0.001
    ul=barenblatt1d(-1,t0,m)
	ur=barenblatt1d(1,t0,m)
	
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])
    boundary_dirichlet!(sys, 1, 1,ul) # Set left Dirichlet boundary conditions
	boundary_dirichlet!(sys, 1, 2,ur)
   
	## Create a solution array
  	inival=unknowns(sys) 
 

	


	## Broadcast the initial value
    inival[1,:].=map(x->barenblatt1d(x,t0,m),grid)

	control=VoronoiFVM.NewtonControl()
	control.Δt_min=0.01*tstep
	control.Δt=tstep
	control.Δt_max=0.1*tend
	control.Δu_opt=0.1
	control.Δt_grow=dtgrowth
	
	tsol=solve(inival,sys,[t0,tend];control=control)
	return grid,tsol
end;

# ╔═╡ e014b546-a855-486b-a5f9-cdda73611bc5
gridd,tsold=Dirichlet(dim=1,n=N);

# ╔═╡ 8931a3fb-43fd-4a38-b760-6e4d8ca6ad77
md"""
Time step number:$(@bind td Slider(1:length(tsold),default=10,show_value=true))
"""

# ╔═╡ bf505926-67b8-4071-8a86-f7dd87850bfc
begin
	
	pdb=GridVisualizer(Plotter=PyPlot,layout=(3,1),fast=true,resolution=(600,600))
    	
		
        	scalarplot!(pdb[1,1],gridd,tsold[1,:,td],title=@sprintf("Numerical solution with Dirichlet Boundary Condition, t=%.3g",tsold.t[td]))
		
        	scalarplot!(pdb[2,1],gridd,map(x->barenblatt1d(x,tsold.t[td],m),gridd),title=@sprintf("Exact Barenblatt Solution, t=%.3g",tsold.t[td]))
		
		scalarplot!(pdb[3,1],gridd,map(x->barenblatt1d(x,tsold.t[t_b],m),grid_b)-tsold[1,:,t_b],title=@sprintf("Error of the Numerical Solution in 1D, t=%.3g",tsold.t[td]))
	
        reveal(pdb)
    	
	gcf()
end

# ╔═╡ 658cb684-7300-4f75-ac86-db32a753aaa0
 function Dirichlet2D(;n=100,dim=1,tstep=1.0e-4,tend=1, dtgrowth=1.1,DiffEq=nothing)
	grid=create_grid(n,dim)
	
	## Diffusion flux between neigboring control volumes
	function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]^m-u[1,2]^m
	end
	

    ## Storage term (under time derivative)
	function storage!(f,u,node)
		f[1]=u[1]
	end
	
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
	t0=0.001
    ul=barenblatt2d(-1,0,t0,2)
	ur=barenblatt2d(1,0,t0,2)
	ud=barenblatt2d(0,-1,t0,2)
	uu=barenblatt2d(0,1,t0,2)
	
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])
    boundary_dirichlet!(sys, 1, 4,ul) # Set left Dirichlet boundary conditions
	#boundary_dirichlet!(sys, 1, 1,ud)
   	boundary_dirichlet!(sys, 1, 2,ur) # Set left Dirichlet boundary conditions
	#boundary_dirichlet!(sys, 1, 3,uu)
	
	## Create a solution array
  	inival=unknowns(sys) 
 

	


	## Broadcast the initial value
 
 	inival[1,:].=map((x,y)->barenblatt2d(x,y,t0,m),grid)
	control=VoronoiFVM.NewtonControl()
	control.Δt_min=0.01*tstep
	control.Δt=tstep
	control.Δt_max=0.1*tend
	control.Δu_opt=0.1
	control.Δt_grow=dtgrowth
	
	tsol=solve(inival,sys,[t0,tend];control=control)
	return grid,tsol
end;

# ╔═╡ 7c4abd9e-f754-451a-8a29-aa9c3e8c712d
gridd2d,tsold2d=Dirichlet2D(dim=2,n=N);

# ╔═╡ 68dd0792-1884-4fb1-bc81-e860fb610701
md"""
Time step number:$(@bind td2d Slider(1:length(tsold2d),default=10,show_value=true))
"""

# ╔═╡ 7f5303ca-3a9f-48af-88f6-5907f2afc5e9
begin
	
	pd2=GridVisualizer(Plotter=PyPlot,layout=(3,1),fast=true,resolution=(600,600))
    	
		
scalarplot!(pd2[1,1],gridd2d,tsold2d[1,:,td2d],title=@sprintf("Numerical Solution with Dirichlet Boundary Condition, t=%.3g",tsold2d.t[td2d]))
		
scalarplot!(pd2[2,1],gridd2d,map((x,y)->barenblatt2d(x,y,tsold2d.t[t_b2d],m),gridd2d),title=@sprintf("Exact Barenblatt Solution, t=%.3g",tsold2d.t[t_b2d]))
		
		scalarplot!(pd2[3,1], gridd2d,map((x,y)->barenblatt2d(x,y,tsold2d.t[td2d],m),gridd2d)-tsold2d[1,:,t_b2d],title=@sprintf("Error of the Numerical Solution in 2D, t=%.3g",tsold2d.t[td2d]))
	
        reveal(pd2)
    	
	gcf()
end

# ╔═╡ 3d60a1e7-c56c-495b-b4c4-5d401db9e3e7
md"
# References
1.[Scientific Computing course winter semester 20/21, Prof. Fruhmann](https://www.wias-berlin.de/people/fuhrmann/SciComp-WS2021/)

2.[VoronoiFVM.jl](https://j-fu.github.io/VoronoiFVM.jl/stable/)

3.[ExtendableGrids.jl](https://j-fu.github.io/ExtendableGrids.jl/dev/)

4.[The Porus Medium Equation, Juan Luis Vazquez](https://oxford.universitypressscholarship.com/view/10.1093/acprof:oso/9780198569039.001.0001/acprof-9780198569039#:~:text=fairly%20well%20understood.-,This%20book%20provides%20a%20presentation%20of%20the%20mathematical%20theory%20of,%2C%20heat%20transfer%2C%20or%20diffusion.)

5.[Iterative methods for sparse linear systems, Y.Saad](http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf)

6.[The Forward Euler Method] (https://www.algorithm-archive.org/contents/forward_euler_method/forward_euler_method.html)

7.[Gas Transport in Porus Media, Clifford K. Ho and Stephen W. Webb](https://www.springer.com/gp/book/9781402039614)

8.[Porous medium equation, A. Wathen and L. Qia](https://people.maths.ox.ac.uk/trefethen/pdectb/porous2.pdf)

9.[]()

"

# ╔═╡ Cell order:
# ╟─584c650e-a819-4c98-809d-1101ec1e83f8
# ╟─1909662f-efa8-4f60-a7d3-33ccbff69897
# ╟─98deb536-382e-42d9-8809-d2bcda17420a
# ╟─a382d48f-835e-470f-8cf4-a8ac35a9d45e
# ╟─f13a6d3b-546b-48c2-8190-97644bc17837
# ╟─ed6f83f7-3db4-4a0e-b867-0cace9aa444e
# ╟─0c59788a-6e6b-4a13-b305-3ceee53fd8f0
# ╟─fc300c0e-12db-4ce3-a26e-0b27591f23a1
# ╟─744983d7-eb8d-4cfe-8d71-06f631e4b13e
# ╟─04cf9a64-3c8d-426e-bae3-93867508d978
# ╟─13f173f3-15ea-4110-861a-4aa663a506ca
# ╟─e1777730-babd-42f6-8022-36ff9e6e7d6e
# ╟─af2facd0-ba67-4821-b83b-2f1e09157770
# ╟─408a53a3-beaa-497c-aa18-c1c2601eda4f
# ╟─e3e56550-825e-4905-bcb8-3c8086671c28
# ╟─d2c895f0-4b44-4260-a2c2-49b06b43f063
# ╟─9c9dcd8a-e79a-45a4-be26-7a5f40348f16
# ╠═1ce4a24c-2be3-4b9b-9045-72b90f329df1
# ╠═14730986-22e5-4c63-a076-2a84686f3edf
# ╠═7f4b9a30-ad19-4e21-8bf1-a443f503ca05
# ╠═9a384177-ab71-4b9e-9a04-263fddb68718
# ╠═7ee60df1-a27f-476b-a40e-9cd1a942d0ba
# ╠═d0f8c2ce-52dd-4451-986a-cea9fd3bd35e
# ╟─eb3b35da-199f-4e5a-b294-a3f428f53535
# ╟─6b373f2e-61e9-4219-aeb4-6f86cc7a2189
# ╟─ae1adee2-cb6f-471c-a9b7-8c67ffd2f844
# ╟─6207ed96-c2d6-4fc3-8a71-29f6f0101b61
# ╟─190d69d8-8980-4ccc-b978-73256cd56de5
# ╟─3d56289e-c1ad-4751-be54-14e35f88dd5a
# ╟─5b923416-0c28-4a10-8d0b-b4cd62e5ffba
# ╟─6e1eb8c9-9e6b-4f5a-b01b-0cca7b10f80e
# ╟─2408a0e9-532c-4be5-ab38-61fb5437aba1
# ╟─8bc22a4f-0def-4d51-b9bf-310a40ea8482
# ╟─ba9f286c-b6ab-4986-976b-761b6c53d26a
# ╟─0f2f3926-34ad-4660-8c13-35dee4ec6474
# ╟─bf1a6b3a-f87f-470e-9074-47f475d760f8
# ╟─998fb1e8-ff11-46a7-bb41-98bbb3eb76f5
# ╟─63a41d85-79ad-49b1-aefd-c7a4fd250cea
# ╟─e6f3a74c-9161-4b83-bd5d-61028b4e81dd
# ╟─041a1492-0aaf-4187-b7ed-cd49b092b727
# ╟─988f0e80-bf87-4c78-bba2-a04d53db519c
# ╟─d9300a5a-0682-4f72-862f-ca90215b553f
# ╟─36f4fd5c-d093-4eb9-8956-4dc4af0d34f7
# ╟─98f441ac-f7fa-4979-b0b6-c4ced36f84b6
# ╟─8ea3bb12-d0f2-487d-a8ad-3d5030a4ffa7
# ╟─46b92be8-c262-4c80-83c7-ae3b12d01c71
# ╟─0a00fa12-0e55-44f3-9685-e5d1cd7adb52
# ╟─2dd023d3-f1b6-404f-b3f0-d8997f121ee2
# ╟─b037aa3b-c4d3-44c5-9b75-a075657fdd81
# ╟─034beddc-3a6b-4ea5-bcf7-73718adf80b8
# ╟─a62b283f-c407-4c8f-911d-486a99b07bec
# ╟─978aca9c-7e19-485d-a943-45e72d4b83a6
# ╟─24a8924f-b9ef-4cce-bf11-d0803f254b95
# ╟─0de74633-25d1-4c06-a38d-f2deb982cd2a
# ╟─7c1eb0a7-ee28-4fb5-b560-8040793202c8
# ╟─86da8092-f771-4ef0-a6c0-4bcb45ec4173
# ╟─276ebe5b-7799-4a41-901c-931972a65e95
# ╟─084cb624-fdeb-42e6-9c1b-c5f03beaa32f
# ╟─a03c19d7-2b7b-456b-8dbe-db2809221f2d
# ╟─944b1ade-abd3-49f8-91d9-b01e8254911e
# ╟─39afbd01-44bd-47b1-9ab1-aadd019bd08d
# ╟─e014b546-a855-486b-a5f9-cdda73611bc5
# ╟─8931a3fb-43fd-4a38-b760-6e4d8ca6ad77
# ╟─bf505926-67b8-4071-8a86-f7dd87850bfc
# ╟─5ab09487-01bf-437d-8e45-53f97fbc9ff0
# ╟─658cb684-7300-4f75-ac86-db32a753aaa0
# ╟─7c4abd9e-f754-451a-8a29-aa9c3e8c712d
# ╟─68dd0792-1884-4fb1-bc81-e860fb610701
# ╟─7f5303ca-3a9f-48af-88f6-5907f2afc5e9
# ╠═f7fd554e-3885-4c58-afa8-26b00dba8b96
# ╠═832fd72b-0d51-44f5-97f0-7ad261d65e13
# ╠═35a21196-41c4-4ed7-8eae-8d559bd5500e
# ╠═84ad06f5-72c9-456a-aefe-f4fff86ce60d
# ╟─16dec684-3cbb-45e6-8761-16734643a2cb
# ╟─1a3f82ef-9ac3-48ae-abea-9b150b59a955
# ╟─c46db835-f4f0-4c09-9e16-b093d73c8312
# ╠═524d5cef-debd-47a3-8ff5-2f25341d57b3
# ╠═c14c98b0-0fb0-4572-a721-a9cffba63e3e
# ╠═1e0c8cbb-58b2-4019-9e30-1589f0ebf92f
# ╟─3d60a1e7-c56c-495b-b4c4-5d401db9e3e7
