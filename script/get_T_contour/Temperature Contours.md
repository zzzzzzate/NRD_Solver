<br><span class='markdown-page-line'>---------------------------------------------<span id='page1' class='markdown-page-text'>[ 第1页 ]</span>---------------------------------------------</span><br><br>

# Temperature Contours and Ghost Surfaces for Chaotic Magnetic Fields  

S. R. Hudson and J. Breslau
Princeton Plasma Physics Laboratory, PO Box 451, Princeton, New Jersey 08543, USA
(Received 5 November 2007; published 3 March 2008)  

Steady state solutions for anisotropic heat transport in a chaotic magnetic field are determined numerically and compared to a set of "ghost surfaces"—surfaces constructed via an action-gradient flow between the minimax and minimizing periodic orbits. The ghost surfaces are in remarkable agreement with the temperature contours.  

DOI: 10.1103/PhysRevLett.100.095001  

PACS numbers: 52.25.Fi, 05.45.Pq  

A variety of transport processes in magnetically confined plasmas are dominated by strong parallel transport along the magnetic field $\mathbf{B}$, with small perpendicular transport. Coordinates adapted to the structure of the magnetic field, magnetic coordinates, therefore provide an elegant theoretical description of plasma dynamics and often enhance numerical accuracy. Magnetic coordinates are analogous to the action-angle coordinates of Hamiltonian systems and may be constructed globally when the magnetic field lines lie on nested, invariant toroidal surfaces, i.e., when the field is integrable. Integrable magnetic fields are, however, the exception rather than the rule. Error fields [1] or internal plasma motions, e.g., microtearing instabilities [2], result in partially chaotic magnetic fields in tokamaks, and chaotic fields are intrinsic to the nonsymmetric stellarator [3].  

Here we present a coordinate framework adapted to the structure of chaotic magnetic fields, which we call chaotic magnetic coordinates, and show that this framework allows a simple description of anisotropic transport. We consider heat transport, as described by  

$$
\begin{array}{r} { \frac { \partial T } { \partial t } = \nabla \cdot ( \kappa _ { \| } \nabla _ { \| } T + \kappa _ { \perp } \nabla _ { \perp } T ) + Q . } \end{array}
$$  

where $T$ is the temperature, $t$ is time, and $\kappa_{\parallel}$, $\kappa_{\perp}$ are the (constant) parallel and perpendicular diffusion coefficients. The parallel derivative, $\nabla_{\parallel}T$, is given $\nabla_{\parallel}T = \mathbf{b}\mathbf{b} \cdot \nabla T$, where $\mathbf{b} = \mathbf{B}/|B|$, and the perpendicular derivative is $\nabla_{\perp}T = \nabla T - \nabla_{\parallel}T$. The term $Q$ allows for heat sources/sinks, but we set this to zero and examine the nontrivial, steady state solutions forced by inhomogeneous boundary conditions.  

For fusion plasmas, the ratio $\kappa_{\parallel}/\kappa_{\perp}$ may exceed $10^{10}$ [4]. Strong anisotropy has different consequences, depending on whether the magnetic field lines lie on nested flux surfaces, whether the field is slightly chaotic, or whether the field is so chaotic that the motion of field lines is effectively random. In the first case the temperature is a surface function, $T=T(\psi)$, where $\psi$ labels flux surfaces, and gradients can be supported. For the opposite case of extreme chaos, where the field lines seem to wander randomly over a volume, the strong parallel transport results in temperature flattening, $T=\mathrm{const}$. It is the intermediate  

case of critical (near-threshold) chaos that is most relevant for toroidal plasma confinement. The temperature is then dominated by the fractal structure of the chaotic magnetic field. How chaotic magnetic coordinates allow this structure to be understood is the topic of this Letter.  

A chaotic magnetic field is a fractal mix of (i) invariant flux (KAM) surfaces [5,6], which are labeled by their irrational rotational transform, (ii) cantori (broken KAM surfaces), in particular, the near-critical cantori which present effective but partial barriers to field-line transport [7], (iii) unstable periodic orbits and their unstable manifolds which constitute the stochastic sea, and (iv) stable periodic orbits and elliptic island chains [5,6].  

The complexity of the field structure dictates that Eq. (1) must be solved numerically [8,9], but this is not an easy task. The temperature must be represented as a scalar field of three-dimensional space, $T = T(\psi, \theta, \phi)$, where $\theta$, $\phi$ are arbitrary poloidal and toroidal angles. The infinitely many irregular field lines in the stochastic sea may come arbitrarily close to each other. For large $\kappa_{\parallel}$ the temperature along the field lines is almost constant, and for small $\kappa_{\perp}$ the cross field interaction is very weak. The temperature becomes a fractal function of position as $\kappa_{\parallel}/\kappa_{\perp}$ increases and the resolution requirements become overwhelming. The challenge is to achieve sufficient accuracy to resolve the near-fractal structure, ensuring that numerical error, "numerical diffusion", does not overwhelm the small perpendicular diffusion.  

It would be of great benefit if some theoretical insight allowed the representation of the temperature to be simplified. For example, on the KAM surfaces, we may expect that the temperature will be constant. We also know [4] that the temperature will flatten inside the island chains when the island width, $\Delta w$, exceeds a critical value, $\Delta w \sim (\kappa_{\perp}/\kappa_{\parallel})^{1/4}$. Within the stochastic sea, it is tempting to conclude that the strong parallel transport results in a flat temperature profile, or that the transport is uniform. For near-threshold chaos, however, this is an oversimplification. Irregular trajectories, with finite Lyapunov exponent, may take an impractically long time to sample the accessible volume. Attempts to determine transport by averaging [10] must take into account that within the stochastic  

<br><span class='markdown-page-line'>---------------------------------------------<span id='page2' class='markdown-page-text'>[ 第2页 ]</span>---------------------------------------------</span><br><br>

sea there exists a finite volume of regular motion (the magnetic islands), and what the relative volume of irregular versus regular motion is remains an open question in nonlinear dynamics. The point is, chaos is not random.  

The key to understanding the structure of the temperature in the stochastic sea is to realize that the most effective barriers to field-line transport are given by the cantori. Cantori are the invariant sets under the field-line flow remaining after a KAM surface has been destroyed by chaos [11–13], but they have an infinity of gaps where field-lines may leak through. In the near-critical case (when the level of chaos just exceeds that required to break the KAM surface) the gaps in the cantorus are small, and the field-line flux across the cantorus is small. As the most robust KAM surfaces have noble rotational transform [14], the most important barriers to field-line transport in chaotic fields are usually the noble cantori. As the level of chaos increases, the gaps in the cantorus enlarge and the field-line flux increases: supercritical cantori have little effect on field-line transport.  

So we have a situation in which regions of local temperature flattening are produced by the significant islands, between which the irrational barriers may support gradients. If coordinate surfaces can be constructed that coincide with the irrational barriers, then the temperature profile will approximate a smoothed devil's staircase [15]. Clearly, coordinate surfaces should coincide with any KAM surfaces that exist, but here we consider a region in which all KAM surfaces are destroyed and the most significant barriers are provided by the noble cantori. To construct a coordinate framework based on cantori we need to "close the gaps", and this can be done by constructing ghost surfaces, as we now describe.  

Cantori are approximated by high-order, action-minimizing periodic orbits [16]. These are conveniently found using the action formalism of magnetic-field-line dynamics [17]. The action formalism is also required for the construction of the ghost surfaces, which are defined using the action-gradient flow. Magnetic field lines are stationary curves $\mathcal{C}$ of the action integral [18],  

$$
\begin{array}{r} { s _ { c } = \int _ { c } \mathbf { A } \cdot d l , } \end{array}
$$  

where $\mathbf{B} = \nabla \times \mathbf{A}$. We use a vector potential in canonical form $\mathbf{A} = \psi \nabla \theta - \chi \nabla \phi$, where $\chi(\psi, \theta, \phi)$ is the field-line Hamiltonian:  

$$
\chi = \psi ^{2} / 2 + \sum \chi _ { m , n } \cos ( m \theta - n \phi ) .
$$  

A piecewise-linear approximation for $\mathcal{C}$ is sufficient, where between $\phi \in [i\Delta\phi, (i+1)\Delta\phi]$ a "trial" curve is given $\theta(\phi) = \theta_i + (\theta_{i+1} - \theta_j)(\phi - \phi_i)/\Delta\phi$ for $\Delta\phi = 2\pi q/N$, and $\psi = \hat{\theta}(\phi)$ [17]. We restrict attention to $(p, q)$ periodic curves, $\theta(\phi + 2\pi q) = \theta(\phi) + 2\pi p$, by constraining $\theta_N = \theta_0 + 2\pi p$. The action integral is now piece  

wise directly solvable and is a rapidly computable function of the $N$ independent parameters, $S(\theta_0, \theta_1, \ldots, \theta_{N-1})$. Periodic orbits are those particular trial curves for which the action gradient, $\nabla S = (\partial S/\partial\theta_1, \partial S/\partial\theta_2, \ldots)^T$, is zero. Finding periodic orbits amounts to a multidimensional root find, and an $N$-dimensional Newton method is suitable. The derivative of the action gradient, the Hessian $D^2S$, is a cyclic, tridiagonal matrix of the second derivatives of $S$. The action extremizing approach allows both the stable (minimax) and unstable (minimizing) orbits to be quickly found, even for orbits with periodicities in the tens of thousands for strongly chaotic fields [17].  

The Hessian at the minimax orbit generically has a single negative eigenvalue, and the associated eigenvector indicates the direction in configuration space along which the action integral decreases. Ghost surfaces are constructed by pushing a trial curve off the minimax orbit in this direction, then allowing the curve to evolve down the gradient flow:  

$$
\frac { d \theta _ { i } } { d \tau } = - \frac { \partial S } { \partial \theta _ { i } } ,
$$  

where $\tau$ is any suitable integration parameter. As the action is decreasing under this flow, and the curves are constrained to be periodic, the trial curve will evolve into the minimizing periodic orbit, and in doing so will trace out a surface, the ghost surface of periodicity $(p, q)$.  

Ghost surfaces were originally introduced for the standard map [19,20] (in this context, they are called ghost circles), and they were found to be nonintersecting: we have not found exceptions to this. Any selection of ghost surfaces may form the framework of the chaotic coordinates, and by choosing rationals $p/q$ that approximate a given irrational we may consider irrational ghost surfaces. To complete the chaotic coordinates, the surfaces can be interpolated radially to provide a continuous foliation of space, and a suitable angle coordinate can be imposed, for example, so that each trial curve comprising the ghost surface is straight.  

Intuition suggests that the irrational ghost surfaces that close the gaps in near-critical cantori would coincide with temperature isocontours. What was unexpected is how closely the ghost surfaces coincide for even the strongly supercritical cantori.  

To compute steady state solutions of Eq. (1), a second-order finite-difference model is employed. The parallel and perpendicular diffusions are separated numerically [21] by locally introducing straight-field-line (Clebsch) coordinates ($\alpha$, $\beta$, $\phi$), where $\mathbf{B} = \nabla\alpha \times \nabla\beta$. The parallel-diffusion operator becomes  

$$
\nabla _ { \parallel } ^{2} T = B ^{\phi} \frac { \partial } { \partial \phi } \left( \frac { B ^{\phi} } { B ^{2} } \frac { \partial T } { \partial \phi } \right) ,
$$  

where the partial derivative with respect to $\phi$ is along a magnetic-field line: for each grid point ($\psi_{i,j}$, $\theta_{i,j}$) on the  

<br><span class='markdown-page-line'>---------------------------------------------<span id='page3' class='markdown-page-text'>[ 第3页 ]</span>---------------------------------------------</span><br><br>

plane $\phi_{k}=k\Delta\phi$, with temperature $T_{i,j,k}$, the parallel gradient on the forward "half-$\phi$" grid is approximated  

$$
\left. \frac { \partial T } { \partial \phi } \right| _ { i , j , k + 1 / 2 } = \frac { T ( \psi , \theta , \phi _ { k + 1 } ) - T _ { i , j , k } } { \Delta \phi } ,
$$  

where ($\psi$, $\theta$, $\phi_{k+1}$) is where the field line starting from ($\psi_{i,j}$, $\theta_{i,j}$, $\phi_{k}$) intersects the $\phi_{k+1}$ plane, which can always be determined by field-line tracing. In general, this point will not coincide with a grid point, so bilinear interpolation is used to estimate $T(\psi, \theta, \phi_{k+1})$. The quantity $\partial_{\phi} T|_{i,j,k-1/2}$ on the backward half-$\phi$ grid is defined similarly. The first partial $\phi$-derivatives on the $k+\frac{1}{2}$ and $k-\frac{1}{2}$ half-grids are combined, along with the factors $B^{\phi}$ and $B^{2}$, to give a centered, finite-difference realization of the second-order, parallel-diffusion operator.  

For $\kappa_{\parallel} \gg \kappa_{\perp}$, the temperature will vary weakly along magnetic field lines. So, $\Delta\phi$ need not be small. We choose $\Delta\phi = 2\pi$ and perform the computation on a single plane. This reduces the computational burden and allows additional resolution within the plane, i.e., in the perpendicular direction, which is required to resolve the small scale of the solution for small $\kappa_{\perp}$.  

The diffusion perpendicular to $\mathbf{B}$ is approximated by a diffusion instead perpendicular to $\phi$. This approximation introduces a negligible error when the field is dominantly toroidal and when $\kappa_{\perp}/\kappa_{\parallel}$ is small, and it eliminates the need to compute the metric elements of the ($\alpha$, $\beta$, $\phi$) coordinates, which in principle are determined from differentiating the field-line integration (i.e., constructing the tangent map). The diffusive operator perpendicular to $\phi$ is given by the Laplacian  

$$
\nabla _ { \perp } ^{2} T = \sqrt { g } ^{- 1} [ \partial _ { \psi } ( \sqrt { g } T ^{\psi} ) + \partial _ { \theta } ( \sqrt { g } T ^{\theta} ) ] ,
$$  

where $T^{\psi} = g^{\psi\theta}T_{\psi} + g^{\psi\theta}T_{\theta}$ and $T^{\theta} = g^{\theta\psi}T_{\psi} + g^{\theta\theta}T_{\theta}$, where $T_{\psi} = \partial T/\partial\psi$ and $T_{\theta} = \partial T/\partial\theta$, and the geometric information is encapsulated in the "raising" metric elements $g^{ab} = \nabla a \cdot \nabla b$ and the Jacobian, $\sqrt{g}$. The Laplacian is discretized using second-order finite differences [22].  

The steady state condition,  

$$
\kappa _ { \parallel } \nabla _ { \parallel } ^{2} T + \kappa _ { \perp } \nabla _ { \perp } ^{2} T = 0 ,
$$  

becomes a sparse linear system which is solved using an iterative Krylov method [biconjugate gradient stabilized method (Bi-CGStab) [23]]. We consider the region between two magnetic islands, namely, the $(p,q)=(1,2)$, $(2,3)$ islands at $\psi=\frac{1}{2}$ and $\psi=\frac{2}{3}$, respectively, which are excited by the $\chi_{2,1}$ and $\chi_{3,2}$ perturbation harmonics in Eq. (3), and we set $2\chi_{2,1}=3\chi_{3,2}=k$, where $k$ is a perturbation parameter. The symmetry of the field allows $T(\psi,-\theta)=T(\psi,\theta)$, so a regular grid in $\psi$, $\theta$ is constructed in the region $\psi\in[\psi_l,\psi_u]$ and $\theta\in[0,\pi]$, where $\psi_l=0.50$ and $\psi_u=0.68$, with grid spacing $\Delta\psi=(\psi_u-$  

$\psi_{l})/N$, $\Delta\theta=2\pi/N$, where $N$ is the grid resolution. It is the chaotic structure of the field that is relevant to the present study, rather than geometry, so we use the simple Cartesian metric, $g^{\psi\psi}=g^{\theta\theta}=1$, $g^{\psi\theta}=0$, and $\sqrt{g}=1$. The most robust KAM surface in this region appears to be the $\tau=0.5607\ldots$ surface, which has a critical perturbation $k=2.039\times10^{-3}$ [17], so here we set the perturbation $k=2.100\times10^{-3}$ to just exceed this critical value to give a field with connected chaos between the $(1,2)$ and $(2,3)$ islands. A Poincaré plot of this field is shown in Fig. 1. The boundary conditions are $T(\psi,\theta)=1$ for $\psi\leq\psi_{l}$, and $T(\psi,\theta)=0$ for $\psi\geq\psi_{u}$. We have confirmed the second-order scaling of the error with respect to grid size, $e\sim\mathcal{O}(N^{-2})$, and the expected scaling of the critical island width $\Delta w\sim(\kappa_{\perp}/\kappa_{\parallel})^{1/4}$. Temperature isocontours are shown in Fig. 1 for the case $\kappa_{\perp}/\kappa_{\parallel}=10^{-10}$, with $N=2^{12}$.  

There is a countable infinity of ghost surfaces that may be selected: the optimal selection is determined by the island widths and $\kappa_{\perp}/\kappa_{\parallel}$. (An island width is not well defined when the separatrix becomes chaotic, but one could instead consider the resonance area [24].) We distinguish three types of surface: (i) low-order surfaces, (ii) high-order surfaces where $p/q$ approximates a noble irrational, and (iii) high-order surfaces where $p/q$ approximates a boundary irrational (an irrational that lies close to a low-order rational [25]). When $\kappa_{\perp}$ is comparatively large, the fine scale structure of the field is overlooked and the low-order islands have the dominant effect on the  

<div style="display: block; width: 100%"><img src="https://storage.simpletex.cn/view/m6687e9e2d2455445eafcd09f3f619a45" style="width: 38%; max-width: 38%" /></div>  
FIG. 1 (color online). For $\theta < 0$: the selected ghost curves (red lines) and cantori (black dots). For $\theta > 0$: Poincaré plot (gray dots), ghost curves (red lines), and the temperature contours (black lines) for $\kappa_{\perp}/\kappa_{\parallel} = 10^{-10}$.  

<br><span class='markdown-page-line'>---------------------------------------------<span id='page4' class='markdown-page-text'>[ 第4页 ]</span>---------------------------------------------</span><br><br>

<div style="display: block; width: 100%"><img src="https://storage.simpletex.cn/view/m547e8efc01678caa7d3e07a9eebe579b" style="width: 39%; max-width: 39%" /></div>  
FIG!  

solution!  

For!  

Given!  

This!  

[!