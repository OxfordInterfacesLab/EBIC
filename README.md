# EBIC Electron Induced Beam Current modelling scripts.
For full scientific description see Zhou, et al, Journal of Applied Physics, 2020.

## Generation of Carriers using Werner's function
Werner’s description is used here:
The generation function is calculated as:
(2)	G(x,y=0,z)=δε/δxδyδz×N_(e-h) 	,
with the depth dose distribution described by δε/δxδyδz and N_(e-h) being the number of photogenerated electron-hole pairs per second. The volumetric generation assumes cylindrical symmetry:

(3)	ϕ(x,y,x)=δε/δxδyδz=∑_(i=1)^2▒〖c_i/(2πσ_ui σ_zi )×exp⁡(-〖(x-x_0)〗^2/(2σ_xi^2 )-〖(y)〗^2/(2σ_xi^2 )-(z-z_i )^2/(2σ_zi^2 )) 〗 ,

with x,y,z representing the coordinates system pictured in Figure 1, and x_0 being the location of the beam in the x axis. The beam was assumed to be centered on the y axis, y=0. 

N_(e-h) is determined by equating the total generation of carriers from the volume integral of the dose as: 
(4)	G_total=N_(e-h) ∫_V▒ϕ(x,y,x) =I_b/q  1/E_ehp  E_b (1-η_bsc  E_b/E_b0 ) ,

where I_b is the beam current in Amps,  E_ehp is the average energy required to create one electron-hole pair in Si, here set to 3.8 eV to use a mean value from those reported in the literature (3.5-4 eV). E_b is the beam energy in eV, such that the product E_b (1-η_bsc E_b/E_b0  ) represents the fraction of electrons which are absorbed rather than backscattered. Such a correction factor for the energy loss due to backscattered electrons (1-η_bsc E_b/E_b0  ) was set to 0.9 after the work of Wu and Wittry.

## Collection of carriers accross a PN junction device

The EBIC signal is given by the solution to the continuity equations for electron and hole dynamics:
(S1)	∇(σ_n ∇ϕ_Fn )=q(G-R)
(S2)	 ∇(σ_p ∇ϕ_Fp )=-q(G-R)
Where q is the fundamental electron charge, ϕ_Fn, ϕ_Fp are the quasi Fermi energies for electrons and holes referred to the equilibrium value ϕ_F0, G is the carrier generation rate, R is the carrier recombination rate given by R=(n-n_0)/τ_eff  =(p-p_0)/τ_eff  ,  with τ_eff denoting the effective carrier lifetime. σ is the carrier conductivity given by:
(S3)	σ_n=qnμ_n
(S4)	σ_p=qpμ_p
Here n,p are the electron and hole concentrations, and μ represents the carrier mobility which was assumed constant as 1110 cm2V-1s-1 and 430 cm2V-1s-1 for electron and holes, respectively. To reduce computational time, only Boltzmann statistics were considered, and thus the carrier concentrations were defined as:
(S5)	n=n_0+δn=n_i  exp⁡((ϕ_Fn-ϕ_i+ψ)/V_t )	
(S6)	p=p_0+δn=n_i  exp⁡((-ϕ_Fp+ϕ_i-ψ)/V_t )
Where n_i is the intrinsic carrier concentration, ψ the electrostatic potential, δn,δp are the excess carrier concentrations,  V_t is the thermal voltage given by kT/q≈25 mV at 300 K, and n_0,p_0 are the electron and hole equilibrium carrier concentrations given by:
(S7)	n_0=n_i  exp⁡((ϕ_F0-ϕ_i)/V_t )
(S8)	p_0=n_i  exp⁡(〖ϕ_i-ϕ〗_F0/V_t )
ϕ_F is the equilibrium Fermi energy. Care was taken to not exceed carrier densities of 8x1017 cm-3, such that Boltzmann statistics were still applicable. Establishing the value of ψ, and the reference ϕ_F0-ϕ_i is achieved via Poisson’s equation:
(S9)	∇(∇ψ)=-q/(ϵ_0 ϵ_Si ) (p-n+N_d^+-N_a^- ).
In Equation (S9), N_d^+ represents the concentration of ionized donors, while N_a^- the concentration of ionized acceptors. Here, complete ionization is also assumed. The ϵ_0 ϵ_Si product is the silicon permittivity. The electrostatic potential ψ is determined either by the development of the depletion region in the pn-junction, by surface band bending, or by the application of a voltage at any given interface. It is useful to note the mass law action follows from the equations above:
(S10)	n_i^2=n_o p_0, and  np=n_i^2  exp⁡((ϕ_Fn-ϕ_Fp)/V_t )
The quasi-neutrality condition p ≈ n + N_a  – N_d was assumed at the metal collecting boundaries (emitter and base), producing a solution for equilibrium carrier concentration as:
(S11)	n_0=-(N_A-n_0)/2+√(((N_A-n_0)/2)^2+n_i^2 )
Such that the equilibrium Fermi energy at the contacts could be determined by:
(S12)	ϕ_F0-ϕ_i=V_t  log⁡( 1/n_i  (-(N_A-N_d)/2+√(((N_A-N_d)/2)^2+n_i^2 )))
The solution to the EBIC model thus reduces to numerically finding ϕ_Fn, ϕ_Fp, and ψ.
The following boundary conditions were applied to find such solution. 
At the cross section surface where the electron beam impacts (Figure 1 in the manuscipt), a non-conductive boundary condition was applied. Carrier flux is determined by a current J_rec, and the total current into the boundary is zero σ_n ∇ϕ_Fn+σ_p ∇ϕ_Fp=0, thus:
(S13)	σ_n ∇ϕ_Fn=-J_rec, and  σ_p ∇ϕ_Fp=J_rec
The recombination current is given by:
(S14)	J_rec=qU_s=q (n_s p_s-n_i^2)/((n_s+n_i)/S_p0+(p_s+n_i)/S_n0 )
Where U_s represents the surface recombination rate, S_n0,S_p0 represents the capture velocities for electrons and holes, and serves as a measure of how much surface recombination takes place at a mid-gap defect [1]. The electrostatic potential was set as a neutral surface ∇ψ=0, or a fixed amount of band bending was introduced by setting ψ=〖ψ_0+ψ〗_sbb at such surface.
At the collecting boundaries, n- and p-type sides, the electrostatic potential was set to the equilibrium Fermi energy difference ψ_0=ϕ_F0-ϕ_i, while the quasi Fermi energies were set to zero to force perfect collection of carriers at the contacts in a short-circuit condition:
(S15)	ϕ_Fn=0, and ϕ_Fp=0
Lastly, at symmetry boundary the material is assumed semi-infinite and no current should flow:
(S16)	σ_n ∇ϕ_Fn=0, σ_p ∇ϕ_Fp=0, and ∇ψ=0
The solution to ϕ_Fn, ϕ_Fp, and ψ was found using Matlab’s Partial Differential Equation Toolbox. The EBIC signal was then found by calculating the total current flowing on the n- and p-type contacts as:
(S17)	J_emitter=1/width ∫_0^width▒dz (σ_n ∇ϕ_Fn+σ_p ∇ϕ_Fp )  |_(x=0)
(S18)	J_base=1/width ∫_0^width▒dz (σ_n ∇ϕ_Fn+σ_p ∇ϕ_Fp )    |_(x=length)
The discrepancy between the current densities in the emitter and the base contacts was used as a measure for model accuracy. When meshing or cell dimensions are chosen to be too coarse these values will differ substantially, thus indicating neutrality was not fully achieved.
