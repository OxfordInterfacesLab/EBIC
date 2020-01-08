# EBIC Electron  Beam Induced Current 

**Modelling Scripts**

For full scientific description see Zhou, et al, Journal of Applied Physics, 2020.

To run simulation of a beam scan accross a pn junction run:

Ebic2D_Func(0,1e6\*1e4,5);   0 band bending , SRV= 1e6 cm/s \*1e4 , Ebeam =5 kV

The location of the beam on the x-axis can be specified in variable *x_beam* in Ebic2D_beamscan.m.

## Generation of Carriers using Werner's function
Werner’s description is used here: U. Werner, F. Koch, and G. Oelgart, “Kilovolt electron energy loss distribution
in Si,” J. Phys. D Appl. Phys. 21, 116–124 (1988).

To test this function un Werner_eBeamCarrierGen.m

The generation function is calculated as:

![image](https://user-images.githubusercontent.com/53188769/71976768-0c457c00-320f-11ea-8680-bb906833b2b3.png)

with the depth dose distribution described by δε/δxδyδz and N_(e-h) being the number of photogenerated electron-hole pairs per second. The volumetric generation assumes cylindrical symmetry:

![image](https://user-images.githubusercontent.com/53188769/71976780-11a2c680-320f-11ea-918e-faba8fdb4556.png)

with x,y,z representing the coordinates system pictured in Figure 1, and x_0 being the location of the beam in the x axis. The beam was assumed to be centered on the y axis, y=0. 

N_(e-h) is determined by equating the total generation of carriers from the volume integral of the dose as: 

![image](https://user-images.githubusercontent.com/53188769/71976792-1798a780-320f-11ea-99f2-eae6daaa2f85.png)

where I_b is the beam current in Amps,  E_ehp is the average energy required to create one electron-hole pair in Si, here set to 3.8 eV to use a mean value from those reported in the literature (3.5-4 eV). E_b is the beam energy in eV, such that the product E_b (1-η_bsc E_b/E_b0  ) represents the fraction of electrons which are absorbed rather than backscattered. Such a correction factor for the energy loss due to backscattered electrons (1-η_bsc E_b/E_b0  ) was set to 0.9 after the work of Wu and Wittry.

## Collection of carriers accross a PN junction device

The EBIC signal is given by the solution to the continuity equations for electron and hole dynamics:

![image](https://user-images.githubusercontent.com/53188769/71976459-482c1180-320e-11ea-81a0-0756d2a5af07.png)


Where q is the fundamental electron charge, ϕ_Fn, ϕ_Fp are the quasi Fermi energies for electrons and holes referred to the equilibrium value ϕ_F0, G is the carrier generation rate, R is the carrier recombination rate given by R=(n-n_0)/τ_eff  =(p-p_0)/τ_eff  ,  with τ_eff denoting the effective carrier lifetime. σ is the carrier conductivity given by:

![image](https://user-images.githubusercontent.com/53188769/71976471-4f531f80-320e-11ea-997e-684602da4d8c.png)

Here n,p are the electron and hole concentrations, and μ represents the carrier mobility which was assumed constant as 1110 cm2V-1s-1 and 430 cm2V-1s-1 for electron and holes, respectively. To reduce computational time, only Boltzmann statistics were considered, and thus the carrier concentrations were defined as:

![image](https://user-images.githubusercontent.com/53188769/71976480-55490080-320e-11ea-9f95-b67d7b98760d.png)

Where n_i is the intrinsic carrier concentration, ψ the electrostatic potential, δn,δp are the excess carrier concentrations,  V_t is the thermal voltage given by kT/q≈25 mV at 300 K, and n_0,p_0 are the electron and hole equilibrium carrier concentrations given by:

![image](https://user-images.githubusercontent.com/53188769/71976493-5f6aff00-320e-11ea-88d3-ef3684ab80e2.png)

ϕ_F is the equilibrium Fermi energy. Care was taken to not exceed carrier densities of 8x1017 cm-3, such that Boltzmann statistics were still applicable. Establishing the value of ψ, and the reference ϕ_F0-ϕ_i is achieved via Poisson’s equation:

![image](https://user-images.githubusercontent.com/53188769/71976515-672aa380-320e-11ea-8bb4-aca7179b36a9.png)

In Equation (S9), N_d^+ represents the concentration of ionized donors, while N_a^- the concentration of ionized acceptors. Here, complete ionization is also assumed. The ϵ_0 ϵ_Si product is the silicon permittivity. The electrostatic potential ψ is determined either by the development of the depletion region in the pn-junction, by surface band bending, or by the application of a voltage at any given interface. It is useful to note the mass law action follows from the equations above:

![image](https://user-images.githubusercontent.com/53188769/71976523-6abe2a80-320e-11ea-9a88-f1667dd20a37.png)

The quasi-neutrality condition p ≈ n + N_a  – N_d was assumed at the metal collecting boundaries (emitter and base), producing a solution for equilibrium carrier concentration as:

![image](https://user-images.githubusercontent.com/53188769/71976530-6f82de80-320e-11ea-9483-c0ff5bb155ce.png)

Such that the equilibrium Fermi energy at the contacts could be determined by:

![image](https://user-images.githubusercontent.com/53188769/71976542-76115600-320e-11ea-9e79-00d42406f52f.png)

The solution to the EBIC model thus reduces to numerically finding ϕ_Fn, ϕ_Fp, and ψ.
The following boundary conditions were applied to find such solution. 
At the cross section surface where the electron beam impacts (Figure 1 in the manuscipt), a non-conductive boundary condition was applied. Carrier flux is determined by a current J_rec, and the total current into the boundary is zero σ_n ∇ϕ_Fn+σ_p ∇ϕ_Fp=0, thus:

![image](https://user-images.githubusercontent.com/53188769/71976548-7b6ea080-320e-11ea-8dce-bcea2271b9af.png)

The recombination current is given by:

![image](https://user-images.githubusercontent.com/53188769/71976561-832e4500-320e-11ea-8eda-e1549b980dfe.png)

Where U_s represents the surface recombination rate, S_n0,S_p0 represents the capture velocities for electrons and holes, and serves as a measure of how much surface recombination takes place at a mid-gap defect [1]. The electrostatic potential was set as a neutral surface ∇ψ=0, or a fixed amount of band bending was introduced by setting ψ=〖ψ_0+ψ〗_sbb at such surface.
At the collecting boundaries, n- and p-type sides, the electrostatic potential was set to the equilibrium Fermi energy difference ψ_0=ϕ_F0-ϕ_i, while the quasi Fermi energies were set to zero to force perfect collection of carriers at the contacts in a short-circuit condition:

![image](https://user-images.githubusercontent.com/53188769/71976573-888b8f80-320e-11ea-990a-c3341eab88ef.png)

Lastly, at symmetry boundary the material is assumed semi-infinite and no current should flow:


![image](https://user-images.githubusercontent.com/53188769/71976583-8fb29d80-320e-11ea-9d28-16d2993f2273.png)

The solution to ϕ_Fn, ϕ_Fp, and ψ was found using Matlab’s Partial Differential Equation Toolbox. The EBIC signal was then found by calculating the total current flowing on the n- and p-type contacts as:

![image](https://user-images.githubusercontent.com/53188769/71976596-97724200-320e-11ea-950b-877897ac2f84.png)


The discrepancy between the current densities in the emitter and the base contacts was used as a measure for model accuracy. When meshing or cell dimensions are chosen to be too coarse these values will differ substantially, thus indicating neutrality was not fully achieved.
