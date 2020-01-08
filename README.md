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
