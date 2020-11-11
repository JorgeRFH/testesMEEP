#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 19:06:16 2020
mmi Si3N4 - optimization- . FUNCIONA!
@author: jorge
"""
from __future__ import division

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# parameters of geometry
# El valor  optimo basado en la formula es 543 um (calculado en mmiSi3N4-w20um-modalanalysys.py)
# Sobre el valor coupler_length y run_until:
#   coupler length 54 y 104, suficiente con run_until 500/fcen
#   coupler length 204, se requiere 900/fcen


coupler_length = 453 
coupler_width = 20  # the right value  is 10
guide_length = 10   # the right value  is 10
waveguide_width = 3
separation_output = 10
circle_diameter = 0.4  # diametro de 400nm
spacing_circles = 0.1  # distancia entre circulos de 100nm


# materials
Si3N4_neff_TM = mp.Medium(index = 1.652312) # effective index TM of Slab Si3N4
SiO2 = mp.Medium(index = 1.457)


resolution = 18 # pixels/um, generalmente potencia de 2
cell_x = 2*guide_length+coupler_length
cell_y = coupler_width*2
cell = mp.Vector3(cell_x,cell_y,0)


dpml = 2.0
pml_layers = [mp.PML(dpml)]
symmetries = [mp.Mirror(mp.Y)] # Posiblemente quitemos la simetria en simulacion de opt

# Idea de la geometria de los cilindros
#       o o o o o  ...  o o 
#   N   ...           ...
#       o o o o o  ...  o o
#      
#            M bolinhas




N = coupler_width//(circle_diameter+spacing_circles)
M = (coupler_length-spacing_circles)//(circle_diameter+spacing_circles)

geometry_circles = []
#matrices = np.zeros((int(N),int(M)))
matrices = np.ones((int(N),int(M)))
i=0
k=0
while (i<M):
    while (k<N):
        geometry_circles = geometry_circles + [mp.Cylinder(radius=circle_diameter/2, 
                                           center=mp.Vector3(-0.5*coupler_length+(i+1)*spacing_circles+circle_diameter*0.5+i*circle_diameter,
                                                             0.5*((N-1)*circle_diameter+(N-1)*spacing_circles)-k*(circle_diameter+spacing_circles),0),
                                           material=Si3N4_neff_TM if matrices[k][i]==1 else SiO2)]
        k = k+1
    i = i+1
    k=0


geometry = [mp.Block(size=mp.Vector3(mp.inf,waveguide_width,mp.inf),
                     center=mp.Vector3(),
                     material=Si3N4_neff_TM)]


# aqui define la fuente gaussiana, pensando en capturar el espectro
fcen = 1/0.632 # pulse center freq (es 1.58 para wavelength de 632nm )
df = 0.1*fcen    # pulse width freq
rot_angle = np.radians(0)
kpoint = mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle)
fsrc = 1/0.632 # frequency of eigenmode or constant-amplitude source
bnum = 1    # band number of eigenmode previo era 1

sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=df),
                                  center=mp.Vector3(-0.5*cell_x+2*dpml,0),# (0,0)
                                  size=mp.Vector3(y=3*waveguide_width),
                                  direction=mp.NO_DIRECTION,
                                  eig_kpoint=kpoint,
                                  eig_band=bnum,
                                  eig_parity=mp.EVEN_Y if rot_angle == 0 else mp.ODD_Z,
                                  eig_match_freq=True)]




sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    default_material=SiO2,
                    symmetries = symmetries,
                    resolution=resolution)


nfreq = 8  # number of freq at which to compute flux. ideal 50 or 100

# reflected flux
refl_fr = mp.FluxRegion(center=mp.Vector3(-0.5*cell_x+2*dpml+1,0,0), size=mp.Vector3(0,2*waveguide_width,0))                            
refl = sim.add_flux(fcen, df, nfreq, refl_fr)

# transmitted flux
tran_fr = mp.FluxRegion(center=mp.Vector3(-0.5*coupler_length-1,0,0), size=mp.Vector3(0,2*waveguide_width,0))
tran = sim.add_flux(fcen, df, nfreq, tran_fr)


pt = mp.Vector3(-0.5*coupler_length-1,0)
sim.run(until_after_sources=mp.stop_when_fields_decayed(20/fcen, mp.Ez, pt, 1e-3))

straight_refl_data = sim.get_flux_data(refl)
straight_tran_flux = mp.get_fluxes(tran)


# Aqui arranca la segunda simulacion, con la geometria completa
sim.reset_meep()
geometry = [mp.Block(mp.Vector3(coupler_length,coupler_width,mp.inf),
                      center=mp.Vector3(),
                      material=Si3N4_neff_TM),
            mp.Block(mp.Vector3(guide_length,waveguide_width,mp.inf),
                      center=mp.Vector3(-(guide_length+coupler_length)/2,0),
                      material=Si3N4_neff_TM),
            mp.Block(mp.Vector3(guide_length,waveguide_width,mp.inf),
                      center=mp.Vector3((guide_length+coupler_length)/2,separation_output/2),
                      material=Si3N4_neff_TM),
            mp.Block(mp.Vector3(guide_length,waveguide_width,mp.inf),
                      center=mp.Vector3((guide_length+coupler_length)/2,-separation_output/2),
                      material=Si3N4_neff_TM)] + geometry_circles

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    default_material=SiO2,
                    symmetries = symmetries,
                    resolution=resolution)

# reflected flux
refl = sim.add_flux(fcen, df, nfreq, refl_fr)

# transmitted flux1
tran_fr1 = mp.FluxRegion(center=mp.Vector3(0.5*cell_x-dpml,0.5*separation_output,0), 
                        size=mp.Vector3(0,1*waveguide_width,0))
tran1 = sim.add_flux(fcen, df, nfreq, tran_fr1)

# transmitted flux2
tran_fr2 = mp.FluxRegion(center=mp.Vector3(0.5*cell_x-dpml,-0.5*separation_output,0), 
                        size=mp.Vector3(0,1*waveguide_width,0))
tran2 = sim.add_flux(fcen, df, nfreq, tran_fr2)

# for normal run, load negated fields to subtract incident from refl. fields
sim.load_minus_flux_data(refl, straight_refl_data)


# Este comando abajo es interesante pues guarda para el gif y tbn actualiza los valores de tx y rx
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.at_every(20*0.632, mp.output_efield_z)),# TM: _z, TE: _y
        until=1500/fcen)

# ESte codigo es para obtener los kpoints de cada modo, debe de salir igual al lumerical:
# Al poner [1,2,3] dara los valores kpoints: TM0 (8 valores), TE0 (8 valores), TM1 (8 valores)
# res = sim.get_eigenmode_coefficients(tran,
#                                      [1,2,3], # Calculará los valores kpoints
#                                      eig_parity=mp.NO_PARITY,#mp.EVEN_Y+mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
#                                      direction=mp.NO_DIRECTION,
#                                      kpoint_func=lambda f,n: kpoint)



mmi_refl_flux = mp.get_fluxes(refl)
mmi_tran1_flux = mp.get_fluxes(tran1)
mmi_tran2_flux = mp.get_fluxes(tran2)

flux_freqs = mp.get_flux_freqs(refl)

eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric) # captura valor dielectrico
ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez) # captura el campo propagante
plt.figure()
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary') # dibuja geometria dielectrico
#plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
plt.axis('off') # quita los ejes con sus numeros 
plt.title('Total length = um')
plt.show()


# Para graficar la curva:
wl = []
Rs = []
T1s = []
T2s = []

for i in range(nfreq):
    wl = np.append(wl, 1/flux_freqs[i])
    Rs = np.append(Rs,-mmi_refl_flux[i]/straight_tran_flux[i])
    T1s = np.append(T1s,mmi_tran1_flux[i]/straight_tran_flux[i])    
    T2s = np.append(T2s,mmi_tran2_flux[i]/straight_tran_flux[i])    

if mp.am_master():
    plt.figure()
    plt.plot(wl,Rs,'bo-',label='reflectance')
    plt.plot(wl,T1s,'ro-',label='transmittance 1 y 0.5*S')
    plt.plot(wl,T2s,'k^-',label='transmittance 2 y -0.5*S')
    plt.plot(wl,1-Rs-T1s-T2s,'go-',label='loss')
    plt.axis([min(wl),max(wl),0,1.2])
    plt.xlabel("wavelength (μm)")
    plt.legend(loc="upper right")
    plt.title("resolution=16")
    plt.show()