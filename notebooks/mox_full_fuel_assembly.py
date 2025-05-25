import openmc
import os

os.makedirs("output", exist_ok=True)
os.chdir("output")

##############################################
                # Materials #
##############################################

mox = openmc.Material(name='mixed oxide fuel')
mox.add_nuclide('U235',0.0023)
mox.add_nuclide('U238',0.9177)
mox.add_nuclide('O16',2.0)
mox.add_nuclide('Pu239',0.052)
mox.add_nuclide('Pu240',0.020)
mox.add_nuclide('Pu241',0.008)
mox.set_density('g/cm3',10.4)

zircaloy = openmc.Material(name='zircaloy')
zircaloy.add_nuclide('Zr90',7.2758e-3)
zircaloy.set_density('g/cm3',6.55)
zircaloy.depletable = False

steel = openmc.Material(name='stainless steel')
steel.add_element('C', 0.08, percent_type='wo')
steel.add_element('Si', 1.00, percent_type='wo')
steel.add_element('P', 0.045, percent_type='wo')
steel.add_element('S', 0.030, percent_type='wo')
steel.add_element('Mn', 2.00, percent_type='wo')
steel.add_element('Cr', 20.0, percent_type='wo')
steel.add_element('Ni', 11.0, percent_type='wo')
steel.add_element('Fe', 65.845, percent_type='wo')
steel.set_density('g/cm3', 8.00)

water = openmc.Material(name="water")
water.add_nuclide('H1',2.0)
water.add_nuclide('O16',1.0)
water.set_density('g/cm3',1.0)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([mox, zircaloy, water, steel])
materials.export_to_xml()

##############################################
                # Geometry #
##############################################

########## Fuel Pins ##########

fuel_outer_radius = openmc.ZCylinder(r=0.4096)
clad_inner_radius = openmc.ZCylinder(r=0.4179)
clad_outer_radius = openmc.ZCylinder(r=0.4750)

fuel_region = -fuel_outer_radius
gap_region = +fuel_outer_radius & -clad_inner_radius
clad_region = +clad_inner_radius & -clad_outer_radius

fuel_cell = openmc.Cell(1, name='fuel')
fuel_cell.fill = mox
fuel_cell.region = fuel_region

gap_cell = openmc.Cell(name='air gap')
gap_cell.region = gap_region

clad_cell = openmc.Cell(name='clad')
clad_cell.fill = zircaloy
clad_cell.region = clad_region

pitch = 1.26

box = openmc.model.RectangularPrism(width=pitch, height=pitch,
                                    boundary_type='reflective')
type(box)

fuel_water_region = +clad_outer_radius & -box

fuel_water_cell = openmc.Cell(name='fuel water')
fuel_water_cell.fill = water
fuel_water_cell.region = fuel_water_region

fuel_pin_universe = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, fuel_water_cell])

########## Guide Pins ##########

waterrod_outer_radius = openmc.ZCylinder(r=0.4179)
wr_clad_outer_radius = openmc.ZCylinder(r=0.4750)

waterrod_inner_region = -waterrod_outer_radius
wr_clad_region = +waterrod_outer_radius & -wr_clad_outer_radius
waterrod_outer_region = +wr_clad_outer_radius & -box

waterrod_inner_cell = openmc.Cell(name='waterrod_inner')
waterrod_inner_cell.fill = water
waterrod_inner_cell.region = waterrod_inner_region

wr_clad_cell = openmc.Cell(name='wr_clad')
wr_clad_cell.fill = zircaloy
wr_clad_cell.region = wr_clad_region

waterrod_outer_cell = openmc.Cell(name='waterrod_outer')
waterrod_outer_cell.fill = water
waterrod_outer_cell.region = waterrod_outer_region

water_pin_universe = openmc.Universe(cells=[waterrod_inner_cell, wr_clad_cell, waterrod_outer_cell])

########## Instrument Pin ##########

ipin_inner_radius = openmc.ZCylinder(r=0.4750)
ipin_outer_radius = openmc.ZCylinder(r=0.5500)

ipin_steel = -ipin_inner_radius
ipin_clad = +ipin_inner_radius & -ipin_outer_radius
ipin_water = +ipin_outer_radius & -box

ipin_steel_cell = openmc.Cell(name='ipin_rod')
ipin_steel_cell.fill = steel
ipin_steel_cell.region = ipin_steel

ipin_clad_cell = openmc.Cell(name='ipin_clad')
ipin_clad_cell.fill = zircaloy
ipin_clad_cell.region = ipin_clad

ipin_water_cell = openmc.Cell(name='ipin_water')
ipin_water_cell.fill = water
ipin_water_cell.region = ipin_water

ipin_universe = openmc.Universe(cells=[ipin_steel_cell, ipin_clad_cell, ipin_water_cell])

##############################################
                # Assembly #
##############################################

full_pitch = pitch * 17
assembly = openmc.RectLattice(name='Full Assembly')
assembly.pitch = (pitch,pitch)
assembly.lower_left = [-full_pitch/2, -full_pitch/2]

assembly.universes = [
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, ipin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe,fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
[fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, water_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe, fuel_pin_universe],
]

assembly_region = -openmc.model.RectangularPrism(width=full_pitch, height=full_pitch, origin=(0,0))
assembly_cell = openmc.Cell(name='full assembly cell', fill=assembly, region=assembly_region)

sleeve_thickness = 0.1
assembly_sleeve = openmc.Cell(name='full assembly sleeve')
assembly_sleeve.region = -openmc.model.RectangularPrism(width=full_pitch+2*sleeve_thickness, height=full_pitch+2*sleeve_thickness) & ~assembly_cell.region
assembly_sleeve.fill = zircaloy

assembly_outer_water = openmc.Cell(name='outer water')
assembly_outer_water.region = ~assembly_sleeve.region & ~assembly_cell.region & -openmc.model.RectangularPrism(width=full_pitch+2*sleeve_thickness+1, height=full_pitch+2*sleeve_thickness+1, boundary_type='reflective')
assembly_outer_water.fill = water

full_assembly_universe = openmc.Universe(cells=[assembly_cell, assembly_sleeve, assembly_outer_water])

geometry = openmc.Geometry(full_assembly_universe)
geometry.export_to_xml()

##############################################
                # Settings #
##############################################

source = openmc.IndependentSource(space=openmc.stats.Point((0,0,0)))

settings = openmc.Settings()
settings.source = source
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.export_to_xml()

##############################################
                # Tallies #
##############################################

mesh = openmc.RegularMesh()
mesh.dimension = [50, 50, 1] 
mesh.lower_left = [-12, -12, -1]
mesh.upper_right = [12, 12, 1]

cell_filter = openmc.CellFilter(fuel_cell)

fuel_tally = openmc.Tally(name='fuel_reactions')
fuel_tally.filters = [cell_filter]
fuel_tally.nuclides = ['U235']
fuel_tally.scores = ['total', 'fission', 'absorption', '(n,gamma)']

flux_tally = openmc.Tally(name='flux')
flux_tally.filters = [openmc.MeshFilter(mesh)]
flux_tally.scores = ['flux']

rr_tally = openmc.Tally(name='reaction_rates')
rr_tally.scores = ['fission', 'absorption']

tallies = openmc.Tallies([fuel_tally, flux_tally, rr_tally])
tallies.export_to_xml()

##############################################
                # Plotting #
##############################################

########### Geometry Vizualization ##########

plot = openmc.Plot()
plot.origin = (0,0,0)
plot.width = (15,15)
plot.pixels = (1000,1000)
plot.color_by = 'material'
plot.basis = 'xy'
plot.filename = 'fuel_assembly_plot'

plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()

##############################################
                # Run #
##############################################

openmc.run()

##############################################
                # Post-Processing #
##############################################

import openmc
import matplotlib.pyplot as plt
import numpy as np

sp = openmc.StatePoint('statepoint.100.h5')

for tally in sp.tallies.values():
    print(tally)

flux_tally = sp.get_tally(name='flux')

df = flux_tally.get_pandas_dataframe()
flux = df['mean'].values.reshape((50, 50))

plt.imshow(flux, origin='lower')
plt.colorbar(label='Flux')
plt.title('Flux Distribution')
plt.show()








