import openmc
import os

os.makedirs("output", exist_ok=True)
os.chdir("output")

##############################################
                # Materials #
##############################################

uo2 = openmc.Material(name="uo2")
uo2.add_nuclide('U235',0.04)
uo2.add_nuclide('U238',0.96)
uo2.add_nuclide('O16',2.0)
uo2.set_density('g/cm3',10.4)
uo2.depletable = True

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

uo2_outer_radius = openmc.ZCylinder(r=0.4096)
mox_outer_radius = openmc.ZCylinder(r=0.4096)
clad_inner_radius = openmc.ZCylinder(r=0.4179)
clad_outer_radius = openmc.ZCylinder(r=0.4750)

uo2_fuel_region = -uo2_outer_radius
mox_fuel_region = -mox_outer_radius
gap_region = +mox_outer_radius & -clad_inner_radius
clad_region = +clad_inner_radius & -clad_outer_radius

uo2_fuel_cell = openmc.Cell(1, name='fuel')
uo2_fuel_cell.fill = uo2
uo2_fuel_cell.region = uo2_fuel_region

mox_fuel_cell = openmc.Cell(1, name='fuel')
mox_fuel_cell.fill = mox
mox_fuel_cell.region = mox_fuel_region

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

uo2_universe = openmc.Universe(cells=[uo2_fuel_cell, gap_cell, clad_cell, fuel_water_cell])
mox_universe = openmc.Universe(cells=[mox_fuel_cell, gap_cell, clad_cell, fuel_water_cell])

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

mod_universe = openmc.Universe(cells=[waterrod_inner_cell, wr_clad_cell, waterrod_outer_cell])

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

########## Assembly ##########

import numpy as np

def create_assembly(fuel_universe, mod_universe, pitch, name):
    layout = np.full((17,17), 'F')
    moderator_positions = [
        (0,0),(0,4),(0,8),(0,12),(0,16),
        (4,0),(4,4),(4,8),(4,12),(4,16),
        (8,0),(8,4),(8,8),(8,12),(8,16),
        (12,0),(12,4),(12,8),(12,12),(12,16),
        (16,0),(16,4),(16,8),(16,12),(16,16),
    ]

    for i, j in moderator_positions:
        layout[i,j] = 'M'

    universe_map = {'F': mox_universe, 'M': mod_universe}
    universes = [[universe_map[cell] for cell in row] for row in layout]

    full_pitch = pitch * 17
    lattice = openmc.RectLattice(name='UO2 Assembly')
    lattice.pitch = (pitch,pitch)
    lattice.lower_left = [-full_pitch/2, -full_pitch/2]
    lattice.universes = universes
    return lattice

uo2_assembly = create_assembly(uo2_universe, mod_universe, pitch, 'UO2 Assembly')
mox_assembly = create_assembly(mox_universe, mod_universe, pitch, 'MOX Assembly')

########## Sleeve and Outer Water ##########

full_pitch = pitch * 17
assembly_region = -openmc.model.RectangularPrism(width=full_pitch, height=full_pitch, origin=(0,0))
uo2_assembly_cell = openmc.Cell(name='uo2 assembly cell', fill=uo2_assembly, region=assembly_region)
mox_assembly_cell = openmc.Cell(name='mox assembly cell', fill=mox_assembly, region=assembly_region)

sleeve_thickness = 0.1

uo2_assembly_sleeve = openmc.Cell(name='full assembly sleeve')
uo2_assembly_sleeve.region = -openmc.model.RectangularPrism(width=full_pitch+2*sleeve_thickness, height=full_pitch+2*sleeve_thickness) & ~assembly_region
uo2_assembly_sleeve.fill = zircaloy

uo2_assembly_outer_water = openmc.Cell(name='outer water')
uo2_assembly_outer_water.region = ~uo2_assembly_sleeve.region & ~assembly_region & -openmc.model.RectangularPrism(width=full_pitch+2*sleeve_thickness+1, height=full_pitch+2*sleeve_thickness+1, boundary_type='reflective')
uo2_assembly_outer_water.fill = water

mox_assembly_sleeve = openmc.Cell(name='full assembly sleeve')
mox_assembly_sleeve.region = -openmc.model.RectangularPrism(width=full_pitch+2*sleeve_thickness, height=full_pitch+2*sleeve_thickness) & ~assembly_region
mox_assembly_sleeve.fill = zircaloy

mox_assembly_outer_water = openmc.Cell(name='outer water')
mox_assembly_outer_water.region = ~mox_assembly_sleeve.region & ~assembly_region & -openmc.model.RectangularPrism(width=full_pitch+2*sleeve_thickness+1, height=full_pitch+2*sleeve_thickness+1, boundary_type='reflective')
mox_assembly_outer_water.fill = water

uo2_assembly_universe = openmc.Universe(cells=[uo2_assembly_cell, uo2_assembly_sleeve, uo2_assembly_outer_water])

mox_assembly_universe = openmc.Universe(cells=[mox_assembly_cell, mox_assembly_sleeve, mox_assembly_outer_water])

core_lattice = openmc.RectLattice(name='2x1 Core Lattice')
core_lattice.pitch = (full_pitch + 2*sleeve_thickness, full_pitch + 2*sleeve_thickness)
core_lattice.lower_left = [-core_lattice.pitch[0], -core_lattice.pitch[1]/2]

core_lattice.universes = [[uo2_assembly_universe, mox_assembly_universe]]

core_width = 2*core_lattice.pitch[0]
core_height = core_lattice.pitch[1]
core_region = -openmc.model.RectangularPrism(width=core_width, height=core_height, boundary_type='reflective')
core_cell = openmc.Cell(name='Core Cell', fill=core_lattice, region=core_region)

root_universe = openmc.Universe(name='root universe', cells=[core_cell])

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

##############################################
                # Settings #
##############################################

source_x_halfwidth = full_pitch + 2*sleeve_thickness
source_y_halfwidth = full_pitch/2 + sleeve_thickness

bounds = [-source_x_halfwidth, -source_y_halfwidth, -1, source_x_halfwidth, source_y_halfwidth, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)

settings = openmc.Settings()
settings.source = openmc.Source(space=uniform_dist)
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.run_mode = 'eigenvalue'
settings.export_to_xml()

##############################################
                # Tallies #
##############################################

mesh = openmc.RegularMesh()
mesh.dimension = [50, 100, 1] 
mesh.lower_left = [-18, -12, -1]
mesh.upper_right = [18, 12, 1]

uo2_cell_filter = openmc.CellFilter(uo2_fuel_cell)

uo2_fuel_tally = openmc.Tally(name='U235_fuel_reactions')
uo2_fuel_tally.filters = [uo2_cell_filter]
uo2_fuel_tally.nuclides = ['U235']
uo2_fuel_tally.scores = ['total', 'fission', 'absorption', '(n,gamma)']

mox_cell_filter = openmc.CellFilter(mox_fuel_cell)

mox_fuel_tally = openmc.Tally(name='Pu239_fuel_reactions')
mox_fuel_tally.filters = [mox_cell_filter]
mox_fuel_tally.nuclides = ['Pu239','Pu240','Pu241']
mox_fuel_tally.scores = ['total', 'fission', 'absorption', '(n,gamma)']

flux_tally = openmc.Tally(name='flux')
flux_tally.filters = [openmc.MeshFilter(mesh)]
flux_tally.scores = ['flux']

rr_tally = openmc.Tally(name='reaction_rates')
rr_tally.scores = ['fission', 'absorption']

tallies = openmc.Tallies([uo2_fuel_tally, mox_fuel_tally, rr_tally, flux_tally])
tallies.export_to_xml()

##############################################
                # Plotting #
##############################################

########### Geometry Vizualization ##########

plot = openmc.Plot()
plot.origin = (0,0,0)
plot.width = (44,22)
plot.pixels = (1000,1000)
plot.color_by = 'material'
plot.basis = 'xy'
plot.filename = 'core_assembly_plot'

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
print(sp.keff)

flux_tally = sp.get_tally(name='flux')

df = flux_tally.get_pandas_dataframe()
flux = df['mean'].values.reshape((50, 100))

plt.imshow(flux, origin='lower')
plt.colorbar(label='Flux')
plt.title('Flux Distribution')
plt.show()








