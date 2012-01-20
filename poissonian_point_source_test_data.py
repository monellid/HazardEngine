#!/usr/bin/python

from point import *
from geo import *
from fmd import *
from msr import *

point = Point(0.0,0.0)

strike = 45.0
dip = 45.0
rake = 0.0
weight_nodal_plane = 1.0
nodal_plane_pmf = [(strike,dip,rake,weight_nodal_plane)]

hypo = 8.0
weight_hypo = 1.0
hypo_depth_pmf = [(hypo,weight_hypo)]

upper_seismo_depth = 2.0
lower_seismo_depth = 16.0
rupt_aspect_ratio = 1.0

tectonic_region_type = 'Active Shallow Crust'

mag_scaling_rel = PeerTestMagAreaScalingRel()

time_span = 50.0

mesh_spacing = 1.0

a_value = 2.0
b_value = 1.0
min_mag = 5.0
max_mag = 6.0
bin_width = 1.0
freq_mag_dist = TruncatedGutenbergRichter(a_value,b_value,min_mag,max_mag,bin_width)

point_source = PoissonianPointSource(point,
									freq_mag_dist,
									nodal_plane_pmf,
									hypo_depth_pmf,
									upper_seismo_depth,
									lower_seismo_depth,
									mag_scaling_rel,
									rupt_aspect_ratio,
									tectonic_region_type,
									time_span,
									mesh_spacing)
									
num_rup = point_source.getNumRuptures()
print 'num ruptures:',num_rup

rupture = point_source.getRupture(0)

mag = rupture['magnitude']
strike = rupture['strike']
dip = rupture['dip']
rake = rupture['rake']
tectonic = rupture['tectonic']
hypocenter = rupture['hypocenter']
surface = rupture['surface']
rate = rupture['rate']
probability = rupture['probability']

print 'magnitude: ',mag
print 'strike: ',strike
print 'dip: ',dip
print 'rake: ',rake
print 'tectonic region type: ',tectonic
print 'hypocenter: %s, %s, %s'%(hypocenter.longitude,hypocenter.latitude,hypocenter.depth)
print 'surface:'
for i in range(surface.shape[0]):
	for j in range(surface.shape[1]):
		print '(%s, %s, %s) ' % (surface[i,j].longitude,surface[i,j].latitude,surface[i,j].depth)
	print '\n'
print 'rate: ',rate
print 'probability: ',probability
