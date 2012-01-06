#!/usr/bin/python

from geo import Point
from fault import *

if __name__=='__main__':
	
	# test 0
	#p1 = Point(0.0,0.0)
	#p2 = Point(0.0,0.0,10.0)
	#distance = 2.0
	#points = p1.getEquallySpacedPoints(p2,distance)
	#for p in points:
	#	print p.longitude, p.latitude, p.depth
	
	# test 1
	#p1 = Point(0.0,0.0,2.0)
	#p2 = Point(0.0,1.0,0.0)
	#distance = 10.0

	#points = p1.getEquallySpacedPoints(p2,distance)
	#for p in points:
	#	print p.longitude, p.latitude, p.depth
	
	# test 2	
	p1 = Point(0.0,0.0)
	p2 = Point(0.0,1.0)
	p3 = Point(0.5,1.5)
	fault_trace = [p1,p2,p3]

	upper_seismo_depth = 2.0
	lower_seismo_depth = 15.0
	dip = 90.0
	mesh_spacing = 1.0

	fault_surface = SimpleFaultSurface(fault_trace,upper_seismo_depth,lower_seismo_depth,dip)

	surface = fault_surface.getFaultSurface(mesh_spacing)
	num_rows = surface.shape[0]
	num_cols = surface.shape[1]
	for i in range(num_rows):
		for j in range(num_cols):
			print surface[i,j].longitude, surface[i,j].latitude, surface[i,j].depth
		print '\n'