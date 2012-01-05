#!/usr/bin/python

from location import *
from fault import *

if __name__=='__main__':
	
	loc1 = Location(0.0,0.0,2.0)
	loc2 = Location(0.0,1.0,0.0)
	distance = 10.0
	
	locations = loc1.getLocations(loc2,distance)
	for loc in locations:
		print loc.longitude, loc.latitude, loc.depth