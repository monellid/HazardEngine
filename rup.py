class ProbEqkRupture:
	"""
	Class defining probabilistic earthquake rupture.
	"""
	
	def __init__(self,magnitude,rupture_surface,hypocenter,rake,tectonic_region_type,probability_occurrence):
		"""
		magnitude: earthquake magnitude
		rupture_surface: 2D array representing rupture surface mesh
		hypocenter: hypocenter location (longitude,latitude,depth)
		rake: rake angle
		tectonic_region_type: tectonic region type
		probability_occurrence: probability of occurrence 
		"""
		self.magnitude = magnitude
		self.rupture_surface = rupture_surface
		self.hypocenter = hypocenter
		self.rake = rake
		self.tectonic_region_type = tectonic_region_type
		self.probability_occurrence = probability_occurrence