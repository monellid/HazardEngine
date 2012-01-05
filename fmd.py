class TruncatedGutenbergRichter:
	"""
	Class defining Truncated Gutenberg-Richter frequency magnitude distribution.
	"""
	
	def __init__(self,a_value,b_value,min_mag,max_mag,bin_width):
		"""
		Define Gutenberg-Richter frequency magnitude distribution.
		a_value: cumulative a value (10**a_value is the rate of earthquakes
		[# of earthquakes / year] with magniutde greater or equal to 0)
		b_value: b value
		min_mag: minimum magnitude
		max_mag: maximum magnitude
		bin_width: bin width for occurrence rate calculation
		"""
		# TODO: check b value is greater than zero
		# TODO: check max_mag is greater or equal than min_mag
		self.a_value = a_value
		self.b_value = b_value
		self.min_mag = min_mag
		self.max_mag = max_mag
		self.bin_width = bin_width

	def getAnnualOccurrenceRates(self):
		"""
		Compute annual occurrence rates.
		"""
		min = round(self.min_mag / self.bin_width) * self.bin_width
		max = round(self.max_mag / self.bin_width) * self.bin_width
		if min != max:
			min = min + self.bin_width / 2
			max = max - self.bin_width / 2
		n_bins = int((max - min) / self.bin_width) + 1
		rates = []
		for i in range(n_bins):
			mag = min + i * self.bin_width
			rate = 10**(self.a_value - self.b_value * (mag - self.bin_width / 2)) - 10**(self.a_value - self.b_value * (mag + self.bin_width / 2))
			point = (mag,rate)
			rates.append(point)
		return rates

	def printAnnualOccurrenceRates(self):
		"""
		Print annual occurrence rates to standard output.
		"""
		for point in self.getAnnualOccurrenceRates():
			print point

