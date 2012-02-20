#!/usr/bin/python

from scipy.stats import truncnorm
from scipy.stats import norm

def getProbabilityOfExceedance(value,mean,standard_deviation,truncation_level):
	"""
	Computes probability of exceeding 'value', by considering
	a normal distribution defined by 'mean', and 'standard_deviation'.
	If truncation_level is not None (must be greater than zero), a
	truncated normal is considered.
	truncation_level represents the number of standard deviations
	used for truncation.
	"""
	#TODO: check truncation_level is None or > 0
	# map 'value' to standard normal distribution
	value = (value - mean) / standard_deviation
	
	if truncation_level is not None:
		rv = truncnorm(-truncation_level,truncation_level)
	else:
		rv = norm()
		
	# computes probability of exceedance, by computing the
	# survival function, that is 1 - CDF (1 minus the cumulative
	# distribution function)
	return rv.sf(value)
	
# test 1
value = -2.995732273553991
mean = -0.7872268528578843
standard_deviation = 0.5962393527251486
truncation_level = 2.0
expected_poe = 1.0
computed_poe = getProbabilityOfExceedance(value,mean,standard_deviation,truncation_level)
print("Test 1. Expected poe: %s, computed poe: %s" % (expected_poe,computed_poe))

# test 2
value = -0.6931471805599453
mean = -0.7872268528578843
standard_deviation = 0.5962393527251486
truncation_level = 2.0
expected_poe = 0.43432352175355504
computed_poe = getProbabilityOfExceedance(value,mean,standard_deviation,truncation_level)
print("Test 2. Expected poe: %s, computed poe: %s" % (expected_poe,computed_poe))

# test 3
value = 0.6931471805599453
mean = -0.7872268528578843
standard_deviation = 0.5962393527251486
truncation_level = 2.0
expected_poe = 0.0
computed_poe = getProbabilityOfExceedance(value,mean,standard_deviation,truncation_level)
print("Test 3. Expected poe: %s, computed poe: %s" % (expected_poe,computed_poe))

# test 4
value = 0.6931471805599453
mean = -0.7872268528578843
standard_deviation = 0.5962393527251486
truncation_level = None
expected_poe = 0.006516701082128207
computed_poe = getProbabilityOfExceedance(value,mean,standard_deviation,truncation_level)
print("Test 4. Expected poe: %s, computed poe: %s" % (expected_poe,computed_poe))