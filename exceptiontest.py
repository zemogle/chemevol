
def test_initial_galaxy_params():
	
	for run,keys in initial_galaxy_params.items():
		for k in init_keys:
			try:	
				dummy = initial_galaxy_params[run][k]
			except KeyError:
				print("Oops key %r is missing in %r" % (k,run))
			else:

			# check gasmass, gamma, inflows and outflows are numbers:
				if (k == 'gasmass_init') or (k == 'inflows') \
					or (k == 'outflows') or (k == 'gamma'):
					try:
						dummy = int(initial_galaxy_params[run][k])
					except ValueError:
						print("Oops we were expecting a number in %r:%r" % (run,k))

			# check SFH is a string :
				if (k == 'SFH'):
					if not isinstance(initial_galaxy_params[run][k], basestring):
						raise TypeError("Oops %r:%r should be a string" % (run,k))

			# check dust_source options correct
				if (k == 'dust_source'):
					dummy = initial_galaxy_params[run][k]
					if not ((dummy == 'SN') or (dummy == 'LIMS') or \
						   (dummy == 'LIMS+SN') or (dummy == 'GG') or\
						   (dummy == 'ALL')):
						raise ValueError("Oops double check %r:%r" % (run,k))

			# check IMF_fn source options correct
				if (k == 'IMF_fn'):
					dummy = initial_galaxy_params[run][k]
					if not ((dummy == 'Chab') or (dummy == 'chab') or (dummy == 'c') or \	
							(dummy == 'TopChab') or (dummy == 'topchab') or (dummy == 'tc') \
						   (dummy == 'Kroup') or (dummy == 'Salp') or \
						   ):
						raise ValueError("Oops check %r in %r" % (k,run))