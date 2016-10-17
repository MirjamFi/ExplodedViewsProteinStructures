# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_movie.py

# ''' load complex from PDB'''
# pdbName = '3oaa'
# pdbFile = pdbName + '.pdb1'
# cmd.fetch(pdbName, type='pdb1')
# cmd.load('C:\Users\Figaschewski\Dropbox\Masterarbeit\TPC2_M484L_wActivator_56ns\TPC2_M484L_wActivator_56ns.pdb', 'test')

'''
	Create an movie of an exploded view for a given protein structure (PDB)
'''
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time


def translate_selection(originXYZ, transXYZ, transname, factor, f):
	''' DESCRIPTION:
		translate an object relative to complex of origin using center of mass
	'''
	## vector between COMs to translate chain
	dist=math.sqrt((transXYZ[0]-originXYZ[0])**2 +
					(transXYZ[1]-originXYZ[1])**2 +
					(transXYZ[2]-originXYZ[2])**2)

	vector=((transXYZ[0]-originXYZ[0])/dist, 
			(transXYZ[1]-originXYZ[1])/dist, 
			(transXYZ[2]-originXYZ[2])/dist)
	trans_vec = [x * factor for x in vector] 

	## translate chain with vector
	cmd.frame(f)
	cmd.translate(trans_vec, object=transname, camera=0)
	cmd.mview('store', object=transname)
	cmd.mview('interpolate', object=transname)
	

def calc_COM(transname):
	'''DESCRIPTION:
		calculate COM of an object 
	'''
	cenma.com(transname, state = 1)
	cmd.zoom(transname + '_COM')
	pos = cmd.get_position(transname + '_COM')
	cmd.delete(transname + '_COM')
	cmd.set_view (
			'0.453708082,   -0.681963921,   -0.573648989,\
			-0.331614792,    0.468286544,   -0.818986893,\
			0.827151597,    0.561811924,   -0.013684246,\
			0.000000000,    0.000000000, -981.596557617,\
			-68.949493408,  -24.481521606,    5.007858276,\
			773.898193359, 1189.294921875,  -20.000000000' )
	return pos
		
	 
def explosion(pdbName):
	''' DESCRIPTION:
		create a movie for an exploded view of given protein
	'''
	
	'''setup of PyMOL'''
	start_time = time.clock()
	## case sensitive for chain ids
	cmd.set('ignore_case', 'off') 
	
	cmd.remove('solvent')
	cmd.show('spheres', 'organic')
	util.cbc(selection= pdbName)
	
	##get chains of complex
	chains = cmd.get_chains(pdbName)

	''' get chains and ligands'''
	##get all ligands (organic) from complex'''
	stored.ligands = []

	## create the original selection of all organic atoms 
	cmd.select('allOrg', 'org')

	##' temporary selection
	cmd.select('temp', 'none')

	## while initial selection is not empty, 'pop' from it and query the atom's 
	## molecule's name
	while cmd.count_atoms('allOrg') != 0:

		## pop--peek, rather
		cmd.select('temp', 'first allOrg')
		
		## store the residue name
		cmd.iterate('temp', 'stored.ligands.append(resn)')
		
		## remove by molecule
		cmd.select('allOrg', 'allOrg and not (bm. temp)')
		
	cmd.delete('allOrg')
	cmd.delete('temp')

	'''initialize movie'''
	cmd.set('matrix_mode', 1)
	cmd.set('movie_panel', 1)
	cmd.set('scene_buttons', 1)
	cmd.set('cache_frames', 1)
	cmd.config_mouse('three_button_motions', 1)

	frameNum = 50
	frames = str(frameNum * (len(cmd.get_chains()) + len(stored.ligands))+50)
	cmd.mset('1 x' + frames)
	cmd.set_view (
		'0.453708082,   -0.681963921,   -0.573648989,\
		-0.331614792,    0.468286544,   -0.818986893,\
		 0.827151597,    0.561811924,   -0.013684246,\
		 0.000000000,    0.000000000, -981.596557617,\
	   -68.949493408,  -24.481521606,    5.007858276,\
	   773.898193359, 1189.294921875,  -20.000000000 ')

			
	''' Explosion: for each state create an object of every chain and for 
		related ligands. In this state calulate COM for state and according 
		chains and translate chain along vector between COMs of state and chain. 
		For every chain do the same with its ligands.
	'''

	## create models for all states of complex
	states = cmd.split_states(pdbName)

	com_num = ''
	for state in range(1,cmd.count_states()+1):
		## state = 0 would calculate the COM for all states at once
		if state == 0: 
			continue

		## name chain object in this state
		s = '00' + str(state)
		if state < 10:
			s = '0' + s
			
		## COM of state (all chains)
		complexXYZ = calc_COM(pdbName + '_' + s)

		## calculate COM for single chains in state
		cNames = []
		chainsCOMS = {}
		## COM for ligands
		ligandsCOMS = {}
		f = 1
		cmd.delete(pdbName)
		for c in chains:	
			## color each chain individually
			if cmd.count_states() > 1:
				cmd.color(get_colors.get_random_color(), pdbName + '_' + s + \
							'& chain ' + c)

			## create object for chain with name
			chainname = pdbName + '_' + s + '_'+c
			cNames.append(chainname)
			origin = pdbName+'_'+s+ ' & chain '+c
			cmd.extract(chainname, origin)
			
			## store object for movie
			cmd.frame(f)
			cmd.show('sticks', chainname)
			cmd.mview('store', object=chainname)
			
			## calculate and save COM of chain
			chainsCOMS[chainname] = calc_COM(chainname)	
			
			## create objects for ligands
			if stored.ligands:
				## select ligand on chain	
				selection = chainname + ' & organic'
				ligandname = chainname + '_org'

				## only translate if chain really has a ligand
				if cmd.count_atoms(selection) > 0:
					cmd.extract(ligandname, selection)
					## calculate COM of ligand
					ligandsCOMS[ligandname] = calc_COM(ligandname)
				
					## color binding site
					binding = 'byres (' + chainname + ' nto. 3 of organic)'
					cmd.select('inter', binding)
					cmd.color('red', 'inter')
					cmd.delete('inter')
		
		f = frameNum
		cmd.frame(f)
		## store chain objects for movie
		for cname in cNames:
			cmd.mview('store', object=cname)
		## store ligand objects for movie
		for ligand in ligandsCOMS.keys():
			cmd.mview('store', object=ligand)
		f = f + frameNum

		## translate chains iteratively 
		for chainname in cNames:
			## get dimensions of chain 
			#([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(chainname)
			## coordinates
			#[x,y,z] = (maxX-minX, maxY-minY, maxZ-minZ)
			## volume
			#vol = ((maxX-minX)*(maxY-minY)*(maxZ-minZ))
			#print [x,y,z], vol
			
			## COM of chain
			chainXYZ = chainsCOMS[chainname]

			## only translate chain if there are multiple chains
			if len(chains) > 1:	
				## translate chains and ligands
				for cname in cNames:
					cXYZ = chainsCOMS[cname]
					translate_selection(complexXYZ, cXYZ, cname, 20, f)
				for ligand in ligandsCOMS.keys():
					ligandXYZ = ligandsCOMS[ligand]
					translate_selection(complexXYZ, ligandXYZ, ligand, 20, f)

			f = f + frameNum
		## translate ligands once more for greater distance to related chains
		for ligand in ligandsCOMS.keys():
			ligandXYZ = ligandsCOMS[ligand]
			translate_selection(complexXYZ, ligandXYZ, ligand, 10, f)

		## delete state object
		cmd.delete(pdbName + '_'+ s)	
		

	print time.clock() - start_time, 'seconds'
	
cmd.extend('explosion', explosion)
