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
'-0.346189618,   -0.600477457,   -0.720818758,\
     0.038613044,    0.758557022,   -0.650460541,\
     0.937369645,   -0.253015876,   -0.239418939,\
     0.000000000,    0.000000000, -1473.033813477,\
   -73.633697510,  -31.334754944,    1.869384766,\
  1161.351074219, 1784.716552734,  -20.000000000'	)
	return pos
		

def isColliding(sel1, sel2):
	''' DESCIPTION:
		calculate if two bounding boxes are colliding
	'''
	box1_min, box1_max = cmd.get_extent(sel1)
	box2_min, box2_max = cmd.get_extent(sel2)

	## min, max vertices of box 1
	x_min1 = box1_min[0]
	y_min1 = box1_min[1]
	z_min1 = box1_min[2]
	x_max1 = box1_max[0]
	y_max1 = box1_max[1]
	z_max1 = box1_max[2]
	
	## min, max vertices of box 2
	x_min2 = box2_min[0]
	y_min2 = box2_min[1]
	z_min2 = box2_min[2]
	x_max2 = box2_max[0]
	y_max2 = box2_max[1]
	z_max2 = box2_max[2]
	
	## test if any of the vertices of box 1 are in any of the faces of box 2 and vice versa
	isColliding = 	(x_min1 <= x_max2 and x_max1 >= x_min2) and \
					(y_min1 <= y_max2 and y_max1 >= y_min2) and \
					(z_min1 <= z_max2 and z_max1 >= z_min2)
	

	return isColliding
	 
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
	
	## calculate translation factor
	protDim_min, protDim_max = cmd.get_extent(pdbName)
	transFac = math.sqrt((protDim_max[0] - protDim_min[0]) ** 2 +
                     (protDim_max[1] - protDim_min[1]) ** 2 +
                     (protDim_max[2] - protDim_min[2]) ** 2)/2
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
	frames = str(160)
	cmd.mset('1 x' + frames)
	cmd.orient('all')
	
	cmd.set_view (
		'-0.346189618,   -0.600477457,   -0.720818758,\
     0.038613044,    0.758557022,   -0.650460541,\
     0.937369645,   -0.253015876,   -0.239418939,\
     0.000000000,    0.000000000, -1473.033813477,\
   -73.633697510,  -31.334754944,    1.869384766,\
  1161.351074219, 1784.716552734,  -20.000000000 ')

			
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
		
		## list all chain and ligand pairs
		chainAndLigand = {}
		ligandAndChain = {}
		
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

				## get pair of ligand and respective chain
				if cmd.count_atoms(selection) > 0:
					cmd.extract(ligandname, selection)
					## calculate COM of ligand
					ligandsCOMS[ligandname] = calc_COM(ligandname)
				
					## color binding site
					binding = 'byres (' + chainname + ' nto. 3 of organic)'
					cmd.select('inter', binding)
					cmd.color('red', 'inter')
					cmd.delete('inter')
					
					## group chains and respective ligands to translate them together at first
					cmd.group('CL_' + chainname + "_" + ligandname, chainname + " and " + ligandname)
					chainAndLigand[chainname] = 'CL_' + chainname + "_" + ligandname
					ligandAndChain[ligandname] = chainname
		
		
		cmd.frame(f+10)
		## store chain objects for movie
		for cname in cNames:
			cmd.mview('store', object=cname)
		for ligand in ligandsCOMS.keys():
			cmd.mview('store', object=ligand)
		
		f = frameNum
		## translate chains iteratively 
		for chainname in cNames:
			
			## COM of chain
			chainXYZ = chainsCOMS[chainname]

			## only translate chain if there are multiple chains
			if len(chains) > 1:	
			
				## translate chains and ligands
				for cname in cNames:
					if cname != chainname:
						while isColliding(chainname, cname):
							for cname in cNames:
								cXYZ = chainsCOMS[cname]
								if cname in chainAndLigand.keys():
									translate_selection(
										complexXYZ, cXYZ, 
										chainAndLigand[cname], transFac, f)
									continue
								translate_selection(
										complexXYZ, cXYZ, cname, transFac, f)
		
		## ungroup chains and ligands to translate ligands seperately from chains
		for cl in chainAndLigand.keys():
			cmd.ungroup(cl)
			cmd.ungroup(cl+'_org')
			cmd.delete(chainAndLigand[cl])

		## store objects for movie
		cmd.frame(f)
		for ligand in ligandsCOMS.keys():
			cmd.mview('store', object=ligand)
		for cname in cNames:
			cmd.mview('store', object=cname)
			
		f = f + frameNum
		cmd.frame(f)
		
		## translate ligands
		for ligand in ligandAndChain.keys():
			cXYZ = chainsCOMS[ligandAndChain[ligand]]
			ligandXYZ = ligandsCOMS[ligand]
			while isColliding(ligand, ligandAndChain[ligand]):
				translate_selection(cXYZ, ligandXYZ, ligand, transFac/2, f)
		for cname in cNames:
			cmd.mview('store', object=cname)
		f = f + frameNum

		## delete state object
		cmd.delete(pdbName + '_'+ s)	
		

	print time.clock() - start_time, 'seconds'
	
cmd.extend('explosion', explosion)
cmd.extend('isColliding', isColliding)