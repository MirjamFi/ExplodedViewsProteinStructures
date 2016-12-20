# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_movie.py

'''
explosion creates a movie of the exploded view of a molecule. 
If there are multiple objects given for explosion, first they get spatially 
separated and then explode individually one after the other. The order of 
explosion is the order of given list of objects.
There are two types of explosion direction:
- 'com' (default): 	the centers of mass (com) of the chains of the object to be 
					tranlated and the object are calulated and the single chains 
					are translated along a vector through the chain's com and 
					object's com.
- 'canonical':		the dimensions of a box around the object are used to select 
					the two longest edges and so the axes to translate along in	
					a consistent distance
					
If only a part of the object shall be translated the object can be given as 
complex to make sure the part is not translated into the object.

DEPENDENCIES:
get_colors.py (https://pymolwiki.org/index.php/Get_colors) and 
center_of_mass.py (https://pymolwiki.org/index.php/Center_of_mass)
in the modules of PyMOL.

'''
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time
from operator import itemgetter

def initialize_movie(selected = None, frames = "100"):
	''' DESCRIPTION: 
		initial setup of a movie'''
	cmd.set('matrix_mode', 1)
	cmd.set('movie_panel', 1)
	cmd.set('scene_buttons', 1)
	cmd.set('cache_frames', 1)
	cmd.config_mouse('three_button_motions', 1)
	# cmd.set('movie_panel', 0)	## hide movie panel
	cmd.set('movie_panel_row_height', 1)
	cmd.set('movie_fps', 5)
	
	cmd.mset('1 x' + frames)
	if selected:
		cmd.orient(selected)
		
def get_ligands():
	''' DESCRIPTION:
		get all ligands (organic) from complex, return them as a set'''
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

	return set(stored.ligands)
	
def calc_COM(transname):
	'''DESCRIPTION:
		calculate CenterOfMass of an object 
	'''
	cenma.com(transname, state = 1)
	cmd.zoom(transname + '_COM')
	pos = cmd.get_position(transname + '_COM')
	cmd.delete(transname + '_COM')
	cmd.zoom('all', complete = 1)
	cmd.mview('store')

	return pos

def calcTransFac(sele):
	'''DESCRIPTION 
		calculate translation factor according to size of selection
	'''
	protDim_min, protDim_max = cmd.get_extent(sele)
	transFac = math.sqrt((protDim_max[0] - protDim_min[0]) ** 2 +
					 (protDim_max[1] - protDim_min[1]) ** 2 +
					 (protDim_max[2] - protDim_min[2]) ** 2)/2
	return transFac

def transAxes(selected):
	'''DESCRIPTION:
		get dimensions of structure to decide along which axes to translate
	'''
	box1_min, box1_max = cmd.get_extent(selected)

	x_axis = abs(box1_min[0] - box1_max[0])
	y_axis = abs(box1_min[1] - box1_max[1])
	z_axis = abs(box1_min[2] - box1_max[2])

	axes = {"x_axis" : x_axis, "y_axis":y_axis, "z_axis":z_axis}
	axes = sorted(axes.items(), key=itemgetter(1), reverse = True)

	transFac = calcTransFac(selected)
	transVec = [0,0,0]
	if axes[0][0] == "x_axis" or axes[1][0] == "x_axis":
		transVec[0] = transFac
	if axes[0][0] == "y_axis" or axes[1][0] == "y_axis":
		transVec[1] = transFac	
	if axes[0][0] == "z_axis" or axes[1][0] == "z_axis":
		transVec[2] = transFac
		
	return transVec
	
def create_objects(chains, selected, storedLigands, chainAndLigand = None, typeOfExplosion = 'com'):

	## names of chains
	cNames = []
	
	## chains and their labels
	chainAndLabel = {}
	
	## chains and according ligands
	if not chainAndLigand:
		chainAndLigand = {}
	
	## ligands and according chains
	ligandAndChain={}
	
	## COMS of ligands/chains
	ligandsCOMS = {}
	chainsCOMS = {}
	
	for c in chains:	
		## color each chain individually
		cmd.color(get_colors.get_random_color(), selected + \
					'& chain ' + c) 
					
		## create object for chain with name
		chainname = selected + '_chain' + c
		cNames.append(chainname)
		cmd.extract(chainname, selected + ' & chain '+c)
		if not storedLigands:
			cmd.show('cartoon',chainname)
		else:
			cmd.show('surface', chainname)
		
		## label chains
		chainAndLabel = label_obj(chainname, chainAndLabel)
		
		## store object for movie
		f = 1
		cmd.frame(f)
		cmd.mview('store', object=chainname)
		cmd.mview('store', object="_label" + chainname)
		
		## calculate and save COM of chain
		if typeOfExplosion == 'com':
			chainsCOMS[chainname] = calc_COM(chainname)
		
		## create objects for ligands
		if storedLigands:
			for l in storedLigands:
			
				## select ligand on chain	
				selection = chainname + ' & resn ' + l
				
				if typeOfExplosion == 'com':
					get_ligand_chain_pair(l, selection, chainname, 
										chainAndLigand, ligandAndChain, 
										ligandsCOMS, COM = True)
				else: 
					get_ligand_chain_pair(l, selection, chainname, 
										chainAndLigand, ligandAndChain)
		
		## if there is no ligand on chain, just keep chain and its label
		if not chainAndLigand[chainname]:
			if chainname in chainAndLigand: del chainAndLigand[chainname]
			cmd.group(chainname + "_", chainname + " and " + "_label" \
						+ chainname)
	f = f + 20
	cmd.frame(f)
	## store chain objects for movie 
	store_view(group = True, all = True)
	f = f + 1
	
	## hide labels
	cmd.hide('labels')
	cmd.scene('off', 'store')
	cmd.frame(f)
	cmd.mview('store', scene='off')
	
	

		
	if typeOfExplosion == 'com':
		return cNames, chainAndLabel, chainAndLigand, ligandAndChain, ligandsCOMS, chainsCOMS, f
	else:
		return cNames, chainAndLabel, chainAndLigand, ligandAndChain, f

def get_ligand_chain_pair(l, selection, chainname, chainAndLigand, 
							ligandAndChain, ligandsCOMS = None, COM = False):
	'''DESCRIPTION:
		get pair of ligand and respective chain''' 
	
	## create object for ligand
	if cmd.count_atoms(selection) > 0:
		ligandname = chainname + '_' + l
		cmd.extract(ligandname, selection)
		cmd.show('spheres', ligandname)
		
		## label ligands			
		label_obj(ligandname)
		
		if COM:
			## calculate COM of ligand
			ligandsCOMS[ligandname] = calc_COM(ligandname)
		
		## color binding site
		color_binding(chainname, ligandname)
		
		## group chains and respective ligands to translate  
		## them together at first (including their labels)
		if COM:
			cmd.group(chainname + "_" + l + "_", 
						chainname + " " + ligandname + " " + "_label" + \
						chainname + " " + "_label" + ligandname)
			chainAndLigand[chainname].append(chainname + "_" + l + "_")
		else:
			## group chains and respective ligands to translate  
			## them together at first (including their labels)
			cmd.group(chainname + "_", 
						chainname + " " + ligandname + " " + "_label" + \
						chainname + " " + "_label" + ligandname)
			chainAndLigand[chainname].append(ligandname)
			
		ligandAndChain[ligandname] = chainname  
		
def label_obj(chainname, chainAndLabel = None):
	'''DESCRIPTION:
		create an label for given chain (if chainAndLabel set) or ligand  
	'''
	## create pseudoatom at corner of boundingbox of object to be labeled
	dim = cmd.get_extent(chainname)
	cmd.pseudoatom("_label" + chainname, pos=dim[0])
	
	## label ligand
	if not chainAndLabel:
		ligName = chainname.split("_")[-1]
		cmd.label("_label" + chainname, "'%s'" %ligName)
		
	## label chain
	else:
		cmd.label("_label" + chainname, "'%s'" %chainname)
		
	## save label of chain in dictonary
	if chainAndLabel:
		chainAndLabel[chainname] = "_label" + chainname
		
	## store ligand label in movie
	if not chainAndLabel:
		cmd.mview('store', object = "_label" + chainname)
		
	## hide pseudoatom represantation
	cmd.hide('nonbonded', "_label" + chainname)
	
	## return dictionary of chain and its label
	if chainAndLabel:
		return chainAndLabel
	
def color_binding(chainname, ligandname, res = False):
	''' DESCRIPTION:
		color binding site '''
	binding = 'byres (' + chainname + ' nto. 3.6 of ' + \
				ligandname +')'
	cmd.select('_inter', binding)
	cmd.color('red', '_inter')
	if not res:
		cmd.delete('_inter')
		
def create_complex(chains, obj):
	'''DESCRIPTION:
		for list of chains create object from chains belonging to obj
	'''
	sel = ''
	for ch in chains[obj]:
		sel = sel + obj+'_chain' + ch + ' '
	cmd.select('_sel' + obj, sel)
	
def isColliding(sel1, sel2):
	''' DESCRIPTION:
		calculate if two bounding boxes are colliding (boxes are axis-aligned)
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
	
	## test if any of the vertices of box 1 are in any of the faces of box 2 
	## and vice versa
	isColliding = 	(x_min1 <= x_max2 and x_max1 >= x_min2) and \
					(y_min1 <= y_max2 and y_max1 >= y_min2) and \
					(z_min1 <= z_max2 and z_max1 >= z_min2)

	return isColliding
	
def translate_selection(originXYZ, transXYZ, transname, factor = 1, f = 1, group = None):
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
	if group:
		cmd.mview('store', object=group)
		cmd.mview('interpolate', object=group)
	else:
		cmd.mview('store', object=transname)
		cmd.mview('interpolate', object=transname)
	
	store_view(group = True, all = True)
	
def com_translation(cname, chainAndLigand, complexXYZ, chainXYZ, transFac, f):
	'''DESCRIPTION:
		translate chain via COM
	'''		
	## if chain contains ligand, translate ligand also
	if cname in chainAndLigand.keys():
		for l in chainAndLigand[cname]:
			translate_selection(complexXYZ, chainXYZ, l, transFac, f, l)
	## only chain has to be translated
	else:
		translate_selection(complexXYZ, chainXYZ, cname + "_", transFac,f, 
								cname + "_")

def canonical_tranlation(ch, i, transVec, chainAndLigand):
	'''DESCRIPTION:
		translate chain ch canonical
	'''
	cmd.translate( [x * i for x in transVec] , object=ch)
	cmd.translate([x * i for x in transVec], 
					object="_label" + ch)
	if ch in chainAndLigand.keys():
		for l in chainAndLigand[ch]:
			cmd.translate([x * i for x in transVec], object=l)
			cmd.translate([x * i for x in transVec], object="_label" + l)
	cmd.mview('store', object = ch+"_")
	cmd.mview('interpolate', object = ch+"_")
		
def store_view(obj = None, group = False, all = True):
	''' DESCRIPTION:
		store selected view '''
		
	## store all objects
	if group:
		for group in cmd.get_names("objects"):
				cmd.mview('store', object = group)
	
	## zoom out on everything
	if all:
		cmd.zoom('all', complete=1)
	
	## store single object
	if obj:
		cmd.mview('store', object = obj)
		
	## store camera position
	cmd.mview('store')
			
def store_res_view(selected, f, resis = None, orient = False, all = False):
	'''DESCRIPTION:
		store movie view for single residue/ligand
	'''
	cmd.frame(f)
	## orient on residue
	if orient:
		cmd.orient('_inter %s' %selected)
	
	## store view
	if resis:
		store_view(resis + "_", all)
	if not resis:
		if all:
			store_view(selected, all)
	
def com_explosion(selected, cNames = None, chainAndLabel = None, 
					chainAndLigand = None, ligandAndChain = None, 
					chainsCOMS = None, complex = None, com = None, 
					storedLigands = None, ligandsCOMS = None, transFac = None,
					frame = 1):
	''' DESCRIPTION:
		create a movie for an exploded view of given protein. Explosion along
		vector of complex' COM and selection's COM.
	'''
	start_time = time.clock()
	
	f = frame 
			
	''' Explosion: Create an object of every chain and for related ligands. 
				Calulate COM for complex and according chains and 
				translate chain along vector between COMs of complex and chain. 
				For every chain do the same with its ligands.
	'''
		
	## COM of complex
	## if only some chains are selected and source complex is given,
	## select source complex for COM calculation
	if com:
		complexXYZ = com
	if not com:
		complexXYZ = calc_COM(selected)
	
	''' translate chains '''
	for chainname in cNames:
		
		## COM of chain
		chainXYZ = chainsCOMS[chainname]

		## only translate chain if there are multiple chains
		if len(cNames) > 1:	
			## translate chains and ligands
			for cname in cNames:
			
				## check if current chain is colliding with any other chain
				if cname != chainname:
					while isColliding(chainname, cname):
					
						## if collision, translate all chains in selection
						for cname in cNames:
							cXYZ = chainsCOMS[cname]
							com_translation(cname, chainAndLigand, complexXYZ, cXYZ, transFac, f)
									
			store_view(group = True, all = True)
									
		## if only one chain is selected, translate it and its ligand
		else:
			cmd.frame(f)
			cmd.zoom(chainname)
			store_view(chainname + '_')
			
			f = f + 30
			cmd.frame(f)
			
			## if source complex is known use its COM to translate
			if complex:
				while isColliding(chainname, complex):
					com_translation(chainname, chainAndLigand, complexXYZ, chainXYZ, transFac, f)
			else:
				if chainname in chainAndLigand.keys():
						cmd.translate([transFac, transFac, transFac], 
										object = chainAndLigand[chainname])
				else:
					cmd.translate([transFac, transFac, transFac], 
										object = chainname + "_")
										
			store_view(chainname + '_', all = True)
					
			## store objects for movie
			f = f + 30
			cmd.frame(f)
			store_view(chainname + '_', all = True)

	## store objects for movie
	f = frame + 30
	if not ligandAndChain:
		## show labels
		cmd.show('labels')
		cmd.scene('on', 'store')
		cmd.frame(f)
		cmd.mview('store', scene='on')
		store_view(group=True, all = True)

	else:
		cmd.frame(f)
		store_view(group=True, all = True)

	''' translate ligands '''
	if ligandAndChain:
		f = f + 30
		cmd.frame(f)

		for ligand in ligandAndChain.keys():
			cXYZ = chainsCOMS[ligandAndChain[ligand]]
			ligandXYZ = ligandsCOMS[ligand]
			
			if complex:
				condition = isColliding(ligand, complex) \
							or isColliding(ligand, ligandAndChain[ligand])
			else:
				condition = isColliding(ligand, ligandAndChain[ligand])
			if condition:
				transFac = transFac
				translate_selection(cXYZ, ligandXYZ, ligand, 
									transFac, f, 
									chainAndLigand[ligandAndChain[ligand]])
				translate_selection(cXYZ, ligandXYZ, 
										"_label" + ligand, 
										transFac, f, 
										chainAndLigand[ligandAndChain[ligand]])
				if complex:
					condition = isColliding(ligand, complex) \
							or isColliding(ligand, ligandAndChain[ligand])
				else:
					condition = isColliding(ligand, ligandAndChain[ligand])

		f = f + 10
		## show labels
		cmd.show('labels')
		cmd.scene('on', 'store')
		cmd.frame(f)
		cmd.mview('store', scene='on')
		cmd.mview('reinterpolate')

		f = f + 20
		cmd.frame(f)
		store_view(group=True, all = True)
	print 'Explosion of', selected, time.clock() - start_time, 'seconds'
	return f

def canonical_explosion(selected, cNames = None, chainAndLabel = None, 
					chainAndLigand = None, ligandAndChain = None, complex = None,
					storedLigands = None, transVec = None, frame = 1):
	'''DESCRIPTION: translate a selection in canonical direction and create
		a movie from it'''
	start_time = time.clock()
	cmd.orient(selected)
	
	f = frame
	''' translate chains'''
	i = 1
	cmd.frame(f)
	if complex:
		for chain in cNames:
			while isColliding(chain, complex):
				for ch in cNames:
					canonical_tranlation(ch, i, transVec, chainAndLigand)
					i += 1
	else:
		for chain in cNames:
			for c in cNames:
				if c != chain:
					while isColliding(chain, c):
						for ch in cNames[1:]:
							canonical_tranlation(ch, i, transVec, chainAndLigand)
							i += 1

				
	store_view(group = True, all = True)
	
	f = f + 30
	if not ligandAndChain:
		## show labels
		cmd.show('labels')
		cmd.scene('on', 'store')
		cmd.frame(f)
		cmd.mview('store', scene='on')
		store_view(group=True, all = True)

	else:
		cmd.frame(f)
		store_view(group=True, all = True)
		
	''' translate ligands'''
	if ligandAndChain:
		f = f + 30
		cmd.frame(f)
		for ligand in ligandAndChain.keys():
			while isColliding(ligand, ligandAndChain[ligand]):
				cmd.translate([x * 1/4  for x in transVec], object=ligand)
				cmd.translate([x * 1/4 for x in transVec], object="_label" + ligand)
			store_view(group = True, all = True)
			
		f = f + 30
		## show labels
		cmd.show('labels')
		cmd.scene('on', 'store')
		cmd.frame(f)
		cmd.mview('store', scene='on')
		cmd.mview('reinterpolate')
		store_view(group=True, all = True)
					
	print 'Explosion of', selected, time.clock() - start_time, 'seconds'
	return f
	
def explosion(selected = [], typeOfExplosion = 'com', complex = None):
	'''DESCRIPTION:
		perform an explosion of selected object(s) given in a list and create	
		a movie. If two objects are given they will be separateted and then 
		exploded individually.
	'''

	if not typeOfExplosion in ['com', 'canonical']:
		sys.exit("Specify explosion: com (default) or canonical.") 
		
	if not selected:
		sys.exit("Please give a list of selected objects (or just one) to translate.") 
		
	'''setup of selected sturcture'''
	## case sensitive for chain ids
	cmd.set('ignore_case', 'off') 
	# cmd.set('ray_trace_mode', 1)
	# cmd.set('ray_trace_frames', 'on')
	
	## remove solvent 
	cmd.remove('solvent')
	cmd.hide('all')
	
	if complex:
		cmd.show('surface', complex)
		
	'''preparation and translation'''
	if len(selected) >= 1:
		## initialize movie
		start_time = time.clock()
		if len(selected) > 1:
			initialize_movie(frames = str(130*len(selected)+100))
		else:
			initialize_movie(frames = str(150))
		
		## get ligands
		storedLigands = get_ligands()

		## lists for further computations
		s = ''
		chains = {}
		cNames = []
		chainAndLabel = {}
		chainAndLigand = {}
		ligandAndChain = {}
		chainsAndObj = {}
		
		if typeOfExplosion == 'com':
			coms = {}
			ligandsCOMS = {}
			chainsCOMS = {}
		
		for obj in selected:
			if typeOfExplosion == 'com':
				## calculate com of obj
				coms[obj] = calc_COM(obj)
			s = s + ' '+ obj
			
		if not complex:	
			cmd.create('_all_obj', s)
			
			if typeOfExplosion == 'com':
				## calc com of complex from objects
				complexXYZ = calc_COM('_all_obj')
				trans = calcTransFac('_all_obj')
			else:
				## calculate translation vector
				transVec = transAxes('_all_obj') 
				print transVec
			cmd.frame(1)
			cmd.orient('_all_obj')
			cmd.mview('store')
			cmd.delete('_all_obj')
				
		else:
			if typeOfExplosion == 'com':
				## calc com of complex from objects
				complexXYZ = calc_COM(complex)
				trans = calcTransFac(complex)
			else:
				## calculate translation vector
				transVec = transAxes(complex) 
				print transVec
			
		for obj in selected:
			##get chains of complex
			chains[obj] = cmd.get_chains(obj)
			for ch in chains[obj]:
				chainsAndObj[obj+ '_chain' + ch] = obj
				chainAndLigand[obj+ '_chain' + ch] = []
				
			if typeOfExplosion == 'com':
				## create objects for all chains and ligands
				cNames_new, chainAndLabel_new, chainAndLigand_new, ligandAndChain_new, ligandsCOMS_new, chainsCOMS_new, f= \
								create_objects(chains[obj], obj, storedLigands, chainAndLigand)
			else:
				cNames_new, chainAndLabel_new, chainAndLigand_new, ligandAndChain_new, f= \
								create_objects(chains[obj], obj, storedLigands, chainAndLigand, typeOfExplosion = 'canonical')
		
			cNames = cNames + cNames_new
			if chainAndLabel_new:
				chainAndLabel.update(chainAndLabel_new)
			if chainAndLigand_new:
				chainAndLigand.update(chainAndLigand_new)
			if ligandAndChain_new:
				ligandAndChain.update(ligandAndChain_new)
				
			if typeOfExplosion == 'com':
				if ligandsCOMS_new:
					ligandsCOMS.update(ligandsCOMS_new)
				if chainsCOMS_new:
					chainsCOMS.update(chainsCOMS_new)

		f = f + 30		
		if len(selected) > 1:
			if typeOfExplosion == 'canonical':
				i = 1
			cmd.frame(f)
			''' separate objects '''
			for obj in chains.keys():
				create_complex(chains, obj)
					
				for obj2 in chains.keys():
					if obj != obj2:
						create_complex(chains, obj2)
						
						## check if obj is colliding with any other obj
						while isColliding('(_sel%s)'%obj, '(_sel%s)'%obj2):
							if typeOfExplosion == 'com':
								## if collision, translate all chains in selection
								for cname in cNames:
									com_translation(cname, chainAndLigand, complexXYZ, coms[chainsAndObj[cname]], trans, f)
								
							else:
								for o in chains.keys():
									for ch in chains[o]:
										ch = obj+ '_chain' + ch
										canonical_tranlation(ch, i, transVec, chainAndLigand)
									i = i + 1
							create_complex(chains, obj)
							create_complex(chains, obj2)
										
				store_view(group = True, all = True)
				
			f = f+ 30
			cmd.frame(f)
			store_view(group = True, all = True)
				
			cmd.delete('(_sel%s)'%obj)
			cmd.delete('(_sel%s)'%obj2)
			
			print 'Preparation:', time.clock() - start_time, 'seconds'
		
		''' translate objects individually '''
		for obj in chains.keys():
			if len(selected) > 1:
				f = f + 30
			sel1=''
			for c in chains[obj]:
				sel1 = sel1 + obj+'_chain' + c + ' '
			cmd.select('_'+obj, sel1)
			if typeOfExplosion == 'com':
				trans = calcTransFac('_'+obj)
			else:
				trans = transAxes('_'+obj)
			
			## get object's chains and ligands
			objChains = [chain for chain, o in chainsAndObj.items() if o == obj]
			ligandAndChain_obj= {}
			for ch in objChains:
				if ch in ligandAndChain.values():
					ligandAndChain_obj.update({ligand: c for ligand, c in ligandAndChain.items() if c == ch})
			
			## translation of object
			if typeOfExplosion == 'com':
				if not complex:
					f = com_explosion('_'+obj, cNames = objChains,
								chainAndLabel = chainAndLabel, 
								chainAndLigand = chainAndLigand, 
								ligandAndChain = ligandAndChain_obj, 
								chainsCOMS= chainsCOMS, com = coms[obj], 
								storedLigands = storedLigands, 
								transFac = trans, ligandsCOMS = ligandsCOMS,
								frame = f)
				else:
					f = com_explosion('_'+obj, cNames = objChains,
								chainAndLabel = chainAndLabel, 
								chainAndLigand = chainAndLigand, 
								ligandAndChain = ligandAndChain_obj, 
								chainsCOMS= chainsCOMS, complex = complex,
								com = complexXYZ, storedLigands = storedLigands, 
								transFac = trans, ligandsCOMS = ligandsCOMS,
								frame = f)
			else:
				if complex:
					f = canonical_explosion('_'+obj, cNames = objChains,
							chainAndLabel = chainAndLabel, 
							chainAndLigand = chainAndLigand, 
							ligandAndChain = ligandAndChain_obj, complex = complex,
							storedLigands = storedLigands, transVec = trans,
							frame = f)
				else:
					f = canonical_explosion('_'+obj, cNames = objChains,
							chainAndLabel = chainAndLabel, 
							chainAndLigand = chainAndLigand, 
							ligandAndChain = ligandAndChain_obj, 
							storedLigands = storedLigands, transVec = trans,
							frame = f)
				
				
			cmd.delete('_'+obj)
			
def relabel(selected, newLabel="new label"):
	'''DESCRIPTION:
		rename an selected object and its label by a new label
	'''
	obj = cmd.get_names("objects")
	x = ''
	
	## relabel label
	cmd.label("_label" + selected , "'%s'" %newLabel )
	
	## rename all objects of selection
	for o in obj:
		if selected in o:
			x = o
			x = x.replace(selected, newLabel)
			cmd.set_name(o, x)

'''TODO: video muss zweimal gespeichert werden (bei erstem mal noch wackeln in 
			letztem gespeicherten frame) '''
def reorient_explosion(frame=1):
	'''DESCRIPTION:
		if called, the orientation of the movie is set to the orientation of 
		given frame from frame till end. '''
	x = cmd.get_view(quiet=1)
	frames = cmd.count_frames()
	for f in range(int(frame),frames+1):
		cmd.frame(f)
		cmd.set_view(x)
		cmd.zoom('all')
		cmd.mview('store')
	
def renew_representation(selection, representation):
	'''DESCRIPTION:
		show selection in given representation
	'''
	cmd.scene(selection+'_'+representation, 'store')
	for f in range(1, cmd.count_frames()):
		cmd.frame(f)
		cmd.mview('store', scene=selection+'_'+representation)
		cmd.show_as(representation, selection)
		cmd.mview('reinterpolate')
		
cmd.extend('reorient_explosion', reorient_explosion)
cmd.extend('explosion', explosion)
cmd.extend('relabel', relabel)
cmd.extend('renew_representation', renew_representation)
