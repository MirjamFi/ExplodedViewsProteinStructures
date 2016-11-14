# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_movie.py

# ''' load complex from PDB'''
# cmd.reinitialize()
# selected = '3oaa'
# pdbFile = selected + '.pdb1'
# cmd.fetch(selected, type='pdb1')
# util.cbc(selection= selected)
# cmd.remove('solvent')
# cmd.load('C:\Users\Figaschewski\Dropbox\Masterarbeit\TPC2_M484L_wActivator_56ns\TPC2_M484L_wActivator_56ns.pdb', 'test')

'''
	Create an movie of an exploded view for a given protein structure (PDB) or 
	parts of it.
	If you use explosion() it will explode parts alon a vector of the center of
	mass of the original complex and the COM of the selection. Here, xou can 
	select single or multiple chains, a residue or ligand to be exploded, 
	just extract them into one object from source structure. If chains contain 
	ligands, those 	will also be exploded.
	To explod specific chains/residue/ligand it is useful to also give the 
	complex for explosion directon, espacially for single chain, so the chain is 
	not	translated into the complex.
	If you use canonical_explosion() the parts of the selection will be 
	translated in the same direction. You can select the whole complex or single/
	multiple chains.
	Exploded parts will be labeled.
'''
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time
from operator import itemgetter

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

def translate_selection(originXYZ, transXYZ, transname, factor, f, group):
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
	cmd.mview('store', object=group)
	cmd.mview('interpolate', object=group)
	
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
		
def isColliding(sel1, sel2):
	''' DESCRIPTION:
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
	
	## test if any of the vertices of box 1 are in any of the faces of box 2 
	## and vice versa
	isColliding = 	(x_min1 <= x_max2 and x_max1 >= x_min2) and \
					(y_min1 <= y_max2 and y_max1 >= y_min2) and \
					(z_min1 <= z_max2 and z_max1 >= z_min2)
	

	return isColliding

def calcTransFac(sele):
	'''DESCRIPTION 
		calculate translation factor according to size of selection
	'''
	protDim_min, protDim_max = cmd.get_extent(sele)
	transFac = math.sqrt((protDim_max[0] - protDim_min[0]) ** 2 +
					 (protDim_max[1] - protDim_min[1]) ** 2 +
					 (protDim_max[2] - protDim_min[2]) ** 2)/2
	return transFac

def label_obj(chainname, chainAndLabel = None):
	'''DESCRIPTION:
		create an label for given chain (if chainAndLabel set) or ligand  '''
	dim = cmd.get_extent(chainname)
	if not chainAndLabel:
		ligName = chainname.split("_")[-1]
	cmd.pseudoatom("_tmpPoint" + chainname, pos=dim[0])
	# cmd.pseudoatom("_tmpPoint" + chainname, pos=dim[1])
	if not chainAndLabel:
		cmd.label("_tmpPoint" + chainname, "'%s'" %ligName)
	else:
		cmd.label("_tmpPoint" + chainname, "'%s'" %chainname)
	if chainAndLabel:
		chainAndLabel[chainname] = "_tmpPoint" + chainname
	if not chainAndLabel:
		cmd.mview('store', object = "_tmpPoint" + chainname)
	cmd.hide('nonbonded', "_tmpPoint" + chainname)
	if chainAndLabel:
		return chainAndLabel
	
def initialize_movie(selected):
	''' DESCRIPTION: 
		initial setup of a movie'''
	cmd.set('matrix_mode', 1)
	cmd.set('movie_panel', 1)
	cmd.set('scene_buttons', 1)
	cmd.set('cache_frames', 1)
	cmd.config_mouse('three_button_motions', 1)
	# cmd.set('movie_panel', 0)	## hide movie panel
	cmd.set('movie_panel_row_height', 3)
	
	frameNum = 50
	frames = str(160)
	cmd.mset('1 x' + frames)
	cmd.orient(selected)
	
def color_binding(chainname, ligandname, res = False):
	''' DESCRIPTION:
		color binding site '''
	binding = 'byres (' + chainname + ' nto. 3 of ' + \
				ligandname +')'
	cmd.select('_inter', binding)
	cmd.color('red', '_inter')
	if not res:
		cmd.delete('_inter')
	
def store_view(obj = None, group = False, all = True):
	''' DESCRIPTION:
		store selected view '''
	if group:
		for group in cmd.get_names("objects"):
					cmd.mview('store', object = group)
	if all:
		cmd.zoom('all', complete=1)
	if obj:
		cmd.mview('store', object = obj)
	cmd.mview('store')

def get_ligand_chain_pair(l, selection, chainname, chainAndLigand, 
							ligandAndChain, ligandsCOMS = None, COM = False):
	'''DESCRIPTION:
		get pair of ligand and respective chain''' 
	if cmd.count_atoms(selection) > 0:
		ligandname = chainname + '_' + l
		cmd.extract(ligandname, selection)							
		
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
							chainname + " " + ligandname \
							+ " " + "_tmpPoint" + chainname \
							+ " " + "_tmpPoint" + ligandname)
			chainAndLigand[chainname] = chainname + \
											"_" + l + "_"
		else:
			## group chains and respective ligands to translate  
			## them together at first (including their labels)
			cmd.group(chainname + "_", 
							chainname + " " + ligandname \
							+ " " + "_tmpPoint" + chainname \
							+ " " + "_tmpPoint" + ligandname)
			chainAndLigand[chainname] = ligandname
			
		ligandAndChain[ligandname] = chainname  
			
def store_res_view(selected, f, resis = None, orient = False, all = False):
	cmd.frame(f)
	if orient:
		cmd.orient('_inter %s' %selected)
	if resis:
		store_view(resis + "_", all)
	if not resis:
		if all:
			store_view(selected, all)


def com_explosion(selected, complex = None):
	''' DESCRIPTION:
		create a movie for an exploded view of given protein. Explosion along
		vector of complex' COM and selection's COM.
	'''
	'''setup of selected sturcture'''
	start_time = time.clock()
	
	## case sensitive for chain ids
	cmd.set('ignore_case', 'off') 
	
	cmd.remove('solvent')
	cmd.show('spheres', 'organic')
	util.cbc(selection= selected)
	
	## calculate translation factor which is the size of the selection
	transFac = calcTransFac(selected)
	
	##get chains of complex
	chains = cmd.get_chains(selected)
	
	# for c in chains:
		# if c == "":
			# chains.remove(c)

	## get ligands
	storedLigands = get_ligands() 

	## initialize movie
	initialize_movie(selected)
			
	''' Explosion: Create an object of every chain and for related ligands. 
				Calulate COM for complex and according chains and 
				translate chain along vector between COMs of complex and chain. 
				For every chain do the same with its ligands.
	'''
		
	## COM of complex
	## if only some chains are selected and source complex is given,
	## select source complex for COM calculation
	if complex:
		complexXYZ = calc_COM(complex)
	else:
		complexXYZ = calc_COM(selected)

	## calculate COM for single chains in state
	cNames = []
	chainsCOMS = {}
	
	## COM for ligands
	ligandsCOMS = {}
	f = 1
	
	## list all chain and ligand pairs
	chainAndLigand = {}
	ligandAndChain = {}
	chainAndLabel = {}
	
	''' chain or complex shall be translated '''
	stored.res=[]
	## test if not only a single residue/ligand is selected
	cmd.iterate('(' + selected + ')',"stored.res.append(resn + resi)")
	if len(set(stored.res)) > 1:
		cmd.hide("(" + selected + ")")
		for c in chains:	
			## color each chain individually
			if cmd.count_states() > 1:
				cmd.color(get_colors.get_random_color(), selected + \
						'& chain ' + c) 
						
			## create object for chain with name
			chainname = 'chain' + c
			cNames.append(chainname)
			cmd.extract(chainname, selected+ ' & chain '+c)
			
			## label chains
			chainAndLabel = label_obj(chainname, chainAndLabel)
			
			## store object for movie
			cmd.frame(f)
			cmd.show('sticks', chainname)
			cmd.mview('store', object=chainname)
			cmd.mview('store', object="_tmpPoint" + chainname)
			
			## calculate and save COM of chain
			chainsCOMS[chainname] = calc_COM(chainname)	
			
			## create objects for ligands
			if stored.ligands:
				for l in storedLigands:
				
					## select ligand on chain	
					selection = chainname + ' & resn ' + l
					
					get_ligand_chain_pair(l, selection, chainname, 
											chainAndLigand, ligandAndChain, 
											ligandsCOMS, COM = True)
			
			## if there is no ligand on chain, just keep chain and its label
			if chainname not in chainAndLigand.keys():
				cmd.group(chainname + "_", chainname + " and " + "_tmpPoint" 
							+ chainname)
		f = f + 20
		cmd.frame(f)
		## store chain objects for movie 
		store_view(group = True, all = True)
		
		f = f + 30
		''' translate chains iteratively '''
		for chainname in cNames:
			
			## COM of chain
			chainXYZ = chainsCOMS[chainname]

			## only translate chain if there are multiple chains
			if len(chains) > 1:	
				## translate chains and ligands
				for cname in cNames:
				
					## check if current chain is colliding with any other chain
					if cname != chainname:
						while isColliding(chainname, cname):
						
							## if collision, translate all chains in selection
							for cname in cNames:
								cXYZ = chainsCOMS[cname]
								
								## if chain contains ligand, translate ligand also
								if cname in chainAndLigand.keys():
									translate_selection(complexXYZ, cXYZ, 
											chainAndLigand[cname], transFac, f, 
											chainAndLigand[cname])
									continue
								## only chain has to be translated
								translate_selection(complexXYZ, cXYZ, 
										cname + "_", transFac,f, cname + "_")
										
				store_view(all = True)
										
			## if only one chain is selected, translate it and its ligand
			else:
				cmd.frame(f)
				cmd.orient(chainname)
				store_view(chainname + '_', all = True)

				f = f + 10
				cmd.frame(f)
				store_view(chainname + '_', all = True)
				
				f = f + 30
				cmd.frame(f)
				
				## if source complex is known use its COM to translate
				if complex:
					while isColliding(chainname, complex):
						if chainname in chainAndLigand.keys():
							translate_selection(complexXYZ, chainXYZ, 
									chainAndLigand[chainname], transFac, f, 
									chainAndLigand[chainname])
						else:
							translate_selection(complexXYZ, chainXYZ, 
									chainname + "_", transFac, f, 
									chainname + "_")
				else:
					if chainname in chainAndLigand.keys():
							cmd.translate([transFac, transFac, transFac], 
											object = chainAndLigand[chainname])
					else:
						cmd.translate([transFac, transFac, transFac], 
											object = chainname + "_")
											
				cmd.orient(chainname)
				store_view(chainname + '_', all = True)
						
		
		## store objects for movie
		f = f + 30
		cmd.frame(f)
		store_view(group=True, all = True)
		
		''' translate ligands '''
		if storedLigands:
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
					translate_selection(cXYZ, ligandXYZ, ligand, 
										transFac/2, f, 
										chainAndLigand[ligandAndChain[ligand]])
					translate_selection(cXYZ, ligandXYZ, 
											"_tmpPoint" + ligand, 
											transFac/2, f, 
											chainAndLigand[ligandAndChain[ligand]])
					if complex:
						condition = isColliding(ligand, complex) \
								or isColliding(ligand, ligandAndChain[ligand])
					else:
						condition = isColliding(ligand, ligandAndChain[ligand])

			f = f + 10
			cmd.frame(f)
			store_view(group=True, all = True)
				
			f = f + 10
			cmd.frame(f)
			
		## store view again to show final explosion for some time
		store_view(group=True, all = True)

		''' translation of single ligand or residue only '''
	else: 
		## label residue
		resis = str(list(set(stored.res))[0])
		dim = cmd.get_extent(selected)
		cmd.pseudoatom("_tmpPoint" + resis, pos=dim[1])
		cmd.label("_tmpPoint" + resis, "'%s'" %resis)
		cmd.hide('nonbonded', "_tmpPoint" + resis)
		cmd.group(resis + "_", selected + " and " + "_tmpPoint" + resis)
		
		color_binding(complex, selected, res = True)
		
		store_res_view(selected, 1, resis, all = True)

		store_res_view(selected, 10, resis, orient = True)

		store_res_view(selected, 20,resis)
		
		## translate residue
		cmd.frame(50)
		ligandXYZ = calc_COM(selected)
		while isColliding(selected, complex):
			translate_selection(complexXYZ, ligandXYZ, resis + "_", 
									transFac*5, 50, resis + "_")
		cmd.zoom('_inter %s_' %resis)
		cmd.mview('store')
		
		store_res_view(selected, 80)
		

	cmd.frame(160)	
	store_view(all = True)
		
	print time.clock() - start_time, 'seconds'

def canonical_explosion(selected):
	'''DESCRIPTION: translate a selection in canonical direction and create
		a movie from it'''
		
	'''initialize PyMOL'''
	## case sensitive for chain ids
	cmd.set('ignore_case', 'off') 

	start_time = time.clock()

	cmd.remove('solvent')
	cmd.show('spheres', 'organic')
	util.cbc(selection= selected)

	cmd.orient(selected)
	
	'''get dimensions of structure to decide along which axes to translate'''
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

	print transVec	
	
	## initialize movie
	initialize_movie(selected)
	f = 1
	
	''' get chains and ligands of complex '''
	chains = cmd.get_chains(selected)

	allChains = []
	chainAndLabel = {}
	for chain in chains:
		## create an object for every chain 
		chainname = "chain" + chain
		allChains.append(chainname)
		cmd.extract(chainname, selected + " & chain " + chain)
		
		## label chains
		chainAndLabel = label_obj(chainname, chainAndLabel)
		
		## store object for movie
		cmd.frame(f)
		cmd.show('sticks', chainname)
		cmd.mview('store', object=chainname)
		cmd.mview('store', object="_tmpPoint" + chainname)
		
	cmd.zoom('all')
	cmd.mview('store')

	## get ligands'''
	storedLigands = get_ligands()
	
			
	chainAndLigand = {}
	ligandAndChain = {}
	## create objects for ligands
	for chainname in allChains:
		if stored.ligands:
			for l in storedLigands:
			
				## select ligand on chain	
				selection = chainname + ' & resn ' + l
				
				get_ligand_chain_pair(l, selection, chainname, chainAndLigand, 
										ligandAndChain, COM = False)
					
		if chainname not in chainAndLigand.keys():
			cmd.group(chainname + "_", 
						chainname + " " + "_tmpPoint" + chainname)
					
	f = f + 20
	cmd.frame(f)
	## store chain objects for movie 
	store_view(group = True, all = True)
					
	
	''' translate chains'''
	i = 1
	f = f + 40
	cmd.frame(f)
	for chain in allChains:
		for c in allChains:
			if c != chain:
				while isColliding(chain, c):
					for ch in allChains[1:]:
						cmd.translate( [x * i for x in transVec] , object=ch)
						cmd.translate([x * i for x in transVec], 
										object="_tmpPoint" + ch)
						if ch in chainAndLigand.keys():
							cmd.translate([x * i for x in transVec], 
											object=chainAndLigand[ch])
							cmd.translate([x * i for x in transVec], 
											object="_tmpPoint" + \
											chainAndLigand[ch])
						cmd.mview('store', object = ch+"_")
						cmd.mview('interpolate', object = ch+"_")
						i = i + 1
			else:
				continue
				
	store_view(group = True, all = True)
	
	f = f + 20
	cmd.frame(f)
	store_view(group=True, all = True)
		
	''' translate ligands'''
	f = f + 30
	cmd.frame(f)
	
	for ligand in ligandAndChain.keys():
		while isColliding(ligand, ligandAndChain[ligand]):
			cmd.translate([x * 1/4  for x in transVec], object=ligand)
			cmd.translate([x * 1/4 for x in transVec], 
							object="_tmpPoint" + ligand)
	store_view(group = True, all = True)

	f = f + 20
	cmd.frame(f)
	store_view(group=True, all = True)
	
	cmd.zoom('all', complete=1)					
	print time.clock() - start_time, 'seconds'
	
cmd.extend('com_explosion', com_explosion)
cmd.extend('canonical_explosion', canonical_explosion)	
