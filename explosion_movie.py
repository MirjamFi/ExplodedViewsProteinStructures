<<<<<<< HEAD
# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_movie.py

# ''' load complex from PDB'''
# selected = '3oaa'
# pdbFile = selected + '.pdb1'
# cmd.fetch(selected, type='pdb1')
# cmd.load('C:\Users\Figaschewski\Dropbox\Masterarbeit\TPC2_M484L_wActivator_56ns\TPC2_M484L_wActivator_56ns.pdb', 'test')

'''
	Create an movie of an exploded view for a given protein structure (PDB) or 
	parts of it.
	You can select single or multiple chains, a residue or ligand to be exploded, 
	just extract them into one object from source structure. If chains contain 
	ligands, those 	will also be exploded.
	To explod specific chains/residue/ligand it is useful to also give the 
	complex for explosion directon, espacially for single chain, so the chain is 
	not	translated into the complex.
	Exploded parts will be labeled.
'''
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time


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
		calculate COM of an object 
	'''
	cenma.com(transname, state = 1)
	cmd.zoom(transname + '_COM')
	pos = cmd.get_position(transname + '_COM')
	cmd.delete(transname + '_COM')
	cmd.zoom('all', complete = 1)
	cmd.mview('store')

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


def calcTransFac(sele):
	'''DESCRIPTION 
		calculate translation factor 
	'''
	protDim_min, protDim_max = cmd.get_extent(sele)
	transFac = math.sqrt((protDim_max[0] - protDim_min[0]) ** 2 +
					 (protDim_max[1] - protDim_min[1]) ** 2 +
					 (protDim_max[2] - protDim_min[2]) ** 2)/2
	return transFac
	
def explosion(selected, complex = None):
	''' DESCRIPTION:
		create a movie for an exploded view of given protein
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

	''' get ligands'''
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
	
	storedLigands = set(stored.ligands) 

	'''initialize movie'''
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
			
	''' Explosion: for each state create an object of every chain and for 
		related ligands. In this state calulate COM for state and according 
		chains and translate chain along vector between COMs of state and chain. 
		For every chain do the same with its ligands.
	'''

	## create models for all states of complex
	states = cmd.split_states(selected)

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
		## if only some chains are selected and source complex is given,
		## select source complex for COM calculation
		if complex:
			complexXYZ = calc_COM(complex)
		else:
			complexXYZ = calc_COM(selected + '_' + s)

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
		
		## chain or complex shall be translated
		stored.res=[]
		cmd.iterate('(' + selected + ')',"stored.res.append(resn + resi)")
		if len(set(stored.res)) > 1:
			cmd.hide("(" + selected + ")")
			for c in chains:	
				## color each chain individually
				if cmd.count_states() > 1:
					cmd.color(get_colors.get_random_color(), selected + '_' + \
								s + '& chain ' + c)

				## create object for chain with name
				chainname = selected + '_' + s + '_'+c
				cNames.append(chainname)
				origin = selected+'_'+s+ ' & chain '+c
				cmd.extract(chainname, origin)
				
				## label chains
				dim = cmd.get_extent(chainname)
				cmd.pseudoatom("_tmpPoint" + chainname, pos=dim[0])
				cmd.pseudoatom("_tmpPoint" + chainname, pos=dim[1])
				cmd.label("_tmpPoint" + chainname, "'chain %s'" %c)
				chainAndLabel[chainname] = "_tmpPoint" + chainname
				
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
						
						## get pair of ligand and respective chain
						if cmd.count_atoms(selection) > 0:
							ligandname = chainname + '_' + l
							cmd.extract(ligandname, selection)
							stored.lig=[]
							cmd.iterate('(' + selected + ')',
											"stored.lig.append(resn + resi)")
							ligres = str(list(set(stored.res))[0])
							
							
							## label ligands			
							dim = cmd.get_extent(ligandname)
							ligName = ligandname.split("_")[-1]
							cmd.pseudoatom("_tmpPoint" + ligandname, pos=dim[0])
							cmd.label("_tmpPoint" + ligandname, "'%s'" %ligres)
							cmd.mview('store', object = "_tmpPoint" + ligandname)
							
							## calculate COM of ligand
							ligandsCOMS[ligandname] = calc_COM(ligandname)
							
							## color binding site
							binding = 'byres (' + chainname + ' nto. 3 of ' + \
										ligandname +')'
							cmd.select('inter', binding)
							cmd.color('red', 'inter')
							cmd.delete('inter')
							
							## group chains and respective ligands to translate  
							## them together at first (including their labels)
							cmd.group(chainname + "_" + ligandname, 
											chainname + " " + ligandname \
											+ " " + "_tmpPoint" + chainname \
											+ " " + "_tmpPoint" + ligandname)
							chainAndLigand[chainname] = chainname + \
															"_" + ligandname
							ligandAndChain[ligandname] = chainname
				
				## if there is no ligand on chain, just keep chain and its label
				if chainname not in chainAndLigand.keys():
					cmd.group(chainname + "_", chainname + " and " + "_tmpPoint" 
								+ chainname)
			f = f + 10
			cmd.frame(f)
			## store chain objects for movie
			for group in cmd.get_names("objects"):
				cmd.mview('store', object = group) 
			cmd.zoom('all', complete=1)
			cmd.mview('store')
			
			f = f + 10
			## translate chains iteratively 
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
					cmd.orient(chainname)
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
											
				## if only one chain is selected, translate it and its ligand
				else:
					cmd.frame(f)
					cmd.orient(chainname)
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
					f = f + 10
					cmd.frame(f)
					cmd.orient(chainname)
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
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
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
							
			
			## store objects for movie
			cmd.frame(f)
			for group in cmd.get_names("objects"):
				cmd.mview('store', object = group)
			cmd.mview('store')
			

			f = f + frameNum/2
			cmd.frame(f)
			cmd.zoom('all', complete = 1 )
			cmd.mview('store')
			
			## translate ligands
			if storedLigands:
				f = f + frameNum/2
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
				## store view
				for group in cmd.get_names("objects"):
					cmd.mview('store', object = group)
				cmd.zoom('all', complete = 1)
				cmd.mview('store')
					
				f = f + 10
				cmd.frame(f)
			## store view again to show final explosion for some time
			for group in cmd.get_names("objects"):
				cmd.mview('store', object = group)	
			cmd.zoom('all', complete = 1)
			cmd.mview('store')
			## delete state object
			cmd.delete(selected + '_'+ s)	
			

		## translation of single ligand only
		else: 
			resis = str(list(set(stored.res))[0])
			dim = cmd.get_extent(selected)
			cmd.pseudoatom("_tmpPoint" + resis, pos=dim[1])
			cmd.label("_tmpPoint" + resis, "'%s'" %resis)
			cmd.group(resis + "_", selected + " and " + "_tmpPoint" + resis)
			
			cmd.delete(selected + '_' + s)
			binding = 'byres (' + complex + ' nto. 3 of ' + \
										selected +')'
			cmd.select('_inter', binding)
			cmd.color('red', '_inter')
			
			f = 1
			cmd.frame(f)
			cmd.zoom('all', complete = 1)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			
			f = f + 10
			cmd.frame(f)
			cmd.orient('_inter %s' %selected)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			
			f = f + 10
			cmd.frame(f)
			cmd.orient('_inter %s' %selected)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			
			f = f + frameNum
			cmd.frame(f)
			ligandXYZ = calc_COM(selected)
			while isColliding(selected, complex):
				translate_selection(complexXYZ, ligandXYZ, resis + "_", 
										transFac*5, f, resis + "_")
			cmd.zoom('_inter %s' %selected)
			cmd.mview('store')
			
			f = f + frameNum
			cmd.frame(f)
			cmd.zoom('all', complete=1)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			

	cmd.frame(160)	
	cmd.zoom('all', complete = 1)
	cmd.mview('store')
		
	print time.clock() - start_time, 'seconds'

cmd.extend('explosion', explosion)
=======
# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_movie.py

# ''' load complex from PDB'''
# selected = '3oaa'
# pdbFile = selected + '.pdb1'
# cmd.fetch(selected, type='pdb1')
# cmd.load('C:\Users\Figaschewski\Dropbox\Masterarbeit\TPC2_M484L_wActivator_56ns\TPC2_M484L_wActivator_56ns.pdb', 'test')

'''
	Create an movie of an exploded view for a given protein structure (PDB) or 
	parts of it.
	You can select single or multiple chains, a residue or ligand to be exploded, 
	just extract them into one object from source structure. If chains contain 
	ligands, those 	will also be exploded.
	To explod specific chains/residue/ligand it is useful to also give the 
	complex for explosion directon, espacially for single chain, so the chain is 
	not	translated into the complex.
	Exploded parts will be labeled.
'''
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time


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
		calculate COM of an object 
	'''
	cenma.com(transname, state = 1)
	cmd.zoom(transname + '_COM')
	pos = cmd.get_position(transname + '_COM')
	cmd.delete(transname + '_COM')
	cmd.zoom('all', complete = 1)
	cmd.mview('store')

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


def calcTransFac(sele):
	'''DESCRIPTION 
		calculate translation factor 
	'''
	protDim_min, protDim_max = cmd.get_extent(sele)
	transFac = math.sqrt((protDim_max[0] - protDim_min[0]) ** 2 +
					 (protDim_max[1] - protDim_min[1]) ** 2 +
					 (protDim_max[2] - protDim_min[2]) ** 2)/2
	return transFac
	
def explosion(selected, complex = None):
	''' DESCRIPTION:
		create a movie for an exploded view of given protein
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

	''' get ligands'''
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
	
	storedLigands = set(stored.ligands) 

	'''initialize movie'''
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
			
	''' Explosion: for each state create an object of every chain and for 
		related ligands. In this state calulate COM for state and according 
		chains and translate chain along vector between COMs of state and chain. 
		For every chain do the same with its ligands.
	'''

	## create models for all states of complex
	states = cmd.split_states(selected)

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
		## if only some chains are selected and source complex is given,
		## select source complex for COM calculation
		if complex:
			complexXYZ = calc_COM(complex)
		else:
			complexXYZ = calc_COM(selected + '_' + s)

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
		
		## chain or complex shall be translated
		stored.res=[]
		cmd.iterate('(' + selected + ')',"stored.res.append(resn + resi)")
		if len(set(stored.res)) > 1:
			cmd.hide("(" + selected + ")")
			for c in chains:	
				## color each chain individually
				if cmd.count_states() > 1:
					cmd.color(get_colors.get_random_color(), selected + '_' + \
								s + '& chain ' + c)

				## create object for chain with name
				chainname = selected + '_' + s + '_'+c
				cNames.append(chainname)
				origin = selected+'_'+s+ ' & chain '+c
				cmd.extract(chainname, origin)
				
				## label chains
				dim = cmd.get_extent(chainname)
				cmd.pseudoatom("_tmpPoint" + chainname, pos=dim[0])
				cmd.pseudoatom("_tmpPoint" + chainname, pos=dim[1])
				cmd.label("_tmpPoint" + chainname, "'chain %s'" %c)
				chainAndLabel[chainname] = "_tmpPoint" + chainname
				
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
						
						## get pair of ligand and respective chain
						if cmd.count_atoms(selection) > 0:
							ligandname = chainname + '_' + l
							cmd.extract(ligandname, selection)
							stored.lig=[]
							cmd.iterate('(' + selected + ')',
											"stored.lig.append(resn + resi)")
							ligres = str(list(set(stored.res))[0])
							
							
							## label ligands			
							dim = cmd.get_extent(ligandname)
							ligName = ligandname.split("_")[-1]
							cmd.pseudoatom("_tmpPoint" + ligandname, pos=dim[0])
							cmd.label("_tmpPoint" + ligandname, "'%s'" %ligres)
							cmd.mview('store', object = "_tmpPoint" + ligandname)
							
							## calculate COM of ligand
							ligandsCOMS[ligandname] = calc_COM(ligandname)
							
							## color binding site
							binding = 'byres (' + chainname + ' nto. 3 of ' + \
										ligandname +')'
							cmd.select('inter', binding)
							cmd.color('red', 'inter')
							cmd.delete('inter')
							
							## group chains and respective ligands to translate  
							## them together at first (including their labels)
							cmd.group(chainname + "_" + ligandname, 
											chainname + " " + ligandname \
											+ " " + "_tmpPoint" + chainname \
											+ " " + "_tmpPoint" + ligandname)
							chainAndLigand[chainname] = chainname + \
															"_" + ligandname
							ligandAndChain[ligandname] = chainname
				
				## if there is no ligand on chain, just keep chain and its label
				if chainname not in chainAndLigand.keys():
					cmd.group(chainname + "_", chainname + " and " + "_tmpPoint" 
								+ chainname)
			f = f + 10
			cmd.frame(f)
			## store chain objects for movie
			for group in cmd.get_names("objects"):
				cmd.mview('store', object = group) 
			cmd.zoom('all', complete=1)
			cmd.mview('store')
			
			f = f + 10
			## translate chains iteratively 
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
					cmd.orient(chainname)
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
											
				## if only one chain is selected, translate it and its ligand
				else:
					cmd.frame(f)
					cmd.orient(chainname)
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
					f = f + 10
					cmd.frame(f)
					cmd.orient(chainname)
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
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
					cmd.zoom('all', complete=1)
					cmd.mview('store', object = chainname + '_')
					cmd.mview('store')
							
			
			## store objects for movie
			cmd.frame(f)
			for group in cmd.get_names("objects"):
				cmd.mview('store', object = group)
			cmd.mview('store')
			

			f = f + frameNum/2
			cmd.frame(f)
			cmd.zoom('all', complete = 1 )
			cmd.mview('store')
			
			## translate ligands
			if storedLigands:
				f = f + frameNum/2
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
				## store view
				for group in cmd.get_names("objects"):
					cmd.mview('store', object = group)
				cmd.zoom('all', complete = 1)
				cmd.mview('store')
					
				f = f + 10
				cmd.frame(f)
			## store view again to show final explosion for some time
			for group in cmd.get_names("objects"):
				cmd.mview('store', object = group)	
			cmd.zoom('all', complete = 1)
			cmd.mview('store')
			## delete state object
			cmd.delete(selected + '_'+ s)	
			

		## translation of single ligand only
		else: 
			resis = str(list(set(stored.res))[0])
			dim = cmd.get_extent(selected)
			cmd.pseudoatom("_tmpPoint" + resis, pos=dim[1])
			cmd.label("_tmpPoint" + resis, "'%s'" %resis)
			cmd.group(resis + "_", selected + " and " + "_tmpPoint" + resis)
			
			cmd.delete(selected + '_' + s)
			binding = 'byres (' + complex + ' nto. 3 of ' + \
										selected +')'
			cmd.select('_inter', binding)
			cmd.color('red', '_inter')
			
			f = 1
			cmd.frame(f)
			cmd.zoom('all', complete = 1)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			
			f = f + 10
			cmd.frame(f)
			cmd.orient('_inter %s' %selected)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			
			f = f + 10
			cmd.frame(f)
			cmd.orient('_inter %s' %selected)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			
			f = f + frameNum
			cmd.frame(f)
			ligandXYZ = calc_COM(selected)
			while isColliding(selected, complex):
				translate_selection(complexXYZ, ligandXYZ, resis + "_", 
										transFac*5, f, resis + "_")
			cmd.zoom('_inter %s' %selected)
			cmd.mview('store')
			
			f = f + frameNum
			cmd.frame(f)
			cmd.zoom('all', complete=1)
			cmd.mview('store', object = resis + "_")
			cmd.mview('store')
			

	cmd.frame(160)	
	cmd.zoom('all', complete = 1)
	cmd.mview('store')
		
	print time.clock() - start_time, 'seconds'

cmd.extend('explosion', explosion)
>>>>>>> origin/Explosion
cmd.extend('isColliding', isColliding)