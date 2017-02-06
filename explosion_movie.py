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
import get_colors ## get random colors
import time
from operator import itemgetter
from viewpoints import best_view
import os.path

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
	cmd.set('movie_fps', 10)
	
	cmd.mset('1 x' + frames)
	if selected:
		cmd.orient(selected)

def remove_solvents(exclude, cutoff, storedLigands = None):
	
	if not os.path.exists('./cc-counts.tdd'):
		import urllib
		urllib.urlretrieve ("http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd", 
								"cc_counts.tdd")
	if not os.path.exists('./cc-counts.tdd'):
		sys.exit('Could not download http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd, please try manually.')
	
	solvents = {}
	aa = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
	cmd.remove('solvent')
	
	import csv
	
	with open('./cc-counts.tdd') as csvfile:
		countreader = csv.reader(csvfile, delimiter='\t')
		for row in countreader:
			if row[0] != 'id':
				ligand = row[0]
				ligcount = int(row[1])
				if ligcount > cutoff and not ligand in aa and \
					not ligand in ['ADP', 'ANP']:
					if not exclude or ligand not in exclude:
						cmd.remove('resn ' + ligand)
						if storedLigands and ligand in storedLigands:
							storedLigands.discard(ligand)
			
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
	
def create_objects(chains, selected, storedLigands, chainsCOMS, complexXYZ, dim, 
					colorBinding, chainAndLigand = None, typeOfExplosion = 'com'):

	label_objec = {}
	
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
	if typeOfExplosion == 'com':
		ligandsCOMS = {}
		chainsCOMS = {}
		
	## store coloring
	objColor = {}
	
	for c in chains:	
		if colorBinding == 'individual':
			## color each chain individually (red is color for binding site)
			col = get_colors.get_random_color()
			while col in objColor.values() or \
					col in ('red','white', 'black', 'marine'):
				col = get_colors.get_random_color()
			cmd.color(col, selected + '& chain ' + c) 
			
		if colorBinding == 'gray':
			## color all chains in gray
			cmd.color('gray', selected + '& chain ' + c)
					
		## create object for chain with name
		chainname = selected + '_chain' + c
		cNames.append(chainname)
		cmd.extract(chainname, selected + ' & chain '+c)
		if colorBinding == 'individual':
			objColor.update({chainname:col})
		# if not storedLigands:
			# cmd.show('cartoon',chainname)
		# else:
			# cmd.show('surface', chainname)

		if typeOfExplosion=='com':
			
			## calculate and save COM of chain
			chainsCOMS[chainname] = calc_COM(chainname)
			chainCOM = chainsCOMS[chainname]
			
			## label chains
			labels = label_obj(chainname, chainCOM, complexXYZ, dim, chainAndLabel)
			if len(labels) > 1:
				chainAndLabel = labels[0]
				label_objects_new = labels[1]
			elif len(labels) == 1:
				label_objects_new = labels
			
			if label_objects_new:
				label_objec.update(label_objects_new)

		## store object for movie
		f = 1
		if typeOfExplosion == 'com':
			cmd.frame(f)
			cmd.mview('store', object=chainname)
			cmd.mview('store', object="label" + chainname)
			cmd.zoom('all', complete = 1)
			cmd.show('dashes')
		
		## create objects for ligands
		if storedLigands:
			for l in storedLigands:
			
				## select ligand on chain	
				selection = chainname + ' & resn ' + l
				
				if typeOfExplosion == 'com':
					label_objec_new = get_ligand_chain_pair(l, selection, 
												chainname, chainAndLigand, 
												ligandAndChain, complexXYZ, dim, 
												label_objec, 'com',ligandsCOMS)
					if label_objects_new:
						label_objec.update(label_objects_new)
				else:
					get_ligand_chain_pair(l, selection, chainname, 
										chainAndLigand, ligandAndChain, 
										complexXYZ, dim, label_objec, 
										typeOfExplosion = 'canonical')
					
					cmd.frame(f)
					cmd.mview('store', object = 'all')
			## if there is no ligand on chain, just keep chain and its label
			if typeOfExplosion == 'com':
				noLigands(chainname, chainAndLigand, label_objec)
			else:
				noLigands(chainname, chainAndLigand)
		## no ligands in structure
		else:
			if typeOfExplosion == 'com':
				noLigands(chainname, chainAndLigand, label_objec)
			else:
				noLigands(chainname, chainAndLigand)
				
	if typeOfExplosion == 'com':	
		return cNames, chainAndLabel, chainAndLigand, ligandAndChain, ligandsCOMS,chainsCOMS, f, label_objec, objColor
	else:
		return cNames, chainAndLabel, chainAndLigand, ligandAndChain, f, label_objec, objColor

def get_ligand_chain_pair(l, selection, chainname, chainAndLigand, 
							ligandAndChain, complexXYZ, dim, label_objects, 
							typeOfExplosion = 'com', ligandsCOMS = None):
	'''DESCRIPTION:
		get pair of ligand and respective chain''' 

	## create object for ligand
	if cmd.count_atoms(selection) > 0:
		ligandname = chainname + '_' + l
		cmd.extract(ligandname, selection)
		# cmd.show('spheres', ligandname)
	
		ligandAndChain[ligandname] = chainname
		
		if typeOfExplosion == 'com':
			## calculate COM of ligand
			ligandsCOMS[ligandname] = calc_COM(ligandname)
			ligandCOM = ligandsCOMS[ligandname]
			## label ligands			
			label_objects_new = label_obj(ligandname, ligandCOM, complexXYZ, dim)
		
			## group chains and respective ligands to translate  
			## them together at first (including their labels)
			cmd.group(chainname + "_" + l + "_", 
							chainname + " " + ligandname + " " + "label" + \
							chainname + " " + "label" + ligandname + " " + \
							label_objects_new[ligandname][0] + " " + \
							label_objects[chainname][0])
			chainAndLigand[chainname].append(chainname + "_" + l + "_")

			label_objects.update(label_objects_new)
			  
			return label_objects
		
		if typeOfExplosion == 'canonical':
			## group chains and respective ligands to translate  
			## them together at first (including their labels)
			cmd.group(chainname + "_" + l + "_", chainname + " " + ligandname)
			chainAndLigand[chainname].append(chainname + "_" + l + "_")			
		
def calc_label_positions_circular(ch,chainCOM, complexXYZ, dim):
	''' DESCRIPTION:
		set labels on sphere around complex object with connection line to com
		of according chain/ligand
	'''
	
	min = dim[0]
	max = dim[1]
	R = math.sqrt((min[0]-max[0])**2 + (min[1]-max[1])**2 + (min[2]-max[2])**2) #euclidean
	R = R/2
	
	P1 = chainCOM
	P2 = complexXYZ
	x0 = P1[0]
	y0 = P1[1]
	z0 = P1[2]
	x1 = P2[0]
	y1 = P2[1]
	z1 = P2[2]
	dx = x1 - x0
	dy = y1 - y0
	dz = z1 -  z0 
	cx = P2[0]
	cy = P2[1]
	cz = P2[2]
	
	a = dx*dx + dy*dy + dz*dz
	b = 2*dx*(x0-cx)+2*dy*(y0-cy)+2*dz*(z0-cz)
	c = cx*cx+cy*cy+cz*cz+x0*x0+y0*y0+z0*z0+-2*(cx*x0+cy*y0+cz*z0)-R*R
	
	discriminant = b**2 - 4*a*c
	
	t = (-b-math.sqrt(discriminant))/(2*a)
	
	x = x0 + t*dx                 
	y = y0 + t*dy                 
	z = z0 + t*dz
	
	n1 = ''.join(map(str,P1))

	cmd.pseudoatom('_pt1' + n1, pos=P1)
	cmd.pseudoatom('_pt2' + n1, pos=[x,y,z])
	cmd.hide('nonbonded', "_pt1" + n1)
	cmd.hide('nonbonded', "_pt2" + n1)
	cmd.distance('_'+ch+'_label', '_pt1' + n1, '_pt2' + n1)
	f = 1
	cmd.frame(f)
	cmd.hide('labels', '_'+ch+'_label')
	cmd.mview('store', object='_'+ch+'_label')
	return [x,y,z], ['_'+ch+'_label']
		
def calc_label_position_flush(chains, transVec, f, chainAndLabel=None):
	chainsComs = {}
	label_objects = {}
	## Determine the extents of empty space regions,
	## Determine the positions of anchor points,
	size = {}
	for ch in chains:
		chainsComs[ch] = calc_COM(ch)
		size[ch] = cmd.count_atoms(ch)
	maxChain = max(size, key=size.get)
	
	maxdim = cmd.get_extent(maxChain)
	cmd.pseudoatom('_min'+maxChain, pos=maxdim[0])
	cmd.hide('nonbonded', '_min'+maxChain)
	cmd.pseudoatom('_max'+maxChain, pos=maxdim[1])
	cmd.hide('nonbonded', '_max'+maxChain)
	labellength = cmd.get_distance('_min'+maxChain,'_max'+maxChain)
	
	## Choose a pivot point
	mean_xpos = sum(chainsComs[ch][0]for ch in chainsComs.keys())/len(chainsComs)
	
	mean_ypos = sum(chainsComs[ch][1] for ch in chainsComs.keys())/len(chainsComs)
	
	mean_zpos = sum(chainsComs[ch][2] for ch in chainsComs.keys())/len(chainsComs)
	
	cmd.pseudoatom('_mean', pos=[mean_xpos,mean_ypos,mean_zpos])
	cmd.hide('nonbonded', '_mean')

	## Assign labels to the left and right region by comparing
	## their anchors with the pivot point.
	for ch in chainsComs:
		cmd.pseudoatom('_com_' + ch, pos = chainsComs[ch])
		cmd.hide('nonbonded', '_com_' + ch)
		
		if transVec[0] > 0 and transVec[1]:
			x_pos = chainsComs[ch][0]
			y_pos = chainsComs[ch][1]
			if chainsComs[ch][2] > mean_zpos:
				z_pos = chainsComs[ch][2] + labellength*1
			else:
				z_pos = chainsComs[ch][2] - labellength*1
			
		elif transVec[1] > 0 and transVec[2] > 0:
			y_pos = chainsComs[ch][1]
			z_pos = chainsComs[ch][2]
			if chainsComs[ch][0] > mean_xpos:
				x_pos = chainsComs[ch][0] + labellength*1
			else:
				x_pos = chainsComs[ch][0] - labellength*1
				
		elif transVec[2] > 0 and transVec[0] > 0:
			x_pos = chainsComs[ch][0]
			z_pos = chainsComs[ch][2]
			if chainsComs[ch][1] > mean_ypos:
				y_pos = chainsComs[ch][1] + labellength*1
			else:
				y_pos = chainsComs[ch][1] - labellength*1
				
		cmd.pseudoatom('_'+ch+'_pos', pos = [x_pos, y_pos, z_pos])
		cmd.distance('_'+ch + '_label','_'+ch+'_pos', '_com_' + ch)
		cmd.hide('nonbonded', '_'+ch+'_pos')	
		cmd.frame(f)
		cmd.pseudoatom('label' + ch, pos = [x_pos, y_pos, z_pos])
		cmd.label("label" + ch, "'%s'" %ch.split('_')[-1])
		cmd.hide('labels', '_'+ch + '_label')
		cmd.hide('nonbonded', 'label' + ch)
		
		chainAndLabel[ch] = "label" + ch
		label_objects.update({ch : '_'+ch + '_label'})
	return [chainAndLabel, label_objects]

def label_obj(chainname, chainCOM, complexXYZ, dim, chainAndLabel = None):
	'''DESCRIPTION:
		create an label for given chain (if chainAndLabel set) or ligand  
	'''
	## create pseudoatom for object to be labeled
	label_pos, lab_obj = calc_label_positions_circular(chainname, chainCOM, 
																complexXYZ, dim)
	
	cmd.pseudoatom("label" + chainname, pos=label_pos)
	
	## label ligand
	if not chainAndLabel:
		ligName = chainname.split("_")[-1]
		cmd.label("label" + chainname, "'%s'" %ligName)
		
	## label chain
	else:
		cmd.label("label" + chainname, "'%s'" %chainname)
		
	## save label of chain in dictonary
	if chainAndLabel:
		chainAndLabel[chainname] = "label" + chainname
		
	## store ligand label in movie
	if not chainAndLabel:
		cmd.mview('store', object = "label" + chainname)
		
	## hide pseudoatom represantation
	cmd.hide('nonbonded', "label" + chainname)
	
	## return dictionary of chain and its label
	if chainAndLabel:
		return [chainAndLabel, {chainname : lab_obj}]
	else:
		return {chainname : lab_obj}

def noLigands(chainname, chainAndLigand, label_objec = None):
	''' DESCRIPTION:
		create group for chain which do not have a ligand
	'''
	if chainname in chainAndLigand.keys() and not chainAndLigand[chainname]:
		del chainAndLigand[chainname]
	if label_objec:
		cmd.group(chainname + "_", chainname + " " + "label" + chainname + " "\
					+ label_objec[chainname][0])
	else:
		cmd.group(chainname + "_", chainname + " " + "label" + chainname)
		
def color_binding(chainname, ligandname, contactColors = None):
	''' DESCRIPTION:
		color binding site in red'''
	binding = 'byres (' + chainname + ' nto. 3.6 of ' + ligandname +')'
	cmd.select('_inter', binding)
	
	if not contactColors:
		cmd.color('red', '_inter')
		
	if contactColors:
		col = get_colors.get_random_color()
		while col in contactColors or col in ('white', 'black', 'marine'):
			col = get_colors.get_random_color()
		cmd.color(col, '_inter') 
		contactColors.append(col)	
	cmd.delete('_inter')

def color_contact(cNames, objColor, colorBinding):
	''' DESCRIPTION:
		color contact site between two chains by color of the other chain
	'''
	if colorBinding == 'gray':
		contactColors = ['gray']
	i = 1
	for chain in cNames:
		for ch in cNames[i:]:
			if ch != chain:
			
				binding = 'byres (' + chain + ' nto. 5 of ' + ch +')'
				cmd.select('_contact1', binding)
				
				binding = 'byres (' + ch + ' nto. 5 of ' + chain +')'
				cmd.select('_contact2', binding)
				
				if colorBinding == 'individual':
					cmd.color(objColor[ch], '_contact1')
					cmd.color(objColor[chain], '_contact2')
					
				if colorBinding == 'gray':
					col = get_colors.get_random_color()
					while col in contactColors or \
									col in ('white', 'black', 'marine'):
						col = get_colors.get_random_color()
					cmd.color(col, '_contact1') 
					cmd.color(col, '_contact2') 
					contactColors.append(col)			
		i+=1
	if colorBinding == 'gray':
		return contactColors

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

def best_view_objects(label = False):
	'''DESCRIPTION: 
		return sting of all objects condsidered for best_view calulation 
	'''
	view_objects = " "
	for obj in cmd.get_names('objects'):
			if not label:
				if not obj.startswith('_') and not obj.startswith('label'):
					view_objects += obj + " "
			else:
				if not obj.startswith('_'):
					view_objects += obj + " "
	return view_objects
				
def translate_selection(originXYZ, transXYZ, transname, factor = 1, f = 1, 
																group = None):
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
	translate_selection(complexXYZ, chainXYZ, cname + '_', transFac, f, cname + "_")
	## if chain contains ligand, translate ligand also
	if cname in chainAndLigand.keys():
		
		for l in chainAndLigand[cname]:
			translate_selection(complexXYZ, chainXYZ, l, transFac, f, l)
	# ## only chain has to be translated
	# else:
		# translate_selection(complexXYZ, chainXYZ, cname + "_", transFac,f, 
								# cname + "_")

def canonical_translation(ch, i, transVec, chainAndLigand, label_objects):
	'''DESCRIPTION:
		translate chain ch canonical
	'''
	## if chain contains ligand, translate ligand also
	if ch in chainAndLigand.keys():
		for l in chainAndLigand[ch]:
			cmd.translate([x * i for x in transVec], object=l, camera=0)
			cmd.mview('interpolate', object=l)
			cmd.mview('store', object = l)	
	## only chain has to be translated
	else:
		cmd.translate([x * i for x in transVec], object=ch + "_", camera=0)
		cmd.mview('interpolate', object = ch+"_")
		cmd.mview('store', object = ch+"_")
		
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

def show_labels(label_objects):
	cmd.show('labels')
	for l in label_objects.values():
		if type(l) is list:
			cmd.hide('labels', l[0])
		else:
			cmd.hide('labels', l)
	cmd.show('dashes')
	
def com_explosion(selected, label_objects, cNames, chainAndLigand, 
					ligandAndChain,	chainsCOMS, ligandsCOMS, transFac, 
					complex = None, com = None,	frame = 1):
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
							com_translation(cname, chainAndLigand, complexXYZ,
												cXYZ, transFac, f)
			store_view(group = True, all = True)
									
		## if only one chain is selected, translate it and its ligand
		else:
			cmd.frame(f)
			cmd.zoom(selected)
			store_view(chainname + '_')
			
			## if source complex is known use its COM to translate
			if complex:
				while isColliding(chainname, complex):
					com_translation(chainname, chainAndLigand, complexXYZ,
												chainXYZ, transFac, f)
			else:
				if chainname in chainAndLigand.keys():
						cmd.translate([transFac, transFac, transFac], 
										object = chainAndLigand[chainname])
				else:
					cmd.translate([transFac, transFac, transFac], 
										object = chainname + "_")
										
			store_view(chainname + '_')

	## store objects for movie
	f = f + 30
	cmd.frame(f)
	if len(cNames) == 1:
		cmd.zoom('all', complete=1)
	store_view(group=True, all = True)

	''' translate ligands '''
	if ligandAndChain:
		f = f + 10
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
				translate_selection(cXYZ, ligandXYZ, ligand, transFac, f, 
									chainAndLigand[ligandAndChain[ligand]])
				translate_selection(cXYZ, ligandXYZ, "label" + ligand, transFac, 
									f, chainAndLigand[ligandAndChain[ligand]])
				translate_selection(cXYZ, ligandXYZ, label_objects[ligand], 
									transFac, f)
				translate_selection(cXYZ, ligandXYZ, label_objects[ligand][0], 
									transFac, f)
				if complex:
					condition = isColliding(ligand, complex) \
							or isColliding(ligand, ligandAndChain[ligand])
				else:
					condition = isColliding(ligand, ligandAndChain[ligand])
		
		f = f + 30
		## show labels
		cmd.frame(f)
		view_objects = " "
		store_view(group=True, all = True)
	return f

def canonical_explosion(selected, label_objects, cNames, chainAndLigand, transVec, 
							ligandAndChain = None, complex = None, frame = 1):
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
					canonical_translation(ch, i, transVec, chainAndLigand, \
											label_objects)
					i += 1
		
	for chain in cNames:
		for c in cNames:
			if c != chain:
				while isColliding(chain, c):
					for ch in cNames[1:]:
						canonical_translation(ch, i, transVec, chainAndLigand,\
												label_objects)
						i += 1
	
	store_view(group = True, all = True)
	
	f = f + 30
	cmd.frame(f)
	store_view(group=True, all = True)
		
	''' translate ligands'''
	if ligandAndChain:
		f = f + 15
		cmd.frame(f)
		for ligand in ligandAndChain.keys():
			
			while isColliding(ligand, ligandAndChain[ligand]):
				cmd.translate([x * 1/4  for x in transVec], object=ligand)
				cmd.translate([x * 1/4 for x in transVec], 
													object="label" + ligand)
				cmd.translate([x * 1/4 for x in transVec], 
													object='_'+ligand+'_label')
			store_view(group = True, all = True)
		
		f = f + 30
		cmd.frame(f)
		store_view(group=True, all = True)

	return f
	
def explosion(selected = ' ', typeOfExplosion = 'com', complex = None, 
				removeSolvents = True, exclude = None, cutoff = 10, 
				colorBinding = 'individual'):
	'''DESCRIPTION:
		perform an explosion of selected object(s) given in a list and create	
		a movie. If two objects are given they will be separateted and then 
		exploded individually.
	'''
	selected = selected.split()
	if exclude:
		exclude = exclude.split()
	if not typeOfExplosion in ['com', 'canonical']:
		sys.exit("Specify explosion: com (default) or canonical.") 
		
	if not selected:
		sys.exit("Please give a list of selected objects (or just one) to translate.") 
		
	if complex:
		if cmd.count_atoms(complex) == 0:
			complex = None
			print 'No need of complex'
	
	if not colorBinding in ['individual', 'gray', 'none']:
		sys.exit("Specify a color for chains: none (do not color), gray, or individual.")
			
	'''setup of selected sturcture'''
	## case sensitive for chain ids
	cmd.set('ignore_case', 'off') 
	# cmd.set('ray_trace_mode', 1)
	# cmd.set('ray_trace_frames', 'on')
	
	cmd.set('dash_color', 'marine')
	cmd.set('dash_round_ends', 'off')
	cmd.set('dash_width', 1)
	cmd.set('dash_gap', 0)
	
	# if complex:
		# cmd.show('surface', complex)
		
	'''preparation and translation'''
	if len(selected) >= 1:
		## initialize movie
		start_time = time.clock()
		if len(selected) > 1:
			initialize_movie(frames = str(100*len(selected)))
		else:
			initialize_movie(frames = str(130))
		
		## get ligands
		storedLigands = get_ligands()
		
		## remove solvent
		if removeSolvents:
			if storedLigands:
				remove_solvents(exclude, cutoff, storedLigands)
			else:
				remove_solvents(exclude, cutoff)

		## lists for further computations
		s = ''
		chains = {}
		cNames = []
		chainAndLabel = {}
		chainAndLigand = {}
		ligandAndChain = {}
		chainsAndObj = {}
		label_objects = {}
		objColor = {}
		
		## calc coms for com-explosion and labeling
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
			dim = cmd.get_extent('_all_obj')
			complexXYZ = calc_COM('_all_obj')
			if typeOfExplosion == 'com':
				## calc com of complex from objects
				trans = calcTransFac('_all_obj')
			else:
				## calculate translation vector
				transVec = transAxes('_all_obj') 
				print "Translation vector:", transVec
			cmd.frame(1)
			cmd.mview('store', object = 'all')
			cmd.delete('_all_obj')
				
		else:
			dim = cmd.get_extent(complex)
			complexXYZ = calc_COM(complex)
			if typeOfExplosion == 'com':
				## calc com of complex from objects
				trans = calcTransFac(complex)
				
			else:
				## calculate translation vector
				transVec = transAxes(complex) 
				print "Translation vector:",transVec
			
		for obj in selected:
			##get chains of complex
			chains[obj] = cmd.get_chains(obj)
			for ch in chains[obj]:
				chainsAndObj[obj+ '_chain' + ch] = obj
				chainAndLigand[obj+ '_chain' + ch] = []
				
			## create objects for all chains and ligands
			if typeOfExplosion == 'com':
				cNames_new, chainAndLabel_new, chainAndLigand_new, \
				ligandAndChain_new, ligandsCOMS_new, chainsCOMS_new, f, \
				label_objects_new, objColor_new = \
					create_objects(chains[obj], obj, storedLigands, chainsCOMS, 
									complexXYZ, dim, colorBinding, chainAndLigand)
				if ligandsCOMS_new:
					ligandsCOMS.update(ligandsCOMS_new)
				if chainsCOMS_new:
					chainsCOMS.update(chainsCOMS_new)
					
			if typeOfExplosion == 'canonical':
				cNames_new, chainAndLabel_new, chainAndLigand_new, \
				ligandAndChain_new, f, label_objects_new, objColor_new=  \
					create_objects(chains[obj], obj, storedLigands, chainsCOMS, 
									complexXYZ, dim, colorBinding, chainAndLigand, 
									typeOfExplosion='canonical')
									
			cNames = cNames + cNames_new
			if chainAndLabel_new:
				chainAndLabel.update(chainAndLabel_new)
			if chainAndLigand_new:
				chainAndLigand.update(chainAndLigand_new)
			if ligandAndChain_new:
				ligandAndChain.update(ligandAndChain_new)
			if label_objects_new:
				label_objects.update(label_objects_new)
			if objColor_new:
				objColor.update(objColor_new)
				
		## color contact sites between chains
		if colorBinding == 'gray':
			contactColors = color_contact(cNames, objColor, colorBinding)	
		elif colorBinding == 'individual':
			color_contact(cNames, objColor, colorBinding)	
			
		## color binding site
		if colorBinding != 'none':
			for cl in chainAndLigand.iteritems():
				for l in cl[1]:
					if colorBinding == 'gray' and contactColors:
						color_binding(cl[0], l[:-1], contactColors)
					elif colorBinding == 'individual':
						color_binding(cl[0], l[:-1])
					
		## set best view
		cmd.frame(f)
		view_objects = best_view_objects()
		best_view(view_objects, 'chain', '10')
		# cmd.orient('all')
		# cmd.zoom('all', complete=1)
		cmd.mview('store', object='all')	
		
		## store chain objects for movie 
		store_view(group = True, all = True)
		
		## labeling for canonical explosion
		if typeOfExplosion == 'canonical':
			chainslist = cNames
			ligandslist = ligandAndChain.keys()
			chainslist.extend(ligandslist)
			chainAndLabel, label_objects = \
					calc_label_position_flush(chainslist, transVec, 1, 
												chainAndLabel)
			for ch in cNames:
				if ch in chainAndLigand.keys() and chainAndLigand[ch]:
					for l in chainAndLigand[ch]:
						cmd.group(l, label_objects[ch]+  " " + \
						label_objects[l[:-1]] + " " + "label" + ch + " " + \
						"label" + l[:-1] + " " + '_'+ch + '_label' + " " +\
						'_'+l[:-1] + '_label', 'add') 
				else:
					cmd.group(ch + '_', "label"+ch + " "+'_'+ch + '_label','add')
			cmd.frame(1)
			cmd.zoom('all', complete=1)
			cmd.mview('store', object='all')
			cmd.scene('on', 'store')
			cmd.mview('store', scene='on')

		## show labels
		if typeOfExplosion == 'com':
			cmd.frame(f-20)
			show_labels(label_objects)
			cmd.scene('on', 'store')
			cmd.mview('store', scene='on')
			cmd.mview('store', object='all')
		
		f = f + 20
		cmd.frame(f)	
		cmd.zoom('all', complete=1)
		cmd.mview('store', object='all')
		cmd.mview('store', scene='on')
		
		## hide labels
		f = f + 1
		cmd.frame(f)
		cmd.hide('labels')
		cmd.hide('dashes')
		cmd.scene('off', 'store')
		cmd.zoom('all',complete=1)
		cmd.mview('store', scene='off')
			
		f = f + 10		
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
									com_translation(cname, chainAndLigand, 
													complexXYZ, 
													coms[chainsAndObj[cname]], 
													trans, f)
								
							else:
								for o in chains.keys():
									for ch in chains[o]:
										ch = obj+ '_chain' + ch
										canonical_translation(ch, i, transVec, 
													chainAndLigand, label_objects)
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
				f = f + 10
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
					ligandAndChain_obj.update({ligand: c for ligand, 
										c in ligandAndChain.items() if c == ch})
			
			## translation of object
			if typeOfExplosion == 'com':
				if not complex:
					f = com_explosion('_'+obj, label_objects, objChains,
								chainAndLigand, ligandAndChain_obj, chainsCOMS, 
								ligandsCOMS, trans, com = coms[obj], frame = f)
								
				else:
					f = com_explosion('_'+obj, label_objects, objChains,
								chainAndLigand,	ligandAndChain_obj, chainsCOMS, 
								ligandsCOMS, trans, complex = complex,
								com = complexXYZ, frame = f)

			else:
				if complex:
					f = canonical_explosion('_'+obj, label_objects, objChains, 
							chainAndLigand, trans,ligandAndChain_obj, complex, f)
							
				else:
					f = canonical_explosion('_'+obj, label_objects, objChains, 
							chainAndLigand,trans, ligandAndChain_obj, frame = f)
				
				
			cmd.delete('_'+obj)

	cmd.frame(f-29)
	show_labels(label_objects)
	cmd.scene('on2', 'store')
	cmd.mview('store', scene='on2')
	cmd.frame(f)
	view_objects =best_view_objects(True)
	best_view(view_objects, 'chain', '10')
	cmd.zoom('all', complete=1)
	cmd.mview('store', scene='on2')
		
	print 'Explosion of', selected, time.clock() - start_time, 'seconds'
			
def relabel(selected, newLabel="new label"):
	'''DESCRIPTION:
		rename an selected object and its label by a new label
	'''
	obj = cmd.get_names("objects")
	x = ''
	
	## relabel label
	cmd.label("label" + selected , "'%s'" %newLabel )
	
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
''' TODO: nach neuer Repräsentation werden labels immer angezeigt '''	
def renew_representation(selection, representation):
	'''DESCRIPTION:
		show selection in given representation
	'''
	cmd.show_as(representation, selection)
	cmd.scene(selection+'_'+representation, 'store')
	for f in range(1, cmd.count_frames()):
		cmd.frame(f)
		cmd.mview('store', scene=selection+'_'+representation)
		cmd.mview('reinterpolate')
		
def hide_labels():
	cmd.hide('labels')
	cmd.hide('dashes')

def labels_on():
	cmd.show('dashes')
	for obj in cmd.get_names('objects'):
		if obj.startswith('label'):
			cmd.show('labels', obj)
		
def remove_solvent(name):
	for obj in cmd.get_names('objects'):
		if name in obj and obj[-1] != '_':
			cmd.delete(obj)
			
cmd.extend('explosion', explosion)
cmd.extend('relabel', relabel)
cmd.extend('reorient_explosion', reorient_explosion)
cmd.extend('renew_representation', renew_representation)
cmd.extend('show_labels', labels_on)
cmd.extend('hide_labels', hide_labels)
cmd.extend('remove_solvent', remove_solvent)
