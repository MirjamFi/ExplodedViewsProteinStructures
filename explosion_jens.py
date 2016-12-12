# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_jens.py
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time
from operator import itemgetter

cmd.reinitialize()
start_time = time.clock()
cmd.load('C:\Users\Figaschewski\Dropbox\Masterarbeit\TPC2_M484L_wActivator_56ns\TPC2_M484L_wActivator_56ns.pdb', 'test')
print 'load: ' + str(time.clock() - start_time), 'seconds'


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
		
cmd.remove('solvent resn  CL')

cmd.select('popc', 'resn POPC')

cmd.show('spheres', 'organic')

list_popc = []
cmd.iterate('popc','list_popc.append((resi,resn))')
popc = set(list_popc)

##get chains of complex
chains = cmd.get_chains('test')
chains = chains[1:]

cmd.set('matrix_mode', 1)
cmd.set('movie_panel', 1)
cmd.set('scene_buttons', 1)
cmd.set('cache_frames', 1)
cmd.config_mouse('three_button_motions', 1)
# cmd.set('movie_panel', 0)	## hide movie panel
cmd.set('movie_panel_row_height', 0.9)

frames = str((len(popc)*2)+100)
cmd.mset('1 x' + frames)

start_time = time.clock()

## calculate which chains are near to each other to create pores
i = 1
pore = {}
chainPOPC = []
for chain in chains:
	if cmd.count_atoms('resn POPC and chain ' + chain) > 0:
		chainPOPC.append(chain)
		chains.remove(chain)
		continue
		
	pore[chain] = []
	chainname = 'chain ' + chain
	for ch in chains[i:]:
	
		## color chains individually
		cmd.color(get_colors.get_random_color(), 'chain ' + ch)
		
		chain2 = "chain " + ch
			
		binding = 'byres (' + chainname + ' nto. 3 of ' + chain2 +')'
		cmd.select('_' + chain + ch + '_', binding)
		
		countAtoms = cmd.count_atoms('_%s'%chain + ch+'_')
		
		if countAtoms > 0:
			inPore = False
			for c in pore.keys():
				if chain in pore[c]:
					if ch not in pore[c]: 
						pore[c].append(ch)
					inPore = True
			if not inPore:
				pore[chain].append(ch)
							
		else:
			cmd.delete(chain + ch + '_')
			
	if len(pore[chain]) == 0:
		del pore[chain]
		
	i += 1
print 'chains: ' + str(time.clock() - start_time), 'seconds'

start_time = time.clock()
cmd.frame(1)
## create pores	
pores = []
for chain in pore.keys():
	if len(pore[chain]) > 1:
		porename = chain
		sel = "chain " + chain
		for c in pore[chain]:
			porename += c
			sel += " chain "+ c
		pores.append(porename)
		cmd.create(porename, sel)
		cmd.mview('store', object = porename)
		
## calculate COM of pores		
poresCOMs = {}
for p in pores:
	poresCOMs[p] = calc_COM(p)
	
	## create pseudoatom for distance calculation
	cmd.pseudoatom("_" + p, pos=poresCOMs[p])
	cmd.hide('nonbonded', "_" + p)

## calculate distance of every pore to every POPC
popcDist = {}
for pop in popc:
	popcDist[str(pop[0])] = 0

plus = []
minus = []
for pop in popc:
	a = 'name O9 and resn POPC and resi ' + str(pop[0]) 
	acoords = cmd.get_coords(a).tolist()[0]
	if acoords[2] < 80:
		plus.append(str(pop[0]))
		
	## create pseudoatom for distance calculation
	cmd.pseudoatom('_'+str(pop[0]), pos = acoords)
	cmd.hide('nonbonded', '_'+str(pop[0]))
	
	for pore in poresCOMs.keys():
		d = cmd.get_distance('_'+str(pop[0]), "_" + pore)
		popcDist[str(pop[0])] = d + popcDist[str(pop[0])]

## sort lipids for ordered explosion
orderedPopc = sorted(popcDist,  key=popcDist.__getitem__, reverse=True)	

print 'pores and distances: ' + str(time.clock() - start_time), 'seconds'

start_time = time.clock()	

## translate popc
i = 20
cmd.group('POPC')
for p in orderedPopc:
	cmd.frame(1)
	cmd.orient(pores[0])
	cmd.zoom('all', complete=1)
	cmd.mview('store')
	cmd.create(p+'_', 'resn POPC and resi '+ p)
	cmd.show_as('mesh', p+'_')
	cmd.group('POPC', p+'_', 'add')
	cmd.mview('store', object = p+'_')
	a = 'name O9 and resn POPC and resi ' + p 
	acoords = cmd.get_coords(a).tolist()[0]
	cmd.frame(i)

	## upper layer 
	if p in plus:	
		cmd.translate([0,100,0], object = p+'_')
	
	## lower layer
	else:
		cmd.translate([0,-100,0], object = p+'_')
	cmd.mview('store', object=p+'_')
	cmd.mview('interpolate', object=p+'_')
	p2 = 0
	
	## store all lipids
	while orderedPopc[p2] != p:
		cmd.mview('store', object=orderedPopc[p2])
		p2 += 1
		
	cmd.orient(pores[0])
	cmd.zoom('all', complete=1)
	cmd.mview('store')
	i += 2
	
print 'translation: ' + str(time.clock() - start_time), 'seconds'




	
			