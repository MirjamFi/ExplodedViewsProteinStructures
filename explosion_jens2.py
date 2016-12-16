# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion_jens2.py
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
		
cmd.remove('solvent resn  CL')

cmd.select('popc', 'resn POPC')

cmd.show('spheres', 'organic')

list_popc = []
cmd.iterate('popc','list_popc.append((resi,resn))')
popc = set(list_popc)

##get chains of complex
chains = cmd.get_chains('test')
chains = chains[1:]

## initialize mmovie
cmd.set('matrix_mode', 1)
cmd.set('movie_panel', 1)
cmd.set('scene_buttons', 1)
cmd.set('cache_frames', 1)
cmd.config_mouse('three_button_motions', 1)
# cmd.set('movie_panel', 0)	## hide movie panel
cmd.set('movie_panel_row_height', 1)

frames = str(100)
cmd.mset('1 x' + frames)

start_time = time.clock()

## calculate which chains are near to each other to create pores
i = 1
pore = {}
chainPOPC = []
for chain in chains:
	cmd.color(get_colors.get_random_color(), 'chain ' + chain) 
	if cmd.count_atoms('resn POPC and chain ' + chain) > 0:
		chainPOPC.append(chain)
		chains.remove(chain)
		continue
		
	pore[chain] = []
	chainname = 'chain ' + chain
	for ch in chains[i:]:
		ligandname = "chain " + ch
			
		binding = 'byres (' + chainname + ' nto. 3 of ' + ligandname +')'
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
		
## separate layers
popc_plus = 'resn POPC and '
popc_minus = 'resn POPC and '
min_dim, max_dim = cmd.get_extent('(popc)')
x = abs(min_dim[0]-max_dim[0])
y = abs(min_dim[1]-max_dim[1])
z = abs(min_dim[2]-max_dim[2])
ax = [x,y,z]
separationLine = ax.index(min(ax))
for pop in popc:
	a = 'name O9 and resn POPC and resi ' + str(pop[0]) 
	acoords = cmd.get_coords(a).tolist()[0]
	p = str(pop[0])
	if acoords[separationLine] < (min(ax)/2 + min_dim[separationLine]):
		popc_plus = popc_plus + 'resi ' + p + ' '
	else:
		popc_minus = popc_minus + 'resi ' + p + ' '

print 'pores and layers: ' + str(time.clock() - start_time), 'seconds'

## create movie of layer translations
cmd.extract('popc_lower', popc_minus)
cmd.show_as('mesh', 'popc_lower')
cmd.extract('popc_upper', popc_plus)
cmd.show_as('mesh', 'popc_upper')
cmd.frame(1)
cmd.mview('store', object='popc_upper')
cmd.mview('store', object='popc_lower')
cmd.orient(pores[0])
cmd.zoom('all', complete=1)
cmd.mview('store')

cmd.frame(20)
cmd.mview('store', object='popc_lower')
cmd.mview('store', object='popc_upper')
cmd.zoom('all', complete=1)
cmd.mview('store')

cmd.frame(50)
cmd.translate([0,-100,0], object = 'popc_lower')
cmd.mview('store', object='popc_lower')
cmd.mview('interpolate', object='popc_upper')
cmd.translate([0,100,0], object = 'popc_upper')
cmd.mview('interpolate', object='popc_lower')
cmd.mview('store', object='popc_upper')
cmd.zoom('all', complete=1)
cmd.mview('store')

cmd.frame(70)
cmd.mview('store', object='popc_upper')
cmd.mview('store', object='popc_lower')
cmd.zoom('all', complete=1)
cmd.mview('store')
	
print 'translation: ' + str(time.clock() - start_time), 'seconds'			