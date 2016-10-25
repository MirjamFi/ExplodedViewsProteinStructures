# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion.py

from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
from pymol import util ## color by chain
import get_colors ## get random colors
import time

start_time = time.clock()

cmd.delete('all')
cmd.set('ignore_case', 'off') ## case sensitive for chain ids

''' load complex from PDB'''
pdbName = '5lmo'
pdbFile = pdbName + '.pdb1'
cmd.fetch(pdbName, type='pdb1')
#cmd.set('all_states', 'on')
cmd.remove('solvent')
cmd.show('spheres', 'organic')
util.cbc(selection= pdbName)

''' get chains of complex'''
chains = cmd.get_chains(pdbName)

''' get all ligands (organic) from complex'''
stored.ligands = []

## create the original selection of all organic atoms 
cmd.select("allOrg", "org")

##" temporary selection
cmd.select("temp", "none")

""" while initial selection is not empty, "pop" from it and query the atom's 
 molecule's name """
while cmd.count_atoms("allOrg") != 0:

	## pop--peek, rather
	cmd.select("temp", "first allOrg")
	
	## store the residue name
	cmd.iterate("temp", "stored.ligands.append(resn)")
	
	## remove by molecule
	cmd.select("allOrg", "allOrg and not (bm. temp)")
	
cmd.delete('allOrg')
cmd.delete('temp')

''' translate an object relative to complex of origin using center of mass'''
def translate_selection(originXYZ, transXYZ, transname, factor):
	
	## vector between COMs to translate chain
	dist=math.sqrt((transXYZ[0]-originXYZ[0])**2 +
					(transXYZ[1]-originXYZ[1])**2 +
					(transXYZ[2]-originXYZ[2])**2)

	vector=((transXYZ[0]-originXYZ[0])/dist, 
			(transXYZ[1]-originXYZ[1])/dist, 
			(transXYZ[2]-originXYZ[2])/dist)
	trans_vec = [x * factor for x in vector] 

	## translate chain with vector
	cmd.translate(trans_vec,transname, camera=0)
	
'''calculate COM of an object '''
def calc_COM(transname):
		cenma.com(transname, state = 1)
		cmd.zoom(transname + '_COM')
		pos = cmd.get_position(transname + '_COM')
		cmd.delete(transname + '_COM')
		return pos
		
''' Explosion'''

## create models for all states of complex
states = cmd.split_states(pdbName)

"""for each state create an object of every chain and for related ligands. 
	In this state calulate COM for state and according chains and translate 
	chain along vector between COMs of state and chain. 
	For every chain do the same with its ligands.
"""
com_num = ""
for state in range(1,cmd.count_states()+1):
	## state = 0 would calculate the COM for all states at once
	if state == 0: 
		continue

	## name chain object in this state
	s = "00" + str(state)
	if state < 10:
		s = "0" + s
		
	## COM of state (all chains)
	complexXYZ = calc_COM(pdbName + "_" + s)

	## calculate COM for single chains in state
	cNames = []
	chainsCOMS = {}
	for c in chains:
		## color each chain individually
		if cmd.count_states() > 1:
			cmd.color(get_colors.get_random_color(), pdbName + "_" + s + "& chain " + c)

		## create object for chain with name
		chainname = s + '_'+c
		cNames.append(chainname)
		origin = pdbName+'_'+s+ ' & chain '+c
		cmd.create(chainname, origin)
		
		## calculate and save COM of chain
		chainsCOMS[chainname] = calc_COM(chainname)		
	
	## translate chains iteratively to asure a distance of 3 Angstrom between chains
	for chainname in cNames:
		## COM of chain
		chainXYZ = chainsCOMS[chainname]
		
		## list of all chains to calculate distances after translation
		chainsDist = dict.fromkeys(cNames, 100)  
		chainsDist[chainname] = 0

		## only translate chain if there are multiple chains
		if len(chains) > 1:	
		
			## translate chain, so that it has at least 3Angstrom distance to
			## all other chains
			while any(x>0 for x in chainsDist.values()):
				for key in chainsDist.keys():
					if key != chainname:
					
						## number of atoms in radius of 3 Angstrom
						nearC = key + ' w. 3 of ' + chainname
						chainsDist[key] = cmd.count_atoms(nearC)
						
						## if number of atoms > 0, translate all chains 
						## (to keep shape of complex)
						if cmd.count_atoms(nearC) > 0:
		
							## translate chains
							for cname in cNames:
								cXYZ = chainsCOMS[cname]
								translate_selection(complexXYZ, cXYZ, cname, 20)
								cmd.delete(cname + '_COM')

	## translate ligands (also distance of 3 Angstrom to chains)
	for chainname in cNames:
			
		## translate ligands on current chain, if ligands are present
		if stored.ligands:
			
			## select ligand on chain	
			selection = chainname + " & organic"
			ligandname = chainname + "_org"

			## only translate if chain really has a ligand
			if cmd.count_atoms(selection) > 0:
				cmd.create(ligandname, selection)
		
				## calculate COM of ligand
				ligandXYZ = calc_COM(ligandname)
				
				## if chain got translated, claculate new COM for translation of ligands
				if len(chains) > 1:
					chainsCOMS[chainname] = calc_COM(chainname)
					chainXYZ = chainsCOMS[chainname]
			
				## number of atoms in radius of 3 Angstrom
				nearC = ligandname + ' w. 3 of ' + chainname
				
				## asure distance of 3 Angstrom between ligand and chain
				while cmd.count_atoms(nearC) > 0:
					## translate ligand
					translate_selection(chainXYZ, ligandXYZ, ligandname, 10)
					nearC = ligandname + ' w. 3 of ' + chainname
			
				## color binding site
				binding = 'byres (' + chainname + ' nto. 3 of organic)'
				cmd.select('inter', binding)
				cmd.color('red', 'inter')
				cmd.delete('inter')
				
				## remove ligand from translated chain
				cmd.remove(selection)	
					
	## delete state object
	cmd.delete(pdbName + "_"+ s)	

cmd.zoom('all')

print time.clock() - start_time, "seconds"
			