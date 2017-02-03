

Einleitung:
	- Visualisierung von Proteinen
		- NMR, x-ray ???
	- Exploded Views in Technik/volume objects

Implementierung:
	- PyMOL:
		- Movies in PyMOL (easing function) ???
	- generell: pseudo code
			- Preparation: remove solvents, generate objects
			- labeling
			- coloring binding site, contact site (gray, individual, none)
			while isColliding
				(- seperate selection)
				- translate chains
				- translate ligands

	- 'canonical'
		- Technik - Bsp.
		- translation vector: transAxes (PCA)
	- 'com'
		- CENTER OF MASS
		- translation vector: translation factor

	- labeling:
		- Canonical: flush right/left
		- Com: circular

	- weitere Methoden:
		- relabel
		- reorient_explosion
		- renew_representation
		- show_labels / hide_labels
		(- remove_solvent)


Results/Discussion/Outlook:
	- Studie
