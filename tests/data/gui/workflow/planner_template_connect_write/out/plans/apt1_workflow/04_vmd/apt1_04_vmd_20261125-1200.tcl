# VMD script written by MDBridge <VERSION> for project apt1, step 04
# run: vmd -e apt1_04_vmd_20261125-1200.tcl

set mol [mol new {../01_amber-build/apt1_01_tleap_20261125-1200.pdb} type pdb waitfor all]

mol delrep 0 $mol
display resize 1280 960
display projection Orthographic
display depthcue off
color Display Background white
axes location Off

color change rgb 23 0.165 0.471 0.839 ;# Aptamer
color change rgb 24 0.922 0.408 0.204 ;# Target

mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 23
mol selection {resid 1 to 12}
mol material AOChalky
mol addrep $mol
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 24
mol selection {resid 13 to 24}
mol material AOChalky
mol addrep $mol
mol representation Licorice 0.200000 12.000000 12.000000
mol color ColorID 23
mol selection {(resid 1 to 12) and nucleic and not backbone and noh}
mol material AOChalky
mol addrep $mol
mol representation Licorice 0.200000 12.000000 12.000000
mol color ColorID 24
mol selection {(resid 13 to 24) and nucleic and not backbone and noh}
mol material AOChalky
mol addrep $mol
mol representation VDW 0.600000 12.000000
mol color ColorID 19
mol selection {resname MG Mg MG2}
mol material AOChalky
mol addrep $mol

display resetview
