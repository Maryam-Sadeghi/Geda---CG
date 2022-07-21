#!/usr/bin/env python3

# The code calculates the angle between two planes inlcuding phenyl rings and evaluates the planiraty of the planes
import MDAnalysis
import numpy.linalg
from numpy.linalg import norm

u = MDAnalysis.Universe("md-anneal2.gro","CG.xtc") 
print(u)


resnames = u.atoms.resnames
resids = u.atoms.resids
residues = list(zip(resnames, resids))
#print(residues)

molecule_flatness = {}

nframes = 0
nres = 0

### initialize hashes
for res in residues:
  if res[0] == "M00":
    molecule_flatness[res[1]] = 0
    nres = nres + 1
#    print('nres ', nres)

angles_t = open('angles_dt.dat','w')

#### single residue properties
for ts in u.trajectory:
  angles_t.write(str(u.trajectory.time)+" ")
  for res in residues:
    if res[0] == "M00":

### angles between two rings
      mol_v1_p1 = u.select_atoms("resname M00 and name B5")
      mol_v1_p2 = u.select_atoms("resname M00 and name B6")
      mol_v1_p3 = u.select_atoms("resname M00 and name B7")
      mol_v2_p1 = u.select_atoms("resname M00 and name B9")
      mol_v2_p2 = u.select_atoms("resname M00 and name B10")
      mol_v2_p3 = u.select_atoms("resname M00 and name B11")

      v1 = mol_v1_p1.positions[0] - mol_v1_p2.positions[0]
      v2 = mol_v1_p1.positions[0] - mol_v1_p3.positions[0]
      v3 = mol_v2_p1.positions[0] - mol_v2_p2.positions[0]
      v4 = mol_v2_p1.positions[0] - mol_v2_p3.positions[0]

      n1 = numpy.cross(v1,v2)
      n2 = numpy.cross(v3,v4)

      angle = numpy.arccos(numpy.dot(n1,n2)/(norm(n1)*norm(n2)))

      angle = numpy.rad2deg(angle)
      molecule_flatness[res[1]] = molecule_flatness[res[1]] + angle

      angles_t.write(str(angle)+" ")

  angles_t.write('\n')
  nframes+=1

print(nframes)

angles_t.close()

##### print results
#angles_ave = open('angles_ave.dat','w')

#for res in residues:
#  molecule_flatness[res] = molecule_flatness[res]/nframes;
#  angles_ave.write(str(res)+" "+str(molecule_flatness[res]))
#  angles_ave.write('\n')

#angles_ave.close()
