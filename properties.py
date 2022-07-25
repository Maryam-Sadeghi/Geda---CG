import MDAnalysis
import numpy.linalg

def normalize(v):
  norm=numpy.linalg.norm(v)
  if norm < 0.001:  
    print ("very small norm")
    sys.exit(1)
  else:
    return v/norm 


u = MDAnalysis.Universe("nvt.gro","traj_im10.xtc")  # always start with a Universe
print(u)


resnames = u.atoms.residues.resnames
resids = u.atoms.residues.resids
residues = zip(resnames, resids)
#print(residues)
#radius_of_gyration();

rads_of_gyration = {}
end_to_end_dist = {}
molecule_flatness = {}

nframes = 0
nres = 0

### initialize hashes
for res in residues:
  if res[0] == "M00":
    rads_of_gyration[res[1]] = 0
    end_to_end_dist[res[1]] = 0
    molecule_flatness[res[1]] = 0
    nres = nres + 1

distance_matrix_firstHalf = numpy.zeros((nres,nres))
orientation_matrix_firstHalf = numpy.zeros((nres,nres))
distance_matrix = numpy.zeros((nres,nres))
orientation_matrix = numpy.zeros((nres,nres))


angles_t = open('angles_dt.dat','w')
ete_t = open('ete_dt.dat','w')
radgyr_t = open('radgyr_dt.dat','w')

#### single residue properties
for ts in u.trajectory:
  angles_t.write(str(u.trajectory.time)+" ")
  ete_t.write(str(u.trajectory.time)+" ")
  radgyr_t.write(str(u.trajectory.time)+" ")
  for res in residues:
    if res[0] == "M00":
### end to end distance
      mol_end_1 = u.select_atoms("resid "+str(res[1])+" and name N1")
      mol_end_2 = u.select_atoms("resid "+str(res[1])+" and name N7")
      r = mol_end_1.positions[0] - mol_end_2.positions[0] ### end to end vector
      d = numpy.linalg.norm(r) ## distance
      end_to_end_dist[res[1]] = end_to_end_dist[res[1]] + d
      ete_t.write(str(d)+" ")
### radius of gyration
      mol = u.select_atoms("resid "+str(res[1]))
      radGyr = mol.radius_of_gyration();
      rads_of_gyration[res[1]] = rads_of_gyration[res[1]] + radGyr 
      radgyr_t.write(str(radGyr)+" ")
### co-planarity of rings
      mol_v1_p1 = u.select_atoms("resid "+str(res[1])+" and name N4")
      mol_v1_p2 = u.select_atoms("resid "+str(res[1])+" and name C14")
      mol_v2_p1 = mol_v1_p1
      mol_v2_p2 = u.select_atoms("resid "+str(res[1])+" and name N3")
      mol_v3_p1 = mol_v2_p2 
      mol_v3_p2 = u.select_atoms("resid "+str(res[1])+" and name C6")
      v1 = mol_v1_p1.positions[0] - mol_v1_p2.positions[0]
      v2 = mol_v2_p1.positions[0] - mol_v2_p2.positions[0]
      v3 = mol_v3_p2.positions[0] - mol_v3_p1.positions[0]

      n1 = numpy.cross(v1,v2)
      n2 = numpy.cross(v2,v3)
  
      n1 = normalize(n1)
      n2 = normalize(n2)
      v2 = normalize(v2)
      m1 = numpy.cross(n1,v2)

      x1 = numpy.dot(n1,n2)
      y1 = numpy.dot(m1,n2)
      angle = numpy.arctan2(x1,y1)
      angle = numpy.rad2deg(angle)
      molecule_flatness[res[1]] = molecule_flatness[res[1]] + angle

      angles_t.write(str(angle)+" ")

  ete_t.write('\n')  
  radgyr_t.write('\n')
  angles_t.write('\n')
  nframes+=1

print(nframes)

ete_t.close()
radgyr_t.close()
angles_t.close()

halfwayframe = nframes/2
print(halfwayframe)
#### matrix properties
for ts in u.trajectory:
  for res1 in residues:
    if res1[0] == "M00":
      mol1 = u.select_atoms("resid "+str(res1[1]))
      centroid_mol1 = mol1.centroid()
      mol1_p1 = u.select_atoms("resid "+str(res1[1])+" and name N4")
      mol1_p2 = u.select_atoms("resid "+str(res1[1])+" and name N7")
      v1 = mol1_p1.positions[0] - mol1_p2.positions[0]
      v1 = normalize(v1)
      for res2 in residues:
      	if res1[1] != res2[1]:
         if res2[0] == "M00":
##### average distances
            mol2 = u.select_atoms("resid "+str(res2[1]))
            centroid_mol2 = mol2.centroid()
            r = centroid_mol1 - centroid_mol2
            d = numpy.linalg.norm(r) 
            index1 = res1[1] - 1
            index2 = res2[1] - 1
            distance_matrix[index1,index2] = distance_matrix[index1,index2] + d
            if ts.frame < halfwayframe:
              distance_matrix_firstHalf[index1,index2] = distance_matrix_firstHalf[index1,index2] + d
###### avaerage ring orientations
            mol2_p1 = u.select_atoms("resid "+str(res2[1])+" and name N4")
            mol2_p2 = u.select_atoms("resid "+str(res2[1])+" and name N7")
            v2 = mol2_p1.positions[0] - mol2_p2.positions[0]
            v2 = normalize(v2)
            costheta = numpy.dot(v1,v2)
            orientation_matrix[index1,index2] = orientation_matrix[index1,index2] + costheta
            if ts.frame < halfwayframe:
            	orientation_matrix_firstHalf[index1,index2] = orientation_matrix_firstHalf[index1,index2] + costheta

distance_matrix = distance_matrix/nframes
distance_matrix_firstHalf = distance_matrix_firstHalf/halfwayframe
orientation_matrix = orientation_matrix/nframes
orientation_matrix_firstHalf = orientation_matrix_firstHalf/halfwayframe

distsFile = open('ave_dists_for_gnuPlot.dat','w')
distsFile_firstHalf = open('ave_distsFirstHalf_for_gnuPlot.dat','w')
oriFile = open('ave_ori_for_gnuPlot.dat','w')
oriFile_firstHalf = open('ave_oriFirstHalf_for_gnuPlot.dat','w')

for (i,j), value in numpy.ndenumerate(distance_matrix):
   s = str(i)+" "+str(j)+" "+str(distance_matrix[i,j])
   distsFile.write(s)
   distsFile.write('\n')

distsFile.close()

for (i,j), value in numpy.ndenumerate(distance_matrix_firstHalf):
   s = str(i)+" "+str(j)+" "+str(distance_matrix_firstHalf[i,j])
   distsFile_firstHalf.write(s)
   distsFile_firstHalf.write('\n')

distsFile_firstHalf.close()

for (i,j), value in numpy.ndenumerate(orientation_matrix):
   s = str(i)+" "+str(j)+" "+str(orientation_matrix[i,j])
   oriFile.write(s)
   oriFile.write('\n')

oriFile.close()

for (i,j), value in numpy.ndenumerate(orientation_matrix_firstHalf):
   s = str(i)+" "+str(j)+" "+str(orientation_matrix_firstHalf[i,j])
   oriFile_firstHalf.write(s)
   oriFile_firstHalf.write('\n')

oriFile_firstHalf.close()


##### print results
angles_ave = open('angles_ave.dat','w')
ete_ave = open('ete_ave.dat','w')
radgyr_ave = open('radgyr_ave.dat','w')

for res in rads_of_gyration:
  rads_of_gyration[res] = rads_of_gyration[res]/nframes;
  end_to_end_dist[res] = end_to_end_dist[res]/nframes;
  molecule_flatness[res] = molecule_flatness[res]/nframes;
  angles_ave.write(str(res)+" "+str(molecule_flatness[res]))
  ete_ave.write(str(res)+" "+str(end_to_end_dist[res]))
  radgyr_ave.write(str(res)+" "+str(rads_of_gyration[res]))
  angles_ave.write('\n')
  ete_ave.write('\n')
  radgyr_ave.write('\n')


angles_ave.close()
ete_ave.close()
radgyr_ave.close()



