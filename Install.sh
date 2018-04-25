# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  # New files that I have added
  cp angle_3spn2_stacking.cpp ..
  cp angle_3spn2_stacking.h ..
  cp dihedral_3spn2.cpp ..
  cp dihedral_3spn2.h ..  
  cp base_pair.h ..
  cp base_pair.cpp ..
  cp pair_3spn2.cpp ..
  cp pair_3spn2.h ..

  # 3SPN.2C
  cp bond_list.cpp ..
  cp bond_list.h ..
  cp angle_list.cpp ..
  cp angle_list.h ..
  cp dihedral_list.cpp ..
  cp dihedral_list.h ..
  cp pair_list.cpp ..
  cp pair_list.h ..
  
  # Patch angle_hybrid.cpp
  cp ../angle_hybrid.cpp .
  patch angle_hybrid.cpp angle_hybrid.patch
  cp angle_hybrid.cpp ../

elif (test $1 = 0) then

  rm -f ../angle_3spn2_stacking.cpp 
  rm -f ../angle_3spn2_stacking.h 
  rm -f ../dihedral_3spn2.cpp
  rm -f ../dihedral_3spn2.h 
  rm -f ../base_pair.cpp
  rm -f ../base_pair.h
  rm -f ../pair_3spn2.cpp
  rm -f ../pair_3spn2.h

  # Files for 3SPN.2C
  rm -f ../bond_list.cpp
  rm -f ../bond_list.h
  rm -f ../angle_list.cpp
  rm -f ../angle_list.h
  rm -f ../dihedral_list.cpp
  rm -f ../dihedral_list.h
  rm -f ../pair_list.cpp
  rm -f ../pair_list.h

  # unpatching angle_hybrid.cpp
  cp ../angle_hybrid.cpp .
  patch -R angle_hybrid.cpp angle_hybrid.patch
  cp angle_hybrid.cpp ../

fi
