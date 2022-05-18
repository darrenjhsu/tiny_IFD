#!/bin/bash

PDB_FILES=$1    # string w/ regex that points to the pdb files of ligands to parameterized
for pdb in "${PDB_FILES[@]}"
do
	temp=${pdb##*/}		# removing any potential directories of a local or global file path; we only want the molecule name for future file naming
	temp2=${temp%.*}	# removing any potential file extension from the remaining string
	echo $temp2		# printing ligand information currently being used to create respective directory and file naming
	echo Adding hydrogens w/ obabel
	obabel -i pdb $pdb -o mol2 -O "$temp2"_H.mol2 -h	# read in the ligand's pdb (missing nonpolar hydrogens) file and output a mol2 file w/ hydrogens
	echo Adding atom numbering to atom_name
	python3 ../../../01_Workflow/utilities/fix_atom_numbering.py "$temp2"_H.mol2		# fix atom names by adding numbers to each atom name; no potential for a bug because mol2 file are fairly flexible
	mv output.mol2 "$temp2"_H.mol2
	obabel -i mol2 "$temp2"_H.mol2 -o pdb -O ${temp2%.*}_H.pdb	# generate the ligand's pdb with hydrogens file from the renamed mol2 file; this is for reference complex usage in the 06_Analysis
	echo Creating charges for the molecule at GAFF2/AM1-BCC lvl of theory
	antechamber -i "$temp2"_H.mol2 -fi mol2 -o "$temp2"_H_bcc.mol2 -fo mol2 -at gaff2 -c bcc -s 2	# read the ligand mol2 file into antechamber to assign gaff2 atom types and calculate AM1-BCC atomic point charges
	echo Checking parameters w/ parmchk2
	parmchk2 -i "$temp2"_H_bcc.mol2 -f mol2 -o "$temp2"_H_bcc.frcmod -s 2	# check the mol2 file for the potential for missing parameters; if any are found, fill in the gaps w/ other parameters that  are potentially adequate. CHECK THIS FILE FOR HIGH PENALTY SCORES
	echo Creating the prmtop file
	sed -e 's,aaa,'"$temp2"'_H_bcc,g' < ../../../01_Workflow/utilities/leap_input.in > temp.in	# sed replace some search strings and output a temporary tleap input file
	tleap -f temp.in #> "$temp2"_tl.out	# read the molecule into tleap to prepare a prmtop file and simulation ready pdb
	# clean up the working directory
	mkdir $temp2
	mv $temp2* $temp2/
	mv leap.log $temp2/
	mv ANTECHAMBER* $temp2/ 
	mv ATOMTYPE* $temp2/ 
	mv sqm* $temp2/ 
	# pass to amber mdgx or openmm script to run simulations
done

