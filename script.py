#!/usr/bin/python
import sys
import math

def calculate_angle(p1, p2, p3, p4):

	b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
	b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
	b3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])

	s1 = (b1[1] * b2[2]) - (b1[2] * b2[1])
	s2 = (b1[2] * b2[0]) - (b1[0] * b2[2])
	s3 = (b1[0] * b2[1]) - (b1[1] * b2[0])
	n1 = (s1, s2, s3)

	s1 = (b2[1] * b3[2]) - (b2[2] * b3[1])
	s2 = (b2[2] * b3[0]) - (b2[0] * b3[2])
	s3 = (b2[0] * b3[1]) - (b2[1] * b3[0])
	n2 = (s1, s2, s3)

	mag_n1 = math.sqrt( (n1[0] * n1[0]) + (n1[1] * n1[1]) + (n1[2] * n1[2]) )
	n1 = (n1[0] / mag_n1, n1[1] / mag_n1, n1[2] / mag_n1)
	
	mag_n2 = math.sqrt( (n2[0] * n2[0]) + (n2[1] * n2[1]) + (n2[2] * n2[2]) )
	n2 = (n2[0] / mag_n2, n2[1] / mag_n2, n2[2] / mag_n2)

	mag_b2 = math.sqrt( (b2[0] * b2[0]) + (b2[1] * b2[1]) + (b2[2] * b2[2]) )
	b2 = (b2[0] / mag_b2, b2[1] / mag_b2, b2[2] / mag_b2)

	s1 = (n1[1] * b2[2]) - (n1[2] * b2[1])
	s2 = (n1[2] * b2[0]) - (n1[0] * b2[2])
	s3 = (n1[0] * b2[1]) - (n1[1] * b2[0])
	m1 = (s1, s2, s3)

	x = (n1[0] * n2[0]) + (n1[1] * n2[1]) + (n1[2] * n2[2])
	y = (m1[0] * n2[0]) + (m1[1] * n2[1]) + (m1[2] * n2[2])

	angle =  math.atan2(y, x) * 180 / math.pi

	return angle



if len(sys.argv) != 2:
	print "./script.py <pdbfilename>"
	sys.exit(0)

filename = sys.argv[1]
outputfilename = ""

if filename.endswith('.pdb'):
	outputfilename = filename[:-4]
	outputfilename += '_output.txt'

else:
	print "Pass a pdb file (having .pdb extension)"
	sys.exit(0)

pdbfile = open(sys.argv[1], 'rU')
output = open(outputfilename, "w")

total_unknowns = 0
angles = []
tuple_of_angles = (0, 0, 0)
index_of_angle = 0
n_1 = {}
n_2 = {}
total_amino_acids = 0
amino_acids = {}
ligands = {}
protein_chains = {}
current_chain = None
n_1['N'] = n_2['N'] = (float('nan'), float('nan'), float('nan'))
n_1['CA'] = n_2['CA'] = (float('nan'), float('nan'), float('nan'))
n_1['C'] = n_2['C'] = (float('nan'), float('nan'), float('nan'))
input_lines = pdbfile.readlines()
pdbfile.close()
i = -1

while i < len(input_lines) - 1:
	i += 1
	present_line = input_lines[i]
	if present_line.startswith('TITLE'):
		name = present_line[5:].strip()

	elif present_line.startswith('HET') and not present_line.startswith('HETATM'):
		input_list = present_line.split()
		ligand = input_list[1].strip()
		if ligand not in ligands:
			ligands[ligand] = 1

	elif present_line.startswith('SEQRES'):
		input_list = present_line.split()
		chain_name = input_list[2].strip()			
		if chain_name not in protein_chains:
			protein_chains[chain_name] = int(input_list[3])
			total_amino_acids += protein_chains[chain_name]
		for it in input_list[4:]:
			it = it.strip()
			if it not in amino_acids:
				amino_acids[it] = 0
			amino_acids[it] += 1
			
	elif present_line.startswith('ATOM'):

		input_list = present_line.split()
		if input_list[4].strip() == current_chain:

			if index_of_angle == 1:
				n_1['N'] = n_2['N']
				n_1['CA'] = n_2['CA']
				n_1['C'] = n_2['C']

				while True:
					if present_line.startswith('ATOM') == False or input_list[4].strip() != current_chain:
						angles.append(('psi', float('nan')))
						angles.append(('omega', float('nan')))
						index_of_angle = (index_of_angle + 2) % 3
						break

					if input_list[2].strip() == 'N':
						n_2['N'] = (float(input_list[6].strip()), float(input_list[7].strip()), float(input_list[8].strip()))
						psi = calculate_angle(n_1['N'], n_1['CA'], n_1['C'], n_2['N'])
						tuple_of_angles = ('psi', psi)
						angles.append(tuple_of_angles)
						index_of_angle = (index_of_angle + 1) % 3
						break
					i += 1
					present_line = input_lines[i]
					input_list = present_line.split()
				continue

			elif index_of_angle == 2:
				n_2['CA'] = (float(input_list[6].strip()), float(input_list[7].strip()), float(input_list[8].strip()))
				omega = calculate_angle(n_1['CA'], n_1['C'], n_2['N'], n_2['CA'])
				tuple_of_angles = ('omega', omega)
				angles.append(tuple_of_angles)
				index_of_angle = (index_of_angle + 1) % 3
				continue

			elif index_of_angle == 0:
				n_2['C'] = (float(input_list[6].strip()), float(input_list[7].strip()), float(input_list[8].strip()))
				phi = calculate_angle(n_1['C'], n_2['N'], n_2['CA'], n_2['C'])
				tuple_of_angles = ('phi', phi)
				angles.append(tuple_of_angles)
				index_of_angle = (index_of_angle + 1) % 3
				continue

		else:
			current_chain = input_list[4].strip()
			n_2['N'] = (float(input_list[6].strip()), float(input_list[7].strip()), float(input_list[8].strip()))
			i += 1
			present_line = input_lines[i]
			input_list = present_line.split()
			n_2['CA'] = (float(input_list[6].strip()), float(input_list[7].strip()), float(input_list[8].strip()))
			i += 1
			present_line = input_lines[i]
			input_list = present_line.split()
			n_2['C'] = (float(input_list[6].strip()), float(input_list[7].strip()), float(input_list[8].strip()))
			
			tuple_of_angles = ('phi', float('nan'))
			angles.append(tuple_of_angles)
			index_of_angle = (index_of_angle + 1) % 3
			


if 'UNK' in amino_acids:
	total_unknowns = amino_acids['UNK']
	del amino_acids['UNK']

output.write(name + '\n')

output.write('LENGTH\t'+  str(total_amino_acids) + '\n')

output.write('CHAINS\t' + str(len(protein_chains.keys())) + '\t' + ','.join(sorted(protein_chains.keys())) + '\n')

total_amino_acids = float(total_amino_acids)	
for k in sorted(amino_acids.keys()):
	output.write(k + '\t' + str(amino_acids[k] / total_amino_acids) + '\n')

output.write('UNKNOWN\t' + str(total_unknowns) + '\n')

output.write('LIGANDS\t' + ','.join(sorted(ligands.keys())) + '\n')

i = 0
for k in sorted(protein_chains.keys()):
	output.write('CHAIN-' + str(k) + '\n')
	while i < len(angles)-1:
		phi = angles[i][1]
		psi = angles[i+1][1]
		omega = angles[i+2][1]
		if math.isnan(phi):
			phi = 'NA'
		if math.isnan(psi):
			psi = 'NA'
		if math.isnan(omega):
			omega = 'NA'
		output.write(str(phi) + '\t' + str(psi) + '\t' + str(omega) + '\n')
		i += 3
		if psi == 'NA' and omega == 'NA':
			break

output.close()