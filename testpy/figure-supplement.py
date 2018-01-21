#! /usr/bin/env/python
#coding: UTF-8
import re
import os
import numpy as np

fileName= "est_chain_60.txt"
const_cut_dist = 10
len_np = 2.5

def importfile(fileName):
	'''
		the imported file is the chain information recording file;
		which recorded the number id of chains in a box of length of len;
		the first line is the parameters of this System(grafted homo- and AB co- polymer)
		the parameters are'"(grafted-space,len-grafted-HP,len-A-block,len-B-block,number-AB-Cp)";
	'''
	fp = open(fileName,'r')
	try:
		tpl_para = tuple(map(eval,fp.readline().strip().split()))
		lst_graftedHomopolymer = []
		lst_ABcopolymer = []
		lst_chains = []
		lst_singleChain = []
		for chain in fp.readlines():
			if chain !='\n':       #if the line of chains file is not the blank line; 
				lst_singleChain = list(map(eval,chain.strip('\n').split()))
				lst_singleChain[0] = lst_singleChain[0] - 1 # subtract 1 for id of chain;
				lst_chains.append(lst_singleChain)
		# print(lst_chains)
	finally:
		fp.close()		#finally we close the import file;
	print(tpl_para)
	return tpl_para,lst_chains

def exportfile():
	pass
def density_cal(monomer,tpl_para,lst_atomInfo,tpl_neighPoint):
	'''
		calculate the density of monomer;
	'''
	len_simuBox = tpl_para[0]
	print(len(tpl_neighPoint[0]))
	typeInitial = 0
	lst_density =[]
	if monomer == 'A' or monomer == 'a':
		typeSelected = 1
	elif monomer == 'B' or monomer == 'b':
		typeSelected = 2
	else:
		typeSelected = 1
		typeInitial = 2
	for atom in lst_atomInfo:			#[0, 0, 0, 0]
		density = 0.0
		id = atom[0]*len_simuBox**2 + atom[1]*len_simuBox +atom[2]
		for order_np in tpl_neighPoint[id]:
			if lst_atomInfo[order_np][-1] == typeSelected or lst_atomInfo[order_np][-1] == typeInitial:
				density = density + 1
		lst_density.append(density/len(tpl_neighPoint[0]))
	# print(tpl_neighPoint[id])
	return lst_density
def box_building(len):
	"""
		the function make that building a simulation box
		and return the API : lst_atomInfo which has record the information of all of atoms;
	"""
	number_totalAtoms = len**3
	lst_atomInfo = []
	singleAtom = [0,0,0,0]
	for x in range(0,len):
		for y in range(0,len):
			for z in range(0,len):
				singleAtom = [x,y,z,0]
				lst_atomInfo.append(singleAtom)
	return lst_atomInfo

def Property_assignment(tpl_para,lst_atomInfo,lst_chains):
	len_simuBox = tpl_para[0]
	num_grafteHP = len_simuBox**2/(tpl_para[1]**2)
	len_HP,len_ASegment,len_BSegment,num_Copolymer =tpl_para[2:]
	for chain in lst_chains:
		if len(chain)<= len_HP+1:
			for atom in chain[1:]:
				lst_atomInfo[atom][-1] = 1
		else:
			for atom in chain[1:len_ASegment+1]:
				lst_atomInfo[atom][-1] = 1
			for atom in chain[-len_BSegment:]:
				lst_atomInfo[atom][-1] = 2
	return lst_atomInfo

def neighbourPoint_building(tpl_para,lst_atomInfo):
	'''
		build the neighbourPoint tuple, and it is would not change;
	'''
	len_simuBox = tpl_para[0]
	lst_NPcoorPlus = []
	tpl_atomInfo = tuple(lst_atomInfo)
	iternum = 0 
	for x in range(-1,2):
		for y in range(-1,2):
			for z in range(-1,2):
				nearCoor_plus = [x,y,z]
				if sum(pow(item,2) for item in nearCoor_plus)<len_np and \
				sum(pow(item,2) for item in nearCoor_plus)>0:
					iternum = iternum + 1
					lst_NPcoorPlus.append(nearCoor_plus)
	lst_NPcoorPlus = tuple(lst_NPcoorPlus)
	lst_neighPoint = []
	if iternum != len(lst_NPcoorPlus):
		print("error in building the list:lst_NPcoorPlus")
	number_totalAtoms = len(tpl_atomInfo)
	for atom in tpl_atomInfo:
		coor_neighP =[]
		for nearCoor_plus in lst_NPcoorPlus:
			x = atom[0] + nearCoor_plus[0]
			y = atom[1] + nearCoor_plus[1]
			z = atom[2] + nearCoor_plus[2]
			x = periodBonCond(x,len_simuBox)
			y = periodBonCond(y,len_simuBox)
			z = periodBonCond(z,len_simuBox)
			id  = x*len_simuBox**2 + y *len_simuBox + z
			if id >number_totalAtoms:
				exit(1)
			coor_neighP.append(id)
		lst_neighPoint.append(coor_neighP)
	lst_neighPoint = tuple(lst_neighPoint)
	return lst_neighPoint

def periodBonCond(value,periodlen):
	if value >=periodlen:
		value = value - periodlen
	elif value < 0:
		value = value + periodlen
	return value

def export2MayaviVTK(monomerName,monomerDensity,tpl_para):
	'''
		export monomer to density profile figure: MayaviVTK file;
	'''
	len_simuBox = tpl_para[0]
	fp = open(monomerName+r'density.vtk','w')
	try:
		fp.write("# vtk DataFile Version 2.0\n")
		fp.write("CT Cylinder interface\n")
		fp.write("ASCII\n\n")
		fp.write("DATASET STRUCTURED_POINTS\n")
		fp.write("DIMENSIONS\t%d\t%d\t%d\n" %(len_simuBox,len_simuBox,len_simuBox))
		fp.write("ORIGIN 0.000000 0.000000 0.000000\n")
		fp.write("SPACING 0.100000 0.100000 0.100000\n\n")
		fp.write("POINT_DATA\t%d\n" %(len_simuBox**3))
		fp.write("SCALARS scalars float\n\n")
		fp.write("LOOKUP_TABLE default\n")
	finally:
		fp.close()
	pass

def main():
	tpl_parameters,lst_chains = importfile(fileName)
	lst_atomInfo = box_building(tpl_parameters[0])
	tpl_neighPoint = neighbourPoint_building(tpl_parameters,lst_atomInfo)
	exportfile()
	lst_atomInfo = Property_assignment(tpl_parameters,lst_atomInfo,lst_chains)
	lst_Adensity = density_cal('A',tpl_parameters,lst_atomInfo,tpl_neighPoint)
	export2MayaviVTK('A',lst_Adensity,tpl_parameters)
	
main()