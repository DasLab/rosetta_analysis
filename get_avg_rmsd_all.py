import numpy as np
import os
from pymol import cmd

##### FOR analysis of DRRAFTER simulated benchmark ######
####### (VERY specific to this use case) #######
## run in pymol, because I was doing some additional visual analysis...

systems = ['1b7f',
	'1dfu',
	'1p6v',
        '1jbs',
        '1wpu',
        '1wsu',
        '2asb',
        '2bh2',
        '2qux',
        '3bx2']

native_pdbs = { '1b7f': '1b7f_bound.pdb',
		'1jbs': '1jbs_bound_orig.pdb',
		'1wpu': '1wpu_bound_NEW.pdb',
		'1wsu': '1wsu_bound.pdb',
		'2asb': '2asb_bound.pdb',
		'2bh2': '2bh2_bound_full.pdb',
		'2qux': '2qux_bound.pdb',
		'3bx2': '3bx2_bound.pdb',
		'1p6v': '1p6v_bound.pdb',
		'1dfu': '1dfu_bound.pdb' }

resolutions = ['3A', '5A', '7A']

start_dir = os.getcwd()

for sys in systems:
	base_name=sys
	for reso in resolutions:
		print sys, reso, '\n\n'
		os.chdir( start_dir )
		# change to the directory
		run_dir = sys
		if sys == '1jbs':
			run_dir += '/add_BP/'
		run_dir += '/%s/' %(reso)

		os.chdir( run_dir )

		rmsd_sele='resn A+U+G+C'
		
		for s1 in range(1,11):
			cmd.load('%s_fold_and_dock.out.%d.pdb' %(base_name,s1), 's%d' %(s1))
		# and load the native
		cmd.load('%s/native_pdbs/%s' %(start_dir, native_pdbs[sys]),'native')

		## get the RMSD to native
		#rmsds_to_native = []
		#for s1 in range(1,11):
		#	cmd.align('not resn A+U+G+C and s%d' %(s1), 'native')
		#	rms = cmd.rms_cur( '%s and s%d' %(rmsd_sele, s1), 'native')
		#	rmsds_to_native.append( rms )
		
		f = open('pairwise_rmsds_top_10.txt','w')
		pairwise_rmsds = []
		for s1 in range(1,11):
			for s2 in range(s1+1,11):
				# align based on the protein
				cmd.align( 'not resn A+U+G+C and s%d' %(s1), 'not resn A+U+G+C and s%d' %(s2) )
				# get only the RNA rmsd (without alignment)
				rms = cmd.rms_cur( '%s and s%d' %(rmsd_sele, s1), '%s and s%d' %(rmsd_sele, s2) )
				pairwise_rmsds.append( rms )
				print s1, s2, rms
				f.write('%d %d %0.3f\n' %(s1,s2,rms))
		
		pairwise_rmsds = np.array( pairwise_rmsds )
		avg_rmsd = np.mean( pairwise_rmsds )
		#rmsds_to_native = np.array( rmsds_to_native )
		#min_rmsd_to_native = np.min( rmsds_to_native )
		print "Mean rmsd ", avg_rmsd
		#print "Min rmsd to native ", min_rmsd_to_native
		f.write('Mean rmsd %0.3f\n' %(avg_rmsd))
		#f.write('Min rmsd to native %0.3f\n' %(min_rmsd_to_native))
		
		f.close()

		# close models
		cmd.delete('native')
		for s1 in range(1,11):
			cmd.delete('s%d' %(s1))
