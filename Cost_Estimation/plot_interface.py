#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Anonymous"
__version__ = "0.0"
__status__ = "Research"

import math
import os 
import matplotlib.pyplot as plt
import TreeHeuristics as TH
import CostFunctions as CF
from  matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import text
import Tools

import numpy as np
from collections import defaultdict

def draw_intersecting_lines(ax, n, log_M, max_depth, circ_desc='Query', leg_font_size = 11):
	"""
	Intersection of plots with
		1. Expected cost of classical enumeration
		2. Quasi-Sqrt(classical cost)
		3. Target security of Kyber
		4. Expected cost of Grover on AES
	:param ax: Target axes object from matplotlib.
	:param n: BKZ Blocksize
	:param log_M: Bound number of combined, randomized bases used in cost estimation
	:param max_depth: MaxDepth constraint.
	:param leg_font_size: Font size of legend.
	:return:
	"""
	treeH = TH.TreeHeuristics()
	costF = CF.CostFunctions(treeH, n=n, bound_coefficient=-1, nbases=2 ** log_M, pprime=1,
							 const_reps_QPE=-1, const_reps_W=-1, num_nodes_branching_dfs=-1, force_DF=-1)
	log_cost_classical = costF.log_enumeration_cost_classical()

	quasi_sqrt = (log_cost_classical + math.log(n, 2)) / 2

	z_coord_line_leg = 35
	if n == 873 and max_depth > 0:
		z_coord_line_leg = 1

	if n == 623 and max_depth != 96 and max_depth != 0:
		z_coord_line_leg = 1

	if circ_desc == 'Minimal' and (n == 623 or n == 873):
		z_coord_line_leg = 1

	# QUASI SQRT ----------------------------------------
	#ax.axhline(y=quasi_sqrt, dashes=[4, 8], color='black', alpha=1, lw=1, label='Quasi-Sqrt(classical cost)'
	ax.axhline(y=quasi_sqrt, dashes=[4, 8], color='black', alpha=1, lw=1) # label='Quasi-Sqrt(classical cost)'

	line_leg_vert_align = 'bottom'
	if n == 623 and max_depth == 40:
		line_leg_vert_align = 'top'
	if circ_desc=='Query' and n == 623 and max_depth == 96:
		z_coord_line_leg = 40
	if n == 406 and max_depth == 64:
		line_leg_vert_align = 'top'

	text(z_coord_line_leg, quasi_sqrt, "Quasi-Sqrt(classical cost)", horizontalalignment='left', verticalalignment=line_leg_vert_align, fontsize=leg_font_size)

	# EXPECTED GROVER COST ----------------------------------------
	#ax.axhline(y=Tools.aes_expected_security[n][max_depth], dashes=[3, 3], color='#B90E0A', alpha=1, lw=1, label='Expected cExpected cost of Grover on AES')
	ax.axhline(y=Tools.aes_expected_security[n][max_depth], dashes=[3, 3], color='#B90E0A', alpha=1, lw=1) #, label='Expected cost of Grover on AES'

	line_leg_vert_align = 'bottom'
	if n == 623 and max_depth == 64:
		line_leg_vert_align = 'top'
	if n == 873 and max_depth > 0:
		line_leg_vert_align = 'top'

	text(z_coord_line_leg, Tools.aes_expected_security[n][max_depth], "Expected cost of Grover on AES", horizontalalignment='left', verticalalignment=line_leg_vert_align, fontsize=leg_font_size, color='#B90E0A')

	# TARGET SECURITY ----------------------------------------
	# ax.axhline(y=Tools.kyber_expected_security[n], dashes=[1, 1], color='#B90E0A', alpha=1, lw=1, label='Target security of Kyber')
	ax.axhline(y=Tools.kyber_expected_security[n], dashes=[1, 1], color='#B90E0A', alpha=1, lw=1) #, label='Target security of Kyber')
	line_leg_vert_align = 'bottom'
	if circ_desc=='Query' and n == 623 and max_depth == 96:
		z_coord_line_leg = 43

	if n == 873 and (max_depth == 64 or max_depth == 96 or max_depth == 0):
		line_leg_vert_align = 'top'
	if circ_desc=='Query' and n == 873:
		z_coord_line_leg = 35
	text(z_coord_line_leg, Tools.kyber_expected_security[n], "Target security of Kyber", horizontalalignment='left', verticalalignment=line_leg_vert_align, color='#B90E0A', fontsize=leg_font_size)

	# EXPECTED CLASSICAL ----------------------------------------
	#ax.axhline(y=log_cost_classical, dashes=[2, 6], color='black', alpha=1, lw=1, label='Expected cost of classical enumeration')
	ax.axhline(y=log_cost_classical, dashes=[2, 6], color='black', alpha=1, lw=1)#, label='Expected cost of classical enumeration')
	text(25, log_cost_classical-1, "Expected cost of classical enumeration", horizontalalignment='left', verticalalignment='top', fontsize=leg_font_size)



def setup_plot(ax, x_coords, y_coords, md, n, log_M, axis_fontsize=14):
	"""
	Plot figure preparation.
	:param ax:	Target axes object from matplotlib
	:param x_coords: List of x_coords in plot.
	:param y_coords: List of y_coords in plot.
	:param md: Maxdepth constraint on circuit depth/
	:param n: BKZ blocksize
	:param log_M: Bound on number of combined bases used in cost estimation
	:param axis_fontsize: Fontsize for axis labels.
	:return:
	"""
	# X Axis
	x_upper = math.ceil(x_coords[-1]) + 5
	major_xticks = range(0, x_upper, 10)
	minor_xticks = range(0, x_upper, 2)
	plt.xlim(0, x_upper)
	ax.set_xticks(major_xticks)
	ax.set_xticks(minor_xticks, minor=True)
	ax.tick_params(labelsize=axis_fontsize)

	# Y Axis
	treeH = TH.TreeHeuristics()
	log_N_kh_M = treeH.log_avg_N_kh_M(k=0, h=n + 1, n=n, nbases=2 ** log_M, pprime=1)
	quasi_sqrt = (log_N_kh_M + math.log(n, 2)) / 2
	y_upper = math.ceil(max(y_coords + [quasi_sqrt, Tools.kyber_expected_security[n], Tools.aes_expected_security[n][md], log_N_kh_M])) + 10

	if md == 0:
		y_lower = int(min(y_coords + [Tools.kyber_expected_security[n], quasi_sqrt, Tools.aes_expected_security[n][md]])) - 1
	else:
		y_lower = 0

	major_yticks = range(max(y_lower, 0), y_upper, 50)
	minor_yticks = range(max(y_lower, 0), y_upper, 10)
	plt.ylim(max(y_lower, 0), y_upper)
	ax.set_yticks(major_yticks)
	ax.set_yticks(minor_yticks, minor=True)

	ax.grid(which='major', alpha=0.5)
	plt.xlabel("Log(Jensen's Gap) z", fontsize=axis_fontsize)
	plt.ylabel("Log", fontsize=axis_fontsize)

def anotate_critical_values(ax, n, max_depth, log_M, z_to_cost, label_fontsize = 12):
	"""
	Anotation of critical intersections.
	:param ax: Axes object from matplotlib
	:param n: BKZ blocksize
	:param max_depth: MaxDepth constraint for depth of quantum circuit
	:param log_M: Bound on number of randomized, combined bases used in cost estimation
	:param z_to_cost: Dictionary [Jenzen's z] --> cost
	:param label_fontsize: Font size of labels.
	:return:
	"""
	def anotate(ax, max_depth, z_value, z_to_cost, inequality = '=', x_offset=0, y_offset=0):
		"""
		Annotation subfunction for all intersections.
		"""
		x_coord = z_value
		y_coord = z_to_cost[z_value].log_gcost_total()

		aLabel = 'z=' + str(z_value) + '\n' + 'Cost' + inequality + str(math.ceil(y_coord))
		if max_depth > 0:
			aLabel += '\n' + 'k=' + str(z_to_cost[z_value].k) + '\n' + 'y=' + str(z_to_cost[z_value].log_y)

		ax.plot(x_coord, y_coord, marker="o", markersize=5, color='black', alpha=.8)
		anot = ax.annotate(aLabel, xy=(x_coord, y_coord), xytext=(x_coord + x_offset, y_coord + y_offset),
						   color='black', ha='center', va='center', fontsize=label_fontsize, alpha=1,
						   arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0", linewidth=0.5))
		anot.set_bbox(dict(alpha=.6, linewidth=0, edgecolor='white', facecolor='white'))
		#x_coords = list(sorted(z_to_cost.keys()))

	min_z = min(z_to_cost.keys())
	y_offset = max(35, max([z_to_cost[min_z].log_gcost_total()]) / 4)

	# First value, potentially z = 0
	anotate(ax, max_depth, min_z, z_to_cost, x_offset=1, y_offset=y_offset)

	# Intersections
	intersections = Tools.get_intersections(n, max_depth, log_M, z_to_cost)

	anotated_z_values = [min_z]
	z_to_ineq = {}

	for k,v in intersections.items():
		if v[0] not in anotated_z_values:
			z_to_ineq[v[0]] = v[1]
			anotated_z_values.append(v[0])

	second_last_z = 0
	x_offset = 0
	for i,z in enumerate(sorted(z_to_ineq.keys())):
		if i % 2 == 1:
			if z - second_last_z < 10:
				x_offset = 4
			else:
				x_offset = 0
			second_last_z = z
		else:
			x_offset = 0

		last_z = z
		if z_to_ineq[z] != '>>':
			# z value, s.t. total <= cost
			y_offset *= -1
			anotate(ax, max_depth, z, z_to_cost, inequality='$\leq$', x_offset=x_offset, y_offset=y_offset)


def plot_cost_estimation_z_to_g_cost(plot_dir, filename, z_to_cost, max_depth, n, log_M, circuit, prefix, postfix, anotate=True, leg_font_size = 12):
	"""

	:param plot_dir: Target directory to save plot to
	:param filename: Target filename for plot
	:param z_to_cost: Dictionary [Jensen's z] --> cost
	:param max_depth: maxDepth constraint on quantum circuit depth as used in cost estimation
	:param n: BKZ blocksize as used in cost estimation
	:param log_M: Bound on number of combined, randomized bases as used in cost estimation
	:param leg_font_size: Legend font size
	:return:
	"""
	fig, ax = plt.subplots()

	# Individual Costs
	x_coords = sorted(list(z_to_cost.keys()))
	y_coords = []
	y_coords_quantum = []
	y_coords_classical = []
	qaracm_values = []
	log_y_values = []
	k_values = []

	for z in x_coords:
		v = z_to_cost[z]
		y_coords.append(float(v.log_gcost_total()))
		y_coords_quantum.append(float(v.log_gcost_quantum))
		y_coords_classical.append(float(v.log_gcost_classical))
		qaracm_values.append(float(v.qracm))
		log_y_values.append(v.log_y)
		k_values.append(v.k)

	# Setup axis, grid and labels
	setup_plot(ax, x_coords, y_coords_quantum, max_depth, n, log_M)

	# Draw aes-grover, target, quasi-quadratic lines
	draw_intersecting_lines(ax, n, log_M, max_depth, circ_desc=circuit.desc)

	# Draw Costs
	if max_depth > 0:
		ax.plot(x_coords, y_coords_quantum, dashes=[24, 2], color='#0080ff', label=f'Quantum GCost', alpha=0.8)
		ax.plot(x_coords, y_coords, color='#669900', label=f'Total GCost', alpha=0.5)
		ax.plot(x_coords, y_coords_classical, dashes=[24, 2], color='#8D4004', label=f'Classical GCost', alpha=0.8)
		ax.plot(x_coords, qaracm_values, linestyle='dashed', color='#55CEFF', label=f'QRACM', linewidth=0.8, alpha=0.8)
	else:
		#ax.plot(x_coords, y_coords_classical, linestyle='dashed', color='#8D4004', label=f'Expected cost of classical enumeration', alpha=0.8)
		ax.plot(x_coords, y_coords_quantum, dashes=[24, 2], color='#0080ff', label=f'Expected cost of quantum enumeration', alpha=0.8)

	# Find and mark intersections
	anotate_critical_values(ax, n, max_depth, log_M, z_to_cost)

	# # Put a legend to the right of the current axis
	#ax.legend(loc='center left', bbox_to_anchor=(1.04, 0.5))   # (horizontal, vertical)
	#ax.legend(loc='upper center', bbox_to_anchor=(.67, 1.15), fontsize=leg_font_size, framealpha=0.5)  # (horizontal, vertical)

	if max_depth == 0 and (n == 873 or n == 623):
		ax.legend(loc='lower left', bbox_to_anchor=(.19, .6), fontsize=leg_font_size, framealpha=0.8)  # (horizontal, vertical)
	elif n == 623:
		ax.legend(loc='lower left', bbox_to_anchor=(0, 0), fontsize=leg_font_size, framealpha=0.8)  # (horizontal, vertical)
	else:
		ax.legend(loc='lower left', bbox_to_anchor=(.19, 0), fontsize=leg_font_size, framealpha=0.8)  # (horizontal, vertical)

	plt.savefig(f"{plot_dir}/{filename}.pdf", dpi=600, transparent=False, bbox_inches='tight')
	plt.close()
	plt.cla()
	plt.clf()

def plot_explanation_contour_y(plot_dir, circuit, y_z_to_cost, n, prefix, postfix):
	"""
		Contour plot for explanation.
	"""
	full_plot_dir = plot_dir
	os.system('mkdir -p ' + full_plot_dir)
	os.system('touch ' + full_plot_dir + '/PLACEHOLDER')

	Z_1 = []
	Z_2 = []
	X = sorted(y_z_to_cost[next(iter(y_z_to_cost))].keys())
	Y = sorted(y_z_to_cost.keys())
	K = []

	for i, (y, z_to_cost) in enumerate(y_z_to_cost.items()):
		Z_1.append([])
		Z_2.append([])
		K.append([])

		print(f"Y {y}")

		for z, cost in z_to_cost.items():
			log_total = cost.log_gcost_total()
			Z_1[-1].append(log_total)
			Z_2[-1].append(cost.k)
			K[-1].append(cost.k)

			print(f"  z={z} --> log_total={log_total}, level k={cost.k}")

	def prepPlot(X, Y, Z, Z_label, full_plot_dir, filename):
		fig, ax = plt.subplots()

		plt.xlabel("Log(Jensen's Gap) z")
		plt.ylabel("Log(Combined Nodes) y)")
		cs = ax.contourf(X, Y, Z, cmap='magma_r', levels=30) # RdYlGn_r, nice but not CB friendly
		cbar = plt.colorbar(cs)
		cbar.set_label(Z_label, rotation=270, labelpad=10)

		fig.savefig(f"{full_plot_dir}/{filename}", dpi=600, transparent=False, bbox_inches='tight')
		plt.cla()
		plt.clf()

	# Plot Z_1: Total GCost
	Z_label = 'Log $\mathbb{E}$(GCost)'
	strings = [prefix, 'Y-to_Z', str(n), postfix, circuit.desc]
	filename = f"{'_'.join(filter(None, strings))}.pdf"
	prepPlot(X, Y, Z_1, Z_label, full_plot_dir, filename)

	# Plot Z_1: Gap Quantum to Classical
	Z_label = 'Log($\mathbb{E}$(Quantum GCost) - $\mathbb{E}$(Classical GCost))'
	strings = [prefix, 'Gap', str(n), postfix, circuit.desc]
	filename = f"{'_'.join(filter(None, strings))}.pdf"
	prepPlot(X, Y, Z_2, Z_label, full_plot_dir, filename)

	# Plot Z_1: Gap Quantum to Classical
	Z_label = 'Level k'
	strings = [prefix, 'k', str(n), postfix, circuit.desc]
	filename = f"{'_'.join(filter(None, strings))}.pdf"
	prepPlot(X, Y, K, Z_label, full_plot_dir, filename)

