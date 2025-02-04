#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research"

import math
import matplotlib.pyplot as plt
import TreeHeuristics as TH
import CostFunctions as CF
import Tools

import configSecurity

def draw_intersecting_lines(ax, n, log_M, max_depth, constants, bound, leg_font_size = 11):
	"""
	Support fucntion to draw intersection of cost with
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
	costF = CF.CostFunctions(n=n, nbases=2 ** log_M, pprime=1, constants=constants, bound=bound)
	log_cost_classical = costF.log_enumeration_cost_classical()

	quasi_sqrt = (log_cost_classical + math.log(n, 2)) / 2

	# QUASI SQRT
	ax.axhline(y=quasi_sqrt, dashes=[4, 8], color='black', alpha=1, lw=1, label='Quasi-Sqrt(classical cost)')

	# EXPECTED GROVER COST
	ax.axhline(y=configSecurity.aes_expected_security[n][max_depth], dashes=[3, 3], color='#B90E0A', alpha=1, lw=1, label='Expected cost of Grover on AES')

	# TARGET SECURITY
	ax.axhline(y=configSecurity.kyber_expected_security[n], dashes=[1, 1], color='black', alpha=1, lw=1, label='Canonical bit security of Kyber')

	# EXPECTED CLASSICAL
	ax.axhline(y=log_cost_classical, dashes=[2, 6], color='black', alpha=1, lw=1, label='Expected cost of classical enumeration')



def setup_plot(ax, x_coords, y_coords, md, n, log_M, axis_fontsize=14):
	"""
	Support function to set plot.

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
	x_upper = math.ceil(x_coords[-1])+1
	major_xticks = range(0, x_upper, 10)
	minor_xticks = range(0, x_upper, 2)
	plt.xlim(0, x_upper)
	ax.set_xticks(major_xticks)
	ax.set_xticks(minor_xticks, minor=True)
	ax.tick_params(labelsize=axis_fontsize)

	# Y Axis
	# ONLY USED FOR AXIS
	log_N_kh_M = TH.log_avg_N_kh_M(k=0, h=n + 1, n=n, nbases=2 ** log_M, pprime=1, bound='LBUB')
	quasi_sqrt = (log_N_kh_M + math.log(n, 2)) / 2
	y_upper = math.ceil(max(y_coords + [quasi_sqrt, configSecurity.kyber_expected_security[n], configSecurity.aes_expected_security[n][md], log_N_kh_M])) + 10

	if md == 0:
		y_lower = int(min(y_coords + [configSecurity.kyber_expected_security[n], quasi_sqrt, configSecurity.aes_expected_security[n][md]])) - 1
	else:
		y_lower = 0

	major_yticks = range(max(y_lower, 0), y_upper, 50)
	minor_yticks = range(max(y_lower, 0), y_upper, 10)
	plt.ylim(max(y_lower, 0), y_upper)
	ax.set_yticks(major_yticks)
	ax.set_yticks(minor_yticks, minor=True)

	#ax.grid(which='major', alpha=0.5)
	plt.xlabel("Log(Jensen's Gap) z", fontsize=axis_fontsize)
	plt.ylabel("Log", fontsize=axis_fontsize)

def anotate_critical_values(ax, z_to_cost, intersections, max_depth, label_fontsize = 14):
	"""
	Support function annotating intersections.

	:param ax: Axes object from matplotlib
	:param n: BKZ blocksize
	:param max_depth: MaxDepth constraint for depth of quantum circuit
	:param log_M: Bound on number of randomized, combined bases used in cost estimation
	:param z_to_cost: Dictionary [Jenzen's z] --> cost
	:param label_fontsize: Font size of labels.
	:return:
	"""
	def anotate(ax, z_to_cost, max_depth, z_value, inequality = '=', x_offset=0, y_offset=0):
		"""
		Annotation subfunction for all intersections.
		"""
		x_coord = z_value
		y_coord = z_to_cost[z_value].log_gcost_total()

		aLabel = 'z=' + str(z_value) + '\n' + 'Cost' + inequality + str(math.ceil(y_coord))
		if max_depth > 0:
			aLabel += '\n' + 'k=' + str(z_to_cost[z_value].k) + '\n' + 'y=' + str(z_to_cost[z_value].log_y)

		anot = ax.annotate(aLabel, xy=(x_coord, y_coord), xytext=(x_coord + x_offset, y_coord + y_offset),
						   color='black', ha='center', va='center', fontsize=label_fontsize, alpha=1,
						   arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0", linewidth=0.5))

		anot.set_bbox(dict(alpha=.6, linewidth=0, edgecolor='white', facecolor='white'))

	min_z = min(z_to_cost.keys())
	y_offset = max(35, max([z_to_cost[min_z].log_gcost_total()]) / 4)

	label_texts = []

	# First value, ussually z = 0
	aLabel = 'z=' + str(min_z) + '\n' + 'Cost' + '$\leq$' + str(math.ceil(z_to_cost[min_z].log_gcost_total()))
	if max_depth > 0:
		aLabel += '\n' + 'k=' + str(z_to_cost[min_z].k) + '\n' + 'y=' + str(z_to_cost[min_z].log_y)

	label_texts.append(plt.text(0, z_to_cost[min_z].log_gcost_total(), aLabel, fontsize=label_fontsize))

	z_to_ineq = {min_z: '='}

	for k,v in intersections.items():
		if v[0] not in z_to_ineq.keys():
			z_to_ineq[v[0]] = v[1]

	adjust_text_present = False
	try:
		# https://github.com/Phlya/adjustText
		from adjustText import adjust_text
		adjust_text_present = True
	except:
		adjust_text_present = False

	for i,z in enumerate(sorted(z_to_ineq.keys())):
		if z_to_ineq[z] != '>>':
			# z value, s.t. total <= cost

			# Dot
			ax.plot(z, z_to_cost[z].log_gcost_total(), marker="o", markersize=3, color='black', alpha=.8)

			# Ugly text
			if not adjust_text_present:
				y_offset *= -1
				anotate(ax, z_to_cost, max_depth, z, inequality='$\leq$', x_offset=0, y_offset=y_offset)

			aLabel = 'z=' + str(z) + '\n' + 'Cost' + '$\leq$' + str(math.ceil(z_to_cost[z].log_gcost_total()))
			if max_depth > 0:
				aLabel += '\n' + 'k=' + str(z_to_cost[z].k) + '\n' + 'y=' + str(z_to_cost[z].log_y)
			label_texts.append(plt.text(z, z_to_cost[z].log_gcost_total(), aLabel, fontsize=label_fontsize))

	# Prettier text
	if adjust_text_present:
		adjust_text(label_texts, only_move={'points': 'y', 'texts': 'y'}, arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0", linewidth=0.6, mutation_scale=12), avoid_self=True)


def plot_cost_estimation_z_to_g_cost(plot_dir, filename, z_to_cost, intersections, max_depth, n, log_M, constants, bound, leg_font_size = 14, log_loc='right'):
	"""
	Main function called for plotting.

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
	draw_intersecting_lines(ax, n, log_M, max_depth, constants, bound)

	# Draw Costs
	if max_depth > 0:
		ax.plot(x_coords, y_coords_quantum, dashes=[20, 10], color='#0080ff', label=f'Quantum GCost', alpha=0.8)
		ax.plot(x_coords, y_coords, color='#669900', label=f'Total GCost', alpha=0.5)
		ax.plot(x_coords, y_coords_classical, dashes=[24, 2], color='#8D4004', label=f'Classical GCost', alpha=0.8)
		ax.plot(x_coords, qaracm_values, linestyle='dashed', color='#55CEFF', label=f'QRACM', linewidth=0.8, alpha=0.8)
	else:
		#ax.plot(x_coords, y_coords_classical, linestyle='dashed', color='#8D4004', label=f'Expected cost of classical enumeration', alpha=0.8)
		ax.plot(x_coords, y_coords_quantum, dashes=[24, 2], color='#0080ff', label=f'Expected cost of quantum enumeration', alpha=0.8)

	# Find and mark intersections
	anotate_critical_values(ax, z_to_cost, intersections, max_depth)

	plt.savefig(f"{plot_dir}/{filename}_nolegend.pdf", dpi=500, transparent=False, bbox_inches='tight')


	if log_loc == 'right':
		ax.legend(loc='center right', bbox_to_anchor=(2, .5), fontsize=leg_font_size, framealpha=0.8)  # (horizontal, vertical)
	else:
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4), fontsize=leg_font_size, framealpha=0.8)  # (horizontal, vertical)

	plt.savefig(f"{plot_dir}/{filename}.pdf", dpi=500, transparent=False, bbox_inches='tight')

	plt.close()
	plt.cla()
	plt.clf()

