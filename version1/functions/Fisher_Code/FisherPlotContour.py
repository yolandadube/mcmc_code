import os, sys
import numpy as np
#---------------------------------------------------------------------------------------------------------------------

def plot_contours(mainPath, dictionary):

	sys.path.append(mainPath+"functions/planckBestFit/")
	#---#
	import seaborn as sns
	from getdist import plots, MCSamples
	from planckCosmologicalParameters import planckValues
	#---#

	# fiducial Planck cosmological parameters
	PlanckData = planckValues(dictionary)
	#---#
	A_s, sigma80, H0, h, omm0, omb0, \
	omcdm0, omega_k0, omega_n0, n_s, gamma, w, fnl = PlanckData

	# load fisher matrices
	Path  = mainPath+'result/contour/FisherforecastingDATA/'

	data1 = np.loadtxt(Path+"cov_matrix_skaOIMBand1.dat")
	data2 = np.loadtxt(Path+"cov_matrix_meerKATUHFBand.dat")
	data3 = np.loadtxt(Path+"cov_matrix_meerKATLBand.dat")
	#--------------------------------------------------------------------------------------------------------------------------


	# define your fisher parameters
	alpha = 0.45
	beta  = 0.30
	P0    = 37300
	k0    = 0.0163

	# fiducial values of parameters considered
	x0    = alpha 
	x1    = beta
	x2    = P0
	x3    = k0
	
	# fiducial values of Fisher parameters
	meanValue = np.array([x0, x1, x2, x3])       

	#--------------------------------------------------------------------------------------------------------------------------

	# covariance matrices
	cov1 = data1
	cov2 = data2
	cov3 = data3
	#--------------------------------------------------------------------------------------------------------------------------

	# axis labels
	axisLabels    = [r'$\alpha$', r'$\beta$', r'$P_{0}$', \
					 r'$k_{0}$']
	axisNames     = ["x%s"%i for i in range(len(axisLabels))]
	
	# creating a normal distribution for the parameters
	numberOfSamplingPoints = 300000
	normalDist1 = np.random.multivariate_normal(meanValue, cov1, size=numberOfSamplingPoints)
	normalDist2 = np.random.multivariate_normal(meanValue, cov2, size=numberOfSamplingPoints)
	normalDist3 = np.random.multivariate_normal(meanValue, cov3, size=numberOfSamplingPoints)

	# plotting
	samplesMC1 = MCSamples(samples=normalDist1, names=axisNames, labels=axisLabels)
	samplesMC2 = MCSamples(samples=normalDist2, names=axisNames, labels=axisLabels)
	samplesMC3 = MCSamples(samples=normalDist3, names=axisNames, labels=axisLabels)
	#--------------------------------------------------------------------------------------------------------------------------

	# plotting
	#g = plots.get_subplot_plotter(rc_sizes=True, width_inch=None)
	g = plots.get_single_plotter(rc_sizes=True, width_inch=None)

	# customizing subplots
	g.settings.scaling                    = True
	g.settings.linewidth                  = 3.0
	#g.settings.line_styles                = ['-', '-']
	#g.settings.line_labels                = False # remove legend
	g.settings.prob_y_ticks               = True
	g.settings.tight_layout               = True
	g.settings.axes_fontsize              = 20
	g.settings.axes_labelsize             = 20
	g.settings.legend_fontsize            = 20
	g.settings.norm_prob_label            = 20
	g.settings.norm_1d_density            = True
	g.settings.num_plot_contours          = 2
	g.settings.figure_legend_loc          = "upper right"
	g.settings.subplot_size_ratio         = 1
	g.settings.figure_legend_frame        = True
	g.settings.axis_tick_max_labels       = 3
	g.settings.alpha_factor_contour_lines = 1

	# colours, line styles, labels
	cases      = [samplesMC3, samplesMC2, samplesMC1]
	colors     = [sns.xkcd_rgb["blue"], sns.xkcd_rgb["green"], sns.xkcd_rgb["red"]]
	linestyles = ['-', '-', '-']
	labelS     = [r'MeerKAT L-BAND', r'MeerKAT UHF-BAND', r'SKA-MID BAND 1']

	g.triangle_plot(
    cases,
    ['x2', 'x3', 'x0', 'x1'],  # Reordered to match the desired axis order
    filled_compare=True,
    filled=False,
    contour_colors=colors,
    line_args=[
        {'ls': linestyles[i], 'color': colors[i]}
        for i in range(len(cases))
    ],
    legend_labels=labelS,
    alpha=1.0
)

	#------------^^^------------^^^------------^^^------------^^^------------^^^------------^^^------------^^^------------

	# Customizing axes limits for each subplot

	# # Subplot (0,0): alpha vs. alpha
	# ax00 = g.subplots[0, 0]
	# ax00.set_xlim([-3, 3])  # Custom limit for alpha
	# #ax00.set_ylim([, 1])  # Custom limit for alpha
	# #ax00.vlines(x0, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# # Subplot (1,0): beta vs. alpha
	# ax10 = g.subplots[1, 0]
	# ax10.set_xlim([-2, 2])  # Custom limit for alpha
	# ax10.set_ylim([-1, 2])    # Custom limit for beta
	# #ax10.vlines(x0, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)
	# #ax10.hlines(x1, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# # Subplot (2,0): P0 vs. alpha
	# ax20 = g.subplots[2, 0]
	# ax20.set_xlim([-3, 3])  # Custom limit for alpha
	# ax20.set_ylim([34500, 41000])  # Custom limit for P0
	# #ax20.vlines(x0, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)
	# #ax20.hlines(x2, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# # Subplot (3,0): k0 vs. alpha
	# ax30 = g.subplots[3, 0]
	# ax30.set_xlim([-3, 4])  # Custom limit for alpha
	# ax30.set_ylim([-1, 1])  # Custom limit for k0
	# # ax30.vlines(x0, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)
	# # ax30.hlines(x3, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# # Subplot (1,1): beta vs. beta
	# ax11 = g.subplots[1, 1]
	# ax11.set_xlim([-2, 2])    # Custom limit for beta
	# #ax11.set_ylim([0, 0.5])    # Custom limit for beta
	# ax11.vlines(x1, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# Subplot (2,1): P0 vs. beta
	#ax21 = g.subplots[2, 1]
	#ax21.set_xlim([-2, 2])    # Custom limit for beta
	#ax21.set_ylim([35000, 40000])  # Custom limit for P0
	#ax21.vlines(x1, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)
	#ax21.hlines(x2, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	#Subplot (3,1): k0 vs. beta
	#ax31 = g.subplots[3, 1]
	#ax31.set_xlim([-3, 3])    # Custom limit for beta
	#ax31.set_ylim([-0.6, 0.6])  # Custom limit for k0
	# ax31.vlines(x1, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)
	# ax31.hlines(x3, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# # Subplot (2,2): P0 vs. P0
	# ax22 = g.subplots[2, 2]
	# ax22.set_xlim([35000, 40000])  # Custom limit for P0
	# #ax22.set_ylim([20000, 38000])  # Custom limit for P0
	# ax22.vlines(x2, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	# # Subplot (3,2): k0 vs. P0
	#ax32 = g.subplots[3, 2]
	#ax32.set_xlim([32000, 42000])  # Custom limit for P0
	# ax32.set_ylim([0.012, 0.018])  # Custom limit for k0
	# ax32.vlines(x2, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)
	# ax32.hlines(x3, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	#Subplot (3,3): k0 vs. k0
	#ax33 = g.subplots[3, 3]
	#ax33.set_xlim([-1.0, 1.0])  # Custom limit for k0
	#ax33.set_ylim([0.012, 0.018])  # Custom limit for k0
	#ax33.set_title(r"$k_0=0.0163$", fontsize=16)
	# ax33.vlines(x3, -1000.0, 1000.0, color=sns.xkcd_rgb["black"], linestyle=':', linewidth=1.5, alpha=1.0)

	
	# saving as pdf
	g.export(Path+'FisherContourNew.pdf')

	return None
#---------------------------------------------------------------------------------------------------------------------