#!/usr/bin/python
from IMP.npctransport import *
import sys
import itertools
import numpy as np

# fetch params
# Usage: <cmd> <outfile>
outfile = sys.argv[1]

config = Configuration()
IMP.npctransport.set_default_configuration(config)
config.statistics_fraction.lower = 1.0
config.interaction_k.lower = 1e-09
config.interaction_range.lower = 10
config.backbone_k.lower = 0.0075
config.time_step_factor.lower = 2.0
config.slack.lower = 10
config.number_of_trials = 1
config.statistics_interval_ns = 1
config.dump_interval_ns = 10
config.output_statistics_interval_ns = 10
config.simulation_time_ns = 1000
config.box_is_on.lower = 1
config.box_side.lower = 1500
config.slab_is_on.lower = 2
config.slab_thickness.lower = 150
config.tunnel_radius.lower = 185
config.is_xyz_hist_stats = 1
config.backbone_tau_ns.lower = 50.0
config.is_backbone_harmonic = 1
config.output_npctransport_version = 4.5
config.xyz_stats_crop_factor = 1
config.xyz_stats_voxel_size_a = 10
config.xyz_stats_max_box_size_a = 700 # -> actually 800
config.nonspecific_k.lower = 0.01
config.excluded_volume_k.lower = 10.0
config.angular_D_factor.lower=0.3
config.is_multiple_hdf5s = True
config.full_output_statistics_interval_factor = 100


anchor_coordinates = [[185.0, 0.0, -75.0], [130.8147545195113, 130.8147545195113, -75.0], [1.1327982892113017e-14, 185.0, -75.0], [-130.81475451951127, 130.8147545195113, -75.0], [-185.0, 2.2655965784226034e-14, -75.0], [-130.81475451951133, -130.81475451951127, -75.0], [-3.398394867633905e-14, -185.0, -75.0], [130.81475451951127, -130.81475451951133, -75.0], [120.0480947161671, 0.0, -37.5], [84.88682184232671, 84.88682184232671, -37.5], [7.350825746894615e-15, 120.0480947161671, -37.5], [-84.8868218423267, 84.88682184232671, -37.5], [-120.0480947161671, 1.470165149378923e-14, -37.5], [-84.88682184232673, -84.8868218423267, -37.5], [-2.2052477240683848e-14, -120.0480947161671, -37.5], [84.88682184232668, -84.88682184232673, -37.5], [120.0480947161671, 0.0, 37.499999999999986], [84.88682184232671, 84.88682184232671, 37.499999999999986], [7.350825746894615e-15, 120.0480947161671, 37.499999999999986], [-84.8868218423267, 84.88682184232671, 37.499999999999986], [-120.0480947161671, 1.470165149378923e-14, 37.499999999999986], [-84.88682184232673, -84.8868218423267, 37.499999999999986], [-2.2052477240683848e-14, -120.0480947161671, 37.499999999999986], [84.88682184232668, -84.88682184232673, 37.499999999999986], [185.0, 0.0, 75.0], [130.8147545195113, 130.8147545195113, 75.0], [1.1327982892113017e-14, 185.0, 75.0], [-130.81475451951127, 130.8147545195113, 75.0], [-185.0, 2.2655965784226034e-14, 75.0], [-130.81475451951133, -130.81475451951127, 75.0], [-3.398394867633905e-14, -185.0, 75.0], [130.81475451951127, -130.81475451951133, 75.0]]
suffix_list = ["_C"] * 12 + ["_N"] * 4
#lowered the amount of beads proportinaly and same with the type: 23,9 -> 12,4

############
# Add FGs: #
############

fgs = []
for i, coordinates in enumerate(anchor_coordinates):
    cur_fg = IMP.npctransport.add_fg_type(config,
                                          type_name=f"fg{i}",
                                          number_of_beads=16,
                                          number=1,
                                          radius=8.0,
                                          interactions=1,
                                          rest_length_factor=1.9,
                                          d_factor=1.0,
                                          interaction_k_factor=1.0,
                                          interaction_range_factor=1.0)
    pos = cur_fg.anchor_coordinates.add()
    pos.x = coordinates[0]
    pos.y = coordinates[1]
    pos.z = coordinates[2]
    fgs.append(cur_fg)
    cur_fg.type_suffix_list.extend(suffix_list)

#################################
# add internal fg cohesiveness: #
#################################
self_k_N = 1.47
self_k_C = 1.32
self_range = 6.00
self_nonspec_k_N = 0.01
self_nonspec_k_C = 0.08
for i in itertools.combinations(range(len(fgs)), 2):
    for suffix0 in ["_N", "_C"]:
        for suffix1 in ["_N", "_C"]:
            interactionFG_FG = IMP.npctransport.add_interaction(config,
                                                        name0=f"fg{i[0]}{suffix0}",
                                                        name1=f"fg{i[1]}{suffix1}",
                                                        interaction_k=0.5*(self_k_N if suffix0=="_N" else self_k_C) + 0.5*(self_k_N if suffix1=="_N" else self_k_C),
                                                        interaction_range=self_range)
            interactionFG_FG.nonspecific_k.lower = np.sqrt((self_nonspec_k_N if suffix0=="_N" else self_nonspec_k_C) * (self_nonspec_k_N if suffix1=="_N" else self_nonspec_k_C))

for i in range(len(fgs)):
    for suffix0 in ["_N", "_C"]:
        for suffix1 in ["_N", "_C"]:
            interactionFG_FG = IMP.npctransport.add_interaction(config,
                                                        name0=f"fg{i}{suffix0}",
                                                        name1=f"fg{i}{suffix1}",
                                                        interaction_k=0.5*(self_k_N if suffix0=="_N" else self_k_C) + 0.5*(self_k_N if suffix1=="_N" else self_k_C),
                                                        interaction_range=self_range)
            interactionFG_FG.nonspecific_k.lower = np.sqrt((self_nonspec_k_N if suffix0=="_N" else self_nonspec_k_C) * (self_nonspec_k_N if suffix1=="_N" else self_nonspec_k_C))


#############
# Add NTRs: #
#############

ntr_vals = [(30,5)]
ntrs = []
for vals in ntr_vals:
    cur_ntr  = IMP.npctransport.add_float_type(config,
                                number=400,
                                radius=vals[0],
                                type_name=f"NTR{vals[0]}",
                                interactions=vals[1],
                                d_factor=1.0,
                                interaction_k_factor=1.0,
                                interaction_range_factor=1.0)
    ntrs.append(cur_ntr)

############################
# Add NTR-FG interactions: #
############################
for i in range(len(fgs)):
    for vals in ntr_vals:
        for suffix in ["_N", "_C"]:
            IMP.npctransport.add_interaction(config,
                                     name0=f"fg{i}{suffix}",
                                     name1=f"NTR{vals[0]}",
                                     interaction_k=2.64,
                                     interaction_range=5.5,
                                     range_sigma0_deg=45.0,
                                     range_sigma1_deg=45.0)

# dump to file
f = open(outfile, "wb")
f.write(config.SerializeToString())
print(config)
