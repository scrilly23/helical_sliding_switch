import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import seaborn as sns
import numpy as np

####USER VARIABLES####
df_path = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20230321_new_backbones/20230321_bgs_sample_z1off_delta_omega1_valine.csv'
socket_df_path = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20230621_socket_test/20230622_0321_bbs_socket_filtered.csv'
outdir = '/Users/stephaniecrilly/Kortemme_lab/helix_sliding/20230621_socket_test'
file_header = '20230718'
#local minima values of z1_offset to show with rectangle on plot
local_min0_min = 8
local_min0_max = 11
local_min1_min = 13
local_min1_max = 16
local_min2_min = 18
local_min2_max = 21

sns.set(font_scale=2)
sns.set_style(style='white')
####
bgs_df = pd.read_csv(df_path)

socket_df = pd.read_csv(socket_df_path)
socket_df = socket_df.rename(columns = {'design_id':'structure'})
socket_df['structure'] = socket_df['structure'] + '.pdb'

df = bgs_df.merge(socket_df, how='left', on='structure')
df['score_for_plot'] = df['score'] * df['socket_call']

delta_omega1_min = df['delta_omega1'].min()
delta_omega1_max = df['delta_omega1'].max()
z1_offset_min = df['z1_offset'].min()
z1_offset_max = df['z1_offset'].max()

#plot original version of e landscape
plt.figure(figsize=(15,6))
ax = sns.scatterplot(data=df, x="z1_offset", y="delta_omega1", marker='s', edgecolor=None, hue="score", zorder=-1)
ax.set_ylim(delta_omega1_min, delta_omega1_max)
ax.set_xlim(z1_offset_min, z1_offset_max)
plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0, title='score')

plt.axvspan(xmin=local_min0_min, xmax=local_min0_max, fill=False, edgecolor='black', linestyle='--', linewidth=3, zorder=1)
plt.axvspan(xmin=local_min1_min, xmax=local_min1_max, fill=False, edgecolor='black', linestyle='--', linewidth=3, zorder=1)
plt.axvspan(xmin=local_min2_min, xmax=local_min2_max, fill=False, edgecolor='black', linestyle='--', linewidth=3, zorder=1)
plt.tight_layout()
#currentAxis = plt.gca()
#someX, someY = 0.5, 1
#currentAxis.add_patch(Rectangle((someX - .5, someY - .5), 4, 6, facecolor="none", ec='k', lw=2))

#plt.show()
plt.savefig(f'{outdir}/{file_header}_0321_bbs_og_e_landscape.png', bbox_inches='tight')
plt.clf()

#plot by socket filter only
#plt.figure(figsize=(15,6))
socket_filtered_df = df.query('socket_call == 0', engine='python').copy()
#ax = sns.scatterplot(data=socket_filtered_df, x="z1_offset", y="delta_omega1", marker='s', edgecolor=None, color="black")
#plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0, title='score')
#ax.set_ylim(delta_omega1_min, delta_omega1_max)
#ax.set_xlim(z1_offset_min, z1_offset_max)

plt.figure(figsize=(15,6))
ax = sns.scatterplot(data=df, x="z1_offset", y="delta_omega1", marker='s', edgecolor=None, hue="score", zorder=-1)
ax.set_ylim(delta_omega1_min, delta_omega1_max)
ax.set_xlim(z1_offset_min, z1_offset_max)
sns.scatterplot(data=socket_filtered_df, x="z1_offset", y="delta_omega1", marker='s', edgecolor=None, color="black", ax=ax)
plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0, title='score')

plt.axvspan(xmin=local_min0_min, xmax=local_min0_max, fill=False, edgecolor='yellow', linestyle='--', linewidth=3, zorder=1)
plt.axvspan(xmin=local_min1_min, xmax=local_min1_max, fill=False, edgecolor='yellow', linestyle='--', linewidth=3, zorder=1)
plt.axvspan(xmin=local_min2_min, xmax=local_min2_max, fill=False, edgecolor='yellow', linestyle='--', linewidth=3, zorder=1)
plt.tight_layout()
#plt.show()
plt.savefig(f'{outdir}/{file_header}_0321_bbs_socket_filtered_e_landscape.png', bbox_inches='tight')

#save df of structures with detected cc's only
ccs_detected_df = df.query('socket_call == 1', engine='python').copy()
ccs_detected_df = ccs_detected_df.drop(['Unnamed: 0_x', 'Unnamed: 0_y'], axis=1)

#ccs_detected_df.to_csv(f'{outdir}/{file_header}_0321_bbs_socket_cc_detected.csv')
