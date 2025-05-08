import matplotlib.pyplot as plt
from propagate import *
# from averaging_effects import *
# import cv2
# from scipy.spatial import Delaunay
import plotly.graph_objects as go
# import  moviepy.editor as mpy
# import io 
# from PIL import Image
# import os
# import matplotlib
# matplotlib.use('Agg')
from plotly_gif import GIF, two_d_time_series
from matplotlib.ticker import MaxNLocator

# # Plot the intensity as a function of time
# plt.figure(figsize=(8, 6))
# plt.plot(system_params.gamma_q1 * times, intensity_results)
# # plt.plot(system_params.gamma_q1 * times, np.load('files/H_strongC.npy') - np.load('files/J_strongC.npy'))
# # plt.plot(system_params.gamma_q1 * times * times, 2 * np.exp(-system_params.gamma_q1 * times * times), label = 'theoretical')
# # plt.plot(gamma_q1*times, 2*np.exp(-2*gamma_q1*times)*(1 + 2*gamma_q1*times), label='theoretical')
# plt.xlabel(r'$\gamma t$', fontsize = 14)
# plt.ylabel(r'$I(t) / I_0$', fontsize = 14)
# # plt.semilogx()
# # plt.xlim([0, 5])  
# # plt.ylim([0, 2.1])
# plt.title('Intensity vs time', fontsize = 18)
# plt.legend(loc='best')
# plt.grid(True)
# plt.show()

# # Plot the intensity as a function of time
# plt.figure(figsize=(8, 6))
# # plt.plot(times, pop_results, label = 'H_noC_dark')
# plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_bright_T0.npy'), label = 'Bright (polaron)')
# plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_bright_T0_fft.npy'), ':', label = 'Bright (polaron) - FFT')
# plt.plot(system_params.gamma_q1 * times, np.load('H_weakC_bright_T0.npy'), '--', label = 'Bright (WCME)')
# plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_dark_T0.npy'), label = 'Dark (polaron)')
# plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_dark_T0_fft.npy'), ':', label = 'Dark (polaron) - FFT')
# plt.plot(system_params.gamma_q1 * times, np.load('H_weakC_dark_T0.npy'), '--', label = 'Dark (WCME)')
# # plt.plot(gamma_q1 * times, 2 * np.exp(-gamma_q1 * times), label = 'theoretical')
# # plt.plot(gamma_q1*times, 2*np.exp(-2*gamma_q1*times)*(1 + 2*gamma_q1*times), label='theoretical')
# plt.xlabel(r'$\gamma t$', fontsize = 14)
# plt.ylabel(r'Population', fontsize = 14)
# # plt.semilogx()
# # plt.xlim([0, 5])  
# # plt.ylim([-0.1, 1])
# plt.title('H-dimer', fontsize = 18)
# plt.legend(loc='best')
# plt.grid(True)
# plt.savefig('H_dimer_weakC.jpg', dpi = 200)
# plt.show()


# # Plot the intensity as a function of time
# plt.figure(figsize=(8, 6))
# # plt.plot(times, pop_results, label = 'H_noC_dark')
# # plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_bright_Tv0.npy'), label = 'Bright (polaron)')
# plt.plot(system_params.gamma_q1 * times, np.load('J_strongC_bright_Tv0_fft.npy'), label = 'Bright (polaron) - FFT')
# plt.plot(system_params.gamma_q1 * times, np.load('J_weakC_bright_Tv0.npy'), '--', label = 'Bright (WCME)')
# # plt.plot(system_params.gamma_q1 * times, np.load('H_strongC_dark_T0.npy'), label = 'Dark (polaron)')
# plt.plot(system_params.gamma_q1 * times, np.load('J_strongC_dark_Tv0_fft.npy'), label = 'Dark (polaron) - FFT')
# plt.plot(system_params.gamma_q1 * times, np.load('J_weakC_dark_Tv0.npy'), '--', label = 'Dark (WCME)')
# # plt.plot(gamma_q1 * times, 2 * np.exp(-gamma_q1 * times), label = 'theoretical')
# # plt.plot(gamma_q1*times, 2*np.exp(-2*gamma_q1*times)*(1 + 2*gamma_q1*times), label='theoretical')
# plt.xlabel(r'$\gamma t$', fontsize = 14)
# plt.ylabel(r'Population', fontsize = 14)
# # plt.semilogx()
# # plt.xlim([0, 5])  
# # plt.ylim([-0.1, 1])
# plt.title('J-dimer', fontsize = 18)
# plt.legend(loc='best')
# plt.grid(True)
# plt.savefig('J_dimer_weakC_Tv0.jpg', dpi = 200)
# plt.show()

'''Population and intensity'''

# # Load data
# # J_weakC_Tinf = np.load('files/intensity/J_dimer/J_weakC_Tinf.npy')
# # J_noC_Tinf = np.load('files/intensity/J_dimer/J_noC_Tinf.npy')
# H_weakC_Tinf = np.load('files/intensity/H_dimer/H_weakC_Tinf.npy')
# H_noC_Tinf = np.load('files/intensity/H_dimer/H_noC_Tinf.npy')

# # Calculate the difference in intensitie
# # J_intensity_diff_Tinf = J_weakC_Tinf - J_noC_Tinf
# H_intensity_diff_Tinf = H_weakC_Tinf - H_noC_Tinf

# # Plot for H dimer
# fig, ax1 = plt.subplots(figsize=(8, 12))  # Keep figure size consistent with previous code


# ax1.set_ylabel('Population', fontsize=35)
# ax1.plot(system_params.gamma_q1 * times, np.load('H_weakC_bright_Tinf.npy'), label='Bright state',
#     linewidth=2.5, linestyle='-', color='blue')
# ax1.plot(system_params.gamma_q1 * times, np.load('H_weakC_dark_Tinf.npy'), label='Dark state',
#     linewidth=2.5, linestyle='--', color='tomato')
# ax1.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax1.tick_params(axis='y', labelsize=30, direction='in', length=6)
# # ax1.set_xscale('log')
# # ax1.set_yscale('log')
# # ax1.legend(fontsize=30, loc='upper right', frameon=False)

# ax2 = ax1.twinx()  # Instantiate a second axes that shares the same x-axis

# ax1.set_xlabel(r'$\gamma t$', fontsize=35)
# ax2.set_ylabel(r'Difference $\Delta I/I_{0}$', fontsize=35)
# ax2.plot(system_params.gamma_q1 * times, H_intensity_diff_Tinf, label=r'Difference $\Delta I/I_{0}$',
#          linewidth=2.5, linestyle=':', color='black')
# ax2.tick_params(axis='y', labelsize=30, direction='in', length=6)
# # ax2.set_xscale('log')
# # ax2.set_yscale('log')
# # ax1.legend(fontsize=30, loc='best', frameon=False)

# # Add (a) label
# ax1.text(0.2, 0.9, '(c)', transform=ax1.transAxes, fontsize=35, verticalalignment='top')


# # Custom legend positioning on the right side of the plot
# lines_labels = [ax1.get_lines()[0], ax1.get_lines()[1], ax2.get_lines()[0]]
# labels = [line.get_label() for line in lines_labels]
# # fig.legend(lines_labels, labels, loc='center left', bbox_to_anchor=(0.45, 0.5), fontsize=30, frameon=False)


# # Adjust layout to ensure text isn't cut off
# fig.tight_layout(pad=3.0)
# # plt.title(r'H Dimer: Population and Intensity Difference', fontsize=35)
# plt.savefig(r'H_pop_int_weak_Tinf.png', dpi=500, bbox_inches='tight')
# plt.show()


# # Plot for J dimer
# fig, ax1 = plt.subplots(figsize=(12, 10))  # Keep figure size consistent with previous code


# ax1.set_ylabel('Population', fontsize=35)
# ax1.plot(system_params.gamma_q1 * times, np.load('J_weakC_bright_T0.npy'), label='Bright state',
#     linewidth=2.5, linestyle='-', color='blue')
# ax1.plot(system_params.gamma_q1 * times, np.load('J_weakC_dark_T0.npy'), label='Dark state',
#     linewidth=2.5, linestyle='--', color='tomato')
# ax1.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax1.tick_params(axis='y', labelsize=30, direction='in', length=6)
# # ax1.set_xscale('log')
# # ax1.set_yscale('log')
# # ax1.legend(fontsize=30, loc='upper right', frameon=False)

# ax2 = ax1.twinx()  # Instantiate a second axes that shares the same x-axis

# ax2.set_xlabel(r'$\gamma t$', fontsize=35)
# ax2.set_ylabel(r'Difference $\Delta I/I_{0}$', fontsize=35)
# ax2.plot(system_params.gamma_q1 * times, J_intensity_diff_T0, label=r'Difference $\Delta I/I_{0}$',
#          linewidth=2.5, linestyle=':', color='black')
# ax2.tick_params(axis='y', labelsize=30, direction='in', length=6)
# # ax1.legend(fontsize=30, loc='best', frameon=False)
# # ax2.set_xscale('log')
# # ax2.set_yscale('log')

# # Custom legend positioning on the right side of the plot
# lines_labels = [ax1.get_lines()[0], ax1.get_lines()[1], ax2.get_lines()[0]]
# labels = [line.get_label() for line in lines_labels]
# fig.legend(lines_labels, labels, loc='center left', bbox_to_anchor=(0.45, 0.5), fontsize=30, frameon=False)


# # Adjust layout to ensure text isn't cut off
# fig.tight_layout(pad=3.0)
# # plt.title(r'J Dimer: Population and Intensity Difference', fontsize=35)
# plt.savefig(r'J_pop_int_weak_T0.png', dpi=500, bbox_inches='tight')
# plt.show()






# # Load data
# J_weakC_T0 = np.load('files/intensity/J_dimer/J_weakC_T0.npy')
# J_noC_T0 = np.load('files/intensity/J_dimer/J_noC_T0.npy')
# H_weakC_T0 = np.load('files/intensity/H_dimer/H_weakC_T0.npy')
# H_noC_T0 = np.load('files/intensity/H_dimer/H_noC_T0.npy')

# # Calculate the difference in intensities
# J_intensity_diff_T0 = J_weakC_T0 - J_noC_T0
# H_intensity_diff_T0 = H_weakC_T0 - H_noC_T0

# # Create a figure with two subplots stacked vertically
# fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(14, 12), sharex=True)  # Share x-axis

# # Plot for H dimer
# ax1.set_ylabel('Population', fontsize=35)
# ax1.plot(system_params.gamma_q1 * times, np.load('H_weakC_bright_T0.npy'), label='Bright state',
#          linewidth=2.5, linestyle='-', color='blue')
# ax1.plot(system_params.gamma_q1 * times, np.load('H_weakC_dark_T0.npy'), label='Dark state',
#          linewidth=2.5, linestyle='--', color='tomato')
# ax1.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax1.tick_params(axis='y', labelsize=30, direction='in', length=6)

# ax2 = ax1.twinx()  # Secondary y-axis for intensity difference
# ax2.set_ylabel(r'Difference $\Delta I/I_{0}$', fontsize=35)
# ax2.plot(system_params.gamma_q1 * times, H_intensity_diff_T0, label=r'Difference $\Delta I/I_{0}$',
#          linewidth=2.5, linestyle=':', color='black')
# ax2.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # Add (a) label
# ax1.text(0.02, 0.9, '(a)', transform=ax1.transAxes, fontsize=35, verticalalignment='top')

# # Plot for J dimer
# ax3.set_ylabel('Population', fontsize=35)
# ax3.plot(system_params.gamma_q1 * times, np.load('J_weakC_bright_T0.npy'), label='Bright state',
#          linewidth=2.5, linestyle='-', color='blue')
# ax3.plot(system_params.gamma_q1 * times, np.load('J_weakC_dark_T0.npy'), label='Dark state',
#          linewidth=2.5, linestyle='--', color='tomato')
# ax3.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax3.tick_params(axis='y', labelsize=30, direction='in', length=6)
# ax3.set_xlabel(r'$\gamma t$', fontsize=35)

# ax4 = ax3.twinx()  # Secondary y-axis for intensity difference
# ax4.set_ylabel(r'Difference $\Delta I/I_{0}$', fontsize=35)
# ax4.plot(system_params.gamma_q1 * times, J_intensity_diff_T0, label=r'Difference $\Delta I/I_{0}$',
#          linewidth=2.5, linestyle=':', color='black')
# ax4.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # Add (b) label
# ax3.text(0.02, 0.9, '(b)', transform=ax3.transAxes, fontsize=35, verticalalignment='top')

# # Combine legends from both plots and display only once
# lines_labels = [ax1.get_lines()[0], ax1.get_lines()[1], ax2.get_lines()[0]]
# labels = [line.get_label() for line in lines_labels]
# fig.legend(lines_labels, labels, loc='center left', bbox_to_anchor=(0.5, 0.7), fontsize=30, frameon=False)

# # Adjust layout to ensure text isn't cut off
# fig.tight_layout(pad=3.0)

# # Save and show the combined figure
# plt.savefig(r'combined_H_J_pop_int_weak_T0.png', dpi=500, bbox_inches='tight')
# plt.show()


'''Distribution of intensty'''


# # Generate R_vals and Theta_vals for the polar plot
# R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# # Create the polar heatmap
# fig = plt.figure(figsize=(12, 8))  # Adjusted figure size
# ax1 = fig.add_subplot(projection='polar')

# # Create the colormap on the polar plot
# c = ax1.pcolormesh(Theta_vals, R_vals, norm_int, cmap='coolwarm', shading='auto', vmin=0, vmax=1)

# # Set labels and formatting
# ax1.set_xlabel(r'$\theta$ (degrees)', fontsize=55, labelpad=0) 
# ax1.set_ylabel(r'$I(\theta)$', fontsize=55, labelpad=80)
# ax1.tick_params(axis='x', labelsize=50, pad=30, direction='in', length=6)  # pad=30 moves angle labels further outside
# ax1.tick_params(axis='y', labelsize=50, direction='in', length=6)

# # # Move the text label to the upper-left corner
# # ax1.text(-0.3, 1.05, r'(c)', size=50, transform=ax1.transAxes)  # Adjusted position to upper-left

# # Remove y-tick labels (optional based on your preference)
# ax1.set_yticklabels([])

# # Adjust layout to prevent cutting off theta labels
# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)


# # Add colorbar with specific size and position adjustments
# cbar = plt.colorbar(c, ax=ax1, shrink=1.0, aspect=10, pad=0.1)  # Shrink and aspect to control size, pad to move it away

# # Adjust text size on colorbar
# cbar.ax.tick_params(labelsize=50)  # Change the number to adjust the colorbar tick label size

# # # Save the figure with the desired format and settings
# # plt.savefig('dist_int_ortho.png', dpi=500, bbox_inches='tight')

# # Display the plot
# plt.show()




# # Read the CSV file
# df = pd.read_csv('files/dipole_positions/dipole_pos_dir_H.csv')

# # Remove leading spaces from column names
# df.columns = df.columns.str.strip()


# # # Normalize the intensity values
# norm_int = norm_int / norm_int.max()

# # Create figure
# fig = go.Figure()

# # Add initial surface plot with dynamic color scaling
# # cmin = norm_int.min()
# # cmax = norm_int.max()
# fig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=norm_int, colorscale='hot', opacity=0.7))

# # Add dipoles (constant for all frames in this example)
# # Plot dipoles
# for index, row in df.iterrows():
#     # Calculate end position of the arrow shaft (half length in each direction)
#     half_length = 0.1  # Half length of the arrow shaft
#     start_x = row['x'] - 0.20 * row['dir_x']
#     start_y = row['y'] - 0.20 * row['dir_y']
#     start_z = row['z'] - 0.20 * row['dir_z']
#     end_x = row['x'] + 0.20 * row['dir_x']
#     end_y = row['y'] + 0.20 * row['dir_y']
#     end_z = row['z'] + 0.20 * row['dir_z']

#     # Plot dipole position
#     fig.add_trace(go.Scatter3d(
#         x=[row['x']],
#         y=[row['y']],
#         z=[row['z']],
#         mode='markers',
#         marker=dict(color='red', size=10),
#         name=''
#     ))

#     # Plot dipole direction shaft
#     fig.add_trace(go.Scatter3d(
#         x=[start_x, end_x],
#         y=[start_y, end_y],
#         z=[start_z, end_z],
#         mode='lines',
#         line=dict(color='blue', width=10),
#         name=''
#     ))

#     # Plot dipole direction arrowhead
#     fig.add_trace(go.Cone(
#         x=[end_x],
#         y=[end_y],
#         z=[end_z],
#         u=[row['dir_x']],
#         v=[row['dir_y']],
#         w=[row['dir_z']],
#         showscale=False,
#         sizemode="absolute",
#         sizeref=0.25,
#         anchor="tail",
#         opacity=1.0,
#         name=''
#     ))


# # # Plot point on the sphere at (1, 0, 0)
# # fig.add_trace(go.Scatter3d(
# #     x=[perp_ortho[0]],
# #     y=[perp_ortho[1]],
# #     z=[perp_ortho[2]],
# #     mode='markers',
# #     marker=dict(color='green', size=10),
# #     name='Sphere Point'
# # ))


# # # Display the animation
# # fig.show()



# # # Create frames with dynamic color scaling
# # frames = [go.Frame(
# #     data=[go.Surface(z=z, surfacecolor=norm_int[:, :, k], colorscale='hot', opacity=0.7)], name=f'frame{k}') for k in range(200)]



# # # Normalize the g2 values
# # g2_corr = np.interp(np.array(g2_corr), (np.array(g2_corr).min(), np.array(g2_corr).max()), (0, 1))  #/ g2_corr.max()

# # # Create figure
# # fig = go.Figure()
# # cmin = g2_corr.min() #0.0 #0.0 * 0.007
# # cmax = g2_corr.max() #1.0 #1.0 * 0.0144


# # # Add initial surface plot with dynamic color scaling
# # fig.add_trace(go.Surface(
# #     x=x, y=y, z=z, surfacecolor=g2_corr, colorscale='hot', opacity=0.6, cmin=cmin, cmax=cmax
# # ))

# # # Add dipoles (constant for all frames in this example)
# # # Plot dipoles
# # for index, row in df.iterrows():
# #     # Calculate end position of the arrow shaft (half length in each direction)
# #     half_length = 0.1  # Half length of the arrow shaft
# #     start_x = row['x'] - 0.20 * row['dir_x']
# #     start_y = row['y'] - 0.20 * row['dir_y']
# #     start_z = row['z'] - 0.20 * row['dir_z']
# #     end_x = row['x'] + 0.20 * row['dir_x']
# #     end_y = row['y'] + 0.20 * row['dir_y']
# #     end_z = row['z'] + 0.20 * row['dir_z']

# #     # Plot dipole position
# #     fig.add_trace(go.Scatter3d(
# #         x=[row['x']],
# #         y=[row['y']],
# #         z=[row['z']],
# #         mode='markers',
# #         marker=dict(color='red', size=10),
# #         name=''
# #     ))

# #     # Plot dipole direction shaft
# #     fig.add_trace(go.Scatter3d(
# #         x=[start_x, end_x],
# #         y=[start_y, end_y],
# #         z=[start_z, end_z],
# #         mode='lines',
# #         line=dict(color='blue', width=10),
# #         name=''
# #     ))

# #     # Plot dipole direction arrowhead
# #     fig.add_trace(go.Cone(
# #         x=[end_x],
# #         y=[end_y],
# #         z=[end_z],
# #         u=[row['dir_x']],
# #         v=[row['dir_y']],
# #         w=[row['dir_z']],
# #         showscale=False,
# #         sizemode="absolute",
# #         sizeref=0.25,
# #         anchor="tail",
# #         opacity=1.0,
# #         name=''
# #     ))

# # # Plot point on the sphere at (1, 0, 0)
# # fig.add_trace(go.Scatter3d(
# #     x=[perp_ortho[0]],
# #     y=[perp_ortho[1]],
# #     z=[perp_ortho[2]],
# #     mode='markers',
# #     marker=dict(color='green', size=10),
# #     name=''
# # ))




# # # Create frames with dynamic color scaling
# # frames = [go.Frame(
# #     data=[go.Surface(z=z, surfacecolor=g2_corr[:, :, k], colorscale='hot', opacity=0.7, cmin=cmin, cmax=cmax)], name=f'frame{k}') for k in range(200)]


# # # Add frames to the figure
# # fig.frames = frames

# # # Update layout with animation settings
# # fig.update_layout(
# #     scene=dict(
# #         xaxis_title='X',
# #         yaxis_title='Y',
# #         zaxis_title='Z',
# #         aspectmode='cube'
# #     ),
# #     title='Normalised Intensity Plot on a Sphere',
# #     updatemenus=[{
# #         'buttons': [
# #             {
# #                 'args': [None, {'frame': {'duration': 100, 'redraw': True}, 'fromcurrent': True}],
# #                 'label': 'Play',
# #                 'method': 'animate'
# #             },
# #             {
# #                 'args': [[None], {'frame': {'duration': 0, 'redraw': True}, 'mode': 'immediate', 'transition': {'duration': 0}}],
# #                 'label': 'Pause',
# #                 'method': 'animate'
# #             }
# #         ],
# #         'direction': 'left',
# #         'pad': {'r': 10, 't': 87},
# #         'showactive': False,
# #         'type': 'buttons',
# #         'x': 0.1,
# #         'xanchor': 'right',
# #         'y': 0,
# #         'yanchor': 'top'
# #     }]
# # )
# # # fig.write_html("g2_my.html")
# # # Display the animation
# fig.show()






# # Create 100 points between [1.0, 1.0, 0.0] and [-1.0, 1.0, 0.0]
# perp_start = np.array([1.0, 1.0, 0.0])
# perp_end = np.array([-1.0, 1.0, 0.0])
# perp_points = np.linspace(perp_start, perp_end, 100)


# # Create directory to store plots
# output_dir = "frames"
# os.makedirs(output_dir, exist_ok=True)

# # Function to compute and save a single frame
# def compute_and_plot_frame(args):
#     i, perp = args
#     r, theta, phi = cart2sph(perp)
#     r_p, theta_p, phi_p = r, theta, phi

#     # Compute g2 correlation
#     g2_corr = quantum_analysis_g2_dist.compute(theta, phi, theta_p, phi_p)

#     # Plot g2 correlation
#     fig, ax = plt.subplots(figsize=(8, 6))
#     ax.plot(
#         np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns),
#         np.append(np.flip(g2_corr[1:]), g2_corr),
#         label=f"Frame {i+1}"
#     )
#     ax.grid(True)
#     ax.set_title(f"g2 Correlation for Frame {i+1}")
#     ax.set_xlabel("Time (ns)")
#     ax.set_ylabel("g2 Correlation")
#     ax.legend()

#     # Add polar plot inset
#     phi_deg = np.rad2deg(phi)  # Convert phi to degrees
#     theta_deg = np.rad2deg(theta)  # Convert theta to degrees

#     inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25], polar=True)  # Create polar inset
#     inset_ax.set_thetamin(45)  # Start angle: π/4 (in degrees)
#     inset_ax.set_thetamax(135)  # End angle: 3π/4 (in degrees)

#     # Plot a dashed line and the red dot for phi and theta
#     inset_ax.plot([np.deg2rad(45), phi], [0, theta], linestyle="--", color="gray")  # Convert 45° to radians
#     inset_ax.plot(phi, theta, "ro", label=r"$\theta, \phi$")  # Use LaTeX for the label

#     # Set radial limits for theta
#     inset_ax.set_ylim(0, np.pi)  # Theta from 0 to π
#     inset_ax.set_yticks([0, np.pi / 2, np.pi])  # Tick marks at 0, π/2, π
#     inset_ax.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"], fontsize=10)  # Use LaTeX for labels

#     # Convert theta ticks to degrees for display
#     inset_ax.set_xticks(np.radians([45, 60, 75, 90, 105, 120, 135]))  # Ticks every 15° within the range
#     inset_ax.set_xticklabels([r"$45^\circ$", r"$60^\circ$", r"$75^\circ$", r"$90^\circ$", 
#                             r"$105^\circ$", r"$120^\circ$", r"$135^\circ$"], fontsize=10)

#     # Add a legend
#     inset_ax.legend(loc="lower left", fontsize="small")

#     # Save the frame
#     plt.savefig(f"{output_dir}/frame_{i:03d}.png")
#     plt.close()

# # Prepare arguments for parallel execution
# args = [(i, perp) for i, perp in enumerate(perp_points)]

# # Use multiprocessing to generate frames in parallel
# if __name__ == "__main__":
#     num_processes = min(cpu_count(), 10)  # Adjust the number of processes if needed
#     with Pool(num_processes) as pool:
#         pool.map(compute_and_plot_frame, args)

#     # Use ffmpeg to create a video
#     os.system("ffmpeg -r 10 -i frames/frame_%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p video.mp4")







# # Plot the steady state second order correlation as a function of time
# plt.figure(figsize=(8, 6))
# plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[1:]), g2_corr), label = 'Measurement-induced')
# plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(np.load('g2_mic_deph.npy')[1:]), np.load('g2_mic_deph.npy')), label = 'Measurement-induced - pure dephasing')
# plt.xlabel(r'$\tau$ (ns)', fontsize = 14)
# plt.ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize = 14)
# # plt.xlim([-5, 5])
# plt.ylim([0.0, 1.1])
# plt.title(r'Photon coincidence with $\gamma_{\mathrm{pump}} = \gamma$', fontsize = 18)
# plt.legend(loc='best')
# plt.grid(True)
# plt.savefig('g2-ortho-parallel.png')
# # plt.show()













# # Generate R_vals and Theta_vals for the polar plot
# R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# ax.plot(Theta_vals, norm_int, color='black', linewidth = 1)
# ax.grid(True)
# ax.set_title("Polar Heatmap of Intensity on a Sphere", va='bottom')
# plt.savefig('2D polar J.png', dpi = 1000)
# # plt.show()


# # Generate R_vals and Theta_vals for the polar plot
# R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# # Create the polar heatmap
# fig = plt.figure(layout='constrained')
# ax1 = fig.add_subplot(projection='polar')
# c = ax1.pcolormesh(Theta_vals, R_vals, g2_corr[:, :, 0], cmap='hot', shading='auto')

# # Add colorbar
# plt.colorbar(c, label='Intensity')
# plt.title('Polar Heatmap of g(2) on a Sphere')


# # Add arrows
# # Arrow facing theta = 0 (along the x-axis)
# ax1.annotate('', xy=(0, 0.5), xytext=(0, 0),
#              arrowprops=dict(facecolor='blue', shrink=0, width=2, headwidth=10))

# # Arrow facing theta = 90 degrees (along the y-axis)
# ax1.annotate('', xy=(np.pi/2, 0.5), xytext=(0, 0),
#              arrowprops=dict(facecolor='blue', shrink=0, width=2, headwidth=10))

# # Save the figure
# plt.savefig('2D polar - 2.png')











# # Function to plot intensity distribution
# def plot_intensity_distribution(intensity_data, label, filename):
   
#     angles = np.linspace(0, 2 * np.pi, len(phi_vals_p)**2)
    
#     pattern = intensity_data * np.sin(theta_vals_p)


#     # Plot the radiation pattern
#     plt.figure(figsize=(8, 6))
#     ax = plt.subplot(111, polar=True)
#     ax.plot(angles, pattern.flatten(), 'b', label=label, linewidth=1)
#     ax.set_xlabel(r'$\theta$ (degrees)', fontsize=18, labelpad=8)  
#     # Remove the y-axis label to prevent overlap
#     ax.set_ylabel(r'$I(\theta, \phi)$', fontsize=18, labelpad=28)
#     ax.set_ylim([0, 1])
#     ax.tick_params(axis='x', labelsize=16, direction='in', length=6)
#     ax.tick_params(axis='y', labelsize=16, direction='in', length=6)
#     # ax.text(0, np.max(pattern), label, size=18, ha='center')
#     # Position the legend outside the plot
#     ax.legend(fontsize=16, frameon=False, bbox_to_anchor=(1.1, 1.1))
#     ax.grid(True)

#     # plt.savefig(f"{filename}.png", dpi=1000, bbox_inches='tight')
#     plt.show()

# # # Example usage (assuming intensity data is provided)
# plot_intensity_distribution(norm_int, 'Anti-parallel', 'dist_int_oppo')






























# # Load g2 data
# H_dimer_g2 = np.load('files/g2/H_dimer/H_weakC_T0.npy')
# J_dimer_g2 = np.load('files/g2/J_dimer/J_weakC_T0.npy')
# ortho_dimer_g2 = np.load('files/g2/ortho_dimer/ortho_weakC_T0.npy')
# # anti_parallel_g2 = np.load('files/g2/oppo_dimer/oppo_weakC_T0.npy')


# #Concatenate times for symmetry
# times_ns = np.linspace(0, 1.0e-3 * tf, 200000)
# times = np.linspace(0, 1.0e-3 * tf, 200000)
# times_sym = np.append(-np.flip(times_ns[1:] - times_ns[0]), times_ns)
# # Create main plot
# fig, ax = plt.subplots(figsize=(12, 8))

# # Plot the data with appropriate labels
# # line1, = ax.plot(times_sym, anti_parallel_g2_sym, '-k', label='Anti-parallel', marker='s', markersize=6, markevery=25000, linewidth=1)




# # Concatenate g2 data for symmetry
# H_dimer_g2_sym = np.append(np.flip(H_dimer_g2[1:]), H_dimer_g2)
# line2, = ax.plot(times_sym, H_dimer_g2_sym, '-b', label='H dimer', linewidth=1)

# J_dimer_g2_sym = np.append(np.flip(J_dimer_g2[1:]), J_dimer_g2)
# line3, = ax.plot(times_sym, J_dimer_g2_sym, '--r', label='J dimer', linewidth=1)

# ortho_dimer_g2_sym = np.append(np.flip(ortho_dimer_g2[1:]), ortho_dimer_g2)
# line4, = ax.plot(times_sym, ortho_dimer_g2_sym, '--g', label='Orthogonal', linewidth=1)
# # anti_parallel_g2_sym = np.append(np.flip(anti_parallel_g2[1:]), anti_parallel_g2)

# # # Create main plot
# # fig, ax = plt.subplots(figsize=(12, 8))

# # # Plot the data with appropriate labels
# # # line1, = ax.plot(times_sym, anti_parallel_g2_sym, '-k', label='Anti-parallel', marker='s', markersize=6, markevery=25000, linewidth=1)
# # line2, = ax.plot(times_sym, H_dimer_g2_sym, '-b', label='H dimer', linewidth=1)
# # line3, = ax.plot(times_sym, J_dimer_g2_sym, '--r', label='J dimer', linewidth=1)
# # line4, = ax.plot(times_sym, ortho_dimer_g2_sym, '--g', label='Orthogonal', linewidth=1)

# # Set labels and title
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
# ax.set_xlim([-1, 1])
# ax.set_ylim([0.0, 1.25])
# # ax.set_xscale('log')
# # Set tick parameters
# ax.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=30, direction='in', length=6)
# # ax.text(-6.7, 1.2, '(a)', size=30)

# # Add legend and grid
# ax.legend(handles=[line2, line3, line4], 
#           labels=['H dimer', 'J dimer', 'Orthogonal'],
#           loc='best', 
#           fontsize=30, 
#           frameon=False)

# # Change number of ticks on y-axis
# ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

# # Inset plot
# ax_inset = plt.axes([0.64, 0.3, 0.25, 0.25])  # [left, bottom, width, height]

# ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), H_dimer_g2_sym, '-b', linewidth=1)
# ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), J_dimer_g2_sym, '-r', linewidth=1)
# ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), ortho_dimer_g2_sym, '--g', linewidth=1)
# ax_inset.set_xlim([-0.005, 0.005])
# ax_inset.set_ylim([0.0, 1.25])
# ax_inset.set_xlabel(r'$\tau$ (ps)', fontsize=27)
# ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=27)
# ax_inset.tick_params(axis='x', labelsize=25, direction='in', length=6)
# ax_inset.tick_params(axis='y', labelsize=25, direction='in', length=6)
# ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))


# # Save the plot
# plt.savefig('orient_g2_weakC_T0.png', dpi=500)

# # Show the plot
# plt.show()


# # Load g2 data
# H_dimer_g2 = np.load('files/g2_con/H_dimer/H_strongC_T0_30ps.npy')
# J_dimer_g2 = np.load('files/g2_con/J_dimer/J_strongC_T0_30ps.npy')
# ortho_dimer_g2 = np.load('files/g2_con/ortho_dimer/ortho_strongC_T0_30ps.npy')


# H_dimer_g2_sym = gaussian_filter1d(np.append(np.flip(H_dimer_g2[1:]), H_dimer_g2), 120) 
# J_dimer_g2_sym = gaussian_filter1d(np.append(np.flip(J_dimer_g2[1:]), J_dimer_g2), 120)
# # ortho_dimer_g2_sym = ortho_dimer_g2

# times_ns = np.linspace(0, 1.0e-3 * tf, 20000000)
# times = np.linspace(0, 1.0e-3 * tf, 20000000)
# times_sym = np.append(-np.flip(times_ns[1:] - times_ns[0]), times_ns)

# plt.plot(times_sym, H_dimer_g2_sym)
# plt.plot(times_sym, J_dimer_g2_sym)
# plt.grid(True)
# # plt.xlim([-2, 2])
# plt.show()

# # Now you can plot without size mismatch
# fig, ax = plt.subplots(figsize=(12, 8))

# # Plot the data with appropriate labels
# line2, = ax.plot(times_sym, H_dimer_g2_sym, '-b', label='H dimer', linewidth=1)
# line3, = ax.plot(times_sym, J_dimer_g2_sym, '--r', label='J dimer', linewidth=1)
# # line4, = ax.plot(times_sym, ortho_dimer_g2_sym, '--g', label='Orthogonal', linewidth=1)

# # Set labels and title
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=35)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=35)
# ax.set_xlim([-10, 10])
# ax.set_ylim([0.0, 1.10])

# # Set tick parameters
# ax.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # Add legend and grid
# # ax.legend(handles=[line2, line3, line4], 
# #           labels=['H dimer', 'J dimer', 'Orthogonal'],
# #           loc='best', 
# #           fontsize=30, 
# #           frameon=False)
# ax.legend(handles=[line2, line3], 
#           labels=['H dimer', 'J dimer'],
#           loc='best', 
#           fontsize=30, 
#           frameon=False)

# # Change number of ticks on y-axis
# ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

# # Inset plot
# ax_inset = plt.axes([0.63, 0.3, 0.25, 0.25])  # [left, bottom, width, height]
# ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), H_dimer_g2_sym, '-b', linewidth=1)
# ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), J_dimer_g2_sym, '-r', linewidth=1)
# # ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), ortho_dimer_g2_sym, '--g', linewidth=1)
# ax_inset.set_xlim([-0.00125, 0.00125])
# ax_inset.set_ylim([0.35, 0.65])
# ax_inset.set_xlabel(r'$\tau$ (ns)', fontsize=27)
# ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=27)
# ax_inset.tick_params(axis='x', labelsize=25, direction='in', length=6)
# ax_inset.tick_params(axis='y', labelsize=25, direction='in', length=6)
# ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))

# # Save the plot
# plt.savefig('orient_g2_con_strongC_T0_140ps.png', dpi=500)

# # Show the plot
# plt.show()





























# # Load data
# J_weakC_T0 = np.load('files/intensity/J_dimer/J_weakC_T0.npy')
# J_noC_T0 = np.load('files/intensity/J_dimer/J_noC_T0.npy')
# H_weakC_T0 = np.load('files/intensity/H_dimer/H_weakC_T0.npy')
# H_noC_T0 = np.load('files/intensity/H_dimer/H_noC_T0.npy')

# # Calculate the difference in intensities
# J_intensity_diff_T0 = J_weakC_T0 - J_noC_T0
# H_intensity_diff_T0 = H_weakC_T0 - H_noC_T0

# # Create a figure with two subplots stacked vertically
# fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(14, 12), sharex=True)  # Share x-axis

# # Plot for H dimer
# ax1.set_ylabel('Population', fontsize=35)
# ax1.plot(system_params.gamma_q1 * times, np.load('H_weakC_bright_T0.npy'), label='Bright state',
#          linewidth=2.5, linestyle='-', color='blue')
# ax1.plot(system_params.gamma_q1 * times, np.load('H_weakC_dark_T0.npy'), label='Dark state',
#          linewidth=2.5, linestyle='--', color='tomato')
# ax1.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax1.tick_params(axis='y', labelsize=30, direction='in', length=6)

# ax2 = ax1.twinx()  # Secondary y-axis for intensity difference
# ax2.set_ylabel(r'Difference $\Delta I/I_{0}$', fontsize=35)
# ax2.plot(system_params.gamma_q1 * times, H_intensity_diff_T0, label=r'Difference $\Delta I/I_{0}$',
#          linewidth=2.5, linestyle=':', color='black')
# ax2.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # Add (a) label
# ax1.text(0.02, 0.9, '(a)', transform=ax1.transAxes, fontsize=35, verticalalignment='top')

# # Plot for J dimer
# ax3.set_ylabel('Population', fontsize=35)
# ax3.plot(system_params.gamma_q1 * times, np.load('J_weakC_bright_T0.npy'), label='Bright state',
#          linewidth=2.5, linestyle='-', color='blue')
# ax3.plot(system_params.gamma_q1 * times, np.load('J_weakC_dark_T0.npy'), label='Dark state',
#          linewidth=2.5, linestyle='--', color='tomato')
# ax3.tick_params(axis='x', labelsize=30, direction='in', length=6)
# ax3.tick_params(axis='y', labelsize=30, direction='in', length=6)
# ax3.set_xlabel(r'$\gamma t$', fontsize=35)

# ax4 = ax3.twinx()  # Secondary y-axis for intensity difference
# ax4.set_ylabel(r'Difference $\Delta I/I_{0}$', fontsize=35)
# ax4.plot(system_params.gamma_q1 * times, J_intensity_diff_T0, label=r'Difference $\Delta I/I_{0}$',
#          linewidth=2.5, linestyle=':', color='black')
# ax4.tick_params(axis='y', labelsize=30, direction='in', length=6)

# # Add (b) label
# ax3.text(0.02, 0.9, '(b)', transform=ax3.transAxes, fontsize=35, verticalalignment='top')

# # Combine legends from both plots and display only once
# lines_labels = [ax1.get_lines()[0], ax1.get_lines()[1], ax2.get_lines()[0]]
# labels = [line.get_label() for line in lines_labels]
# fig.legend(lines_labels, labels, loc='center left', bbox_to_anchor=(0.5, 0.7), fontsize=30, frameon=False)

# # Adjust layout to ensure text isn't cut off
# fig.tight_layout(pad=3.0)

# # Save and show the combined figure
# plt.savefig(r'combined_H_J_pop_int_weak_T0.png', dpi=500, bbox_inches='tight')
# plt.show()




# n_theta_p = 20
# n_phi_p = 20

# # Define parameters
# theta_vals_p = np.linspace(0, np.pi, n_theta_p)
# phi_vals_p = np.linspace(0,  2*np.pi, n_phi_p)
# r_p = 1  # radius of sphere

# # Generate x, y, z coordinates of points on the sphere
# x = r_p * np.outer(np.cos(phi_vals_p), np.sin(theta_vals_p))
# y = r_p * np.outer(np.sin(phi_vals_p), np.sin(theta_vals_p))
# z = r_p * np.outer(np.ones(np.size(phi_vals_p)), np.cos(theta_vals_p))


# # Create a Bloch-Redfield calculator instance
# BR_calculator = BlochRedfieldCalculator(H, a_ops_opt, Gamma_opt, rates, a_ops_vib = a_ops_vib, Gamma_vib = Gamma_vib, a_ops_coup = a_ops_coup, Gamma_coup = Gamma_coup, P_vib = P_vib, mode=mode, pump = pump)

# # Initialize the DistributionAnalysis class for calculating intensity
# quantum_analysis_intensity_dist = DistributionAnalysis(BR_calculator, rates,quantity='intensity distribution', rho=rho0, final_time=tf, steps=steps)

# # Compute intensity distribution
# norm_int = np.zeros((n_phi_p, n_theta_p))
# for l, phi_p in enumerate(phi_vals_p):
#     for k, theta_p in enumerate(theta_vals_p):
#         # theta_p = np.pi / 2
#         norm_int[l, k] = quantum_analysis_intensity_dist.compute(theta_p, phi_p)

# print(norm_int.max())
# # plt.plot(system_params.gamma_q1 * times, np.mean(norm_int, axis = (0, 1)))
# # plt.show()









# # Function to plot data
# def plot_g2(filenames, save_name, colors, linestyles):
#     fig, ax = plt.subplots(figsize=(14, 10))

#     for filename, color, linestyle in zip(filenames, colors, linestyles):
#         # Load data
#         g2_data = np.append(
#             np.flip(np.load(f'files/g2/{filename}.npy')[1:]),
#             np.load(f'files/g2/{filename}.npy')
#         )

#         times = np.append(
#             -np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]),
#             np.linspace(0, 1.0e-3 * tf, 20000)
#         )

#         # Plot the dataset
#         ax.plot(times, g2_data, linestyle, color=color, linewidth=3, label=filename)

#     # Set axis limits
#     ax.set_xlim([-10, 10])
#     ax.set_ylim([0, 1])
#     ax.tick_params(axis='x', labelsize=45, direction='in', length=6)
#     ax.tick_params(axis='y', labelsize=45, direction='in', length=6)

#     # Set axis labels with very large font size
#     ax.set_xlabel(r'Time Delay $\tau$ (ns)', fontsize=55, labelpad=20)
#     ax.set_ylabel(r'Photon Coincidence $g^{(2)}$', fontsize=55, labelpad=20)

#     # # Add legend
#     # ax.legend(fontsize=30)

#     # Save the plot
#     plt.savefig(save_name, dpi=300, bbox_inches='tight')
#     plt.close()

# # Define file names, output names, and colors
# file_names = ['g2_sing', 'g2_dist', 'g2_coop']
# save_names = ['plot_1.png', 'plot_2.png', 'plot_3.png']
# colors = ['green', 'red', 'blue']
# linestyles = ['--', '-','-' ]

# # Generate the three plots incrementally
# for i in range(1, 4):
#     plot_g2(file_names[:i], save_names[i-1], colors[:i], linestyles[:i])





# # Plot the downsampled data
# fig, ax = plt.subplots(figsize=(16, 12)) 

# # Load data
# h_dimer = np.append(
#     np.flip(np.load('files/g2/H_strongC_2nm_symm.npy')[1:]),
#     np.load('files/g2/H_strongC_2nm_symm.npy')
# )
# j_dimer = np.append(
#     np.flip(np.load('files/g2/J_strongC_2nm_symm.npy')[1:]),
#     np.load('files/g2/J_strongC_2nm_symm.npy')
# )
# times = np.append(
#     -np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]),
#     np.linspace(0, 1.0e-3 * tf, 100000)
# )

# # Main plot
# line0, = ax.plot(times, h_dimer, '-b', linewidth=3, label='H dimer')
# # line1, = ax.plot(times, j_dimer, '-r', linewidth=3, label='J dimer')

# # Set axis labels and limits with large font size
# ax.set_xlabel(r'Time Delay $\tau$ (ns)', fontsize=55, labelpad=20)
# ax.set_ylabel(r'Photon Coincidence $g^{(2)}$', fontsize=55, labelpad=20)
# ax.set_xlim([-10, 10])
# ax.set_ylim([0.4, 1.3])
# ax.tick_params(axis='x', labelsize=40, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=40, direction='in', length=6)
# ax.legend(handles=[line0], loc='upper center', fontsize=40, frameon=False, ncol=2, bbox_to_anchor=(0.15, 0.95), handlelength=1.5)

# # Save and show the plot
# plt.savefig('plot_H_dimer.png', dpi=300, bbox_inches='tight')
# plt.show()

# fig, ax = plt.subplots(figsize=(16, 12)) 
# # Main plot
# line0, = ax.plot(times, h_dimer, '-b', linewidth=3, label='H dimer')
# line1, = ax.plot(times, j_dimer, '-r', linewidth=3, label='J dimer')

# # Set axis labels and limits with large font size
# ax.set_xlabel(r'Time Delay $\tau$ (ns)', fontsize=55, labelpad=20)
# ax.set_ylabel(r'Photon Coincidence $g^{(2)}$', fontsize=55, labelpad=20)
# ax.set_xlim([-10, 10])
# ax.set_ylim([0.4, 1.3])
# ax.tick_params(axis='x', labelsize=40, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=40, direction='in', length=6)
# # Create the first legend for H dimer at the same position as in the first plot:
# legend_h = ax.legend([line0], ['H dimer'], loc='upper center', fontsize=40,
#                        frameon=False, handlelength=1.5, bbox_to_anchor=(0.15, 0.95))
# ax.add_artist(legend_h)  # Ensure it is added to the axes

# legend_j = ax.legend([line1], ['J dimer'], loc='upper right', fontsize=40,
#                        frameon=False, bbox_to_anchor=(0.85, 0.95), handlelength=1.5)
# ax.add_artist(legend_j) 
# # Save and show the plot
# plt.savefig('plot_HJ_dimer.png', dpi=300, bbox_inches='tight')
# plt.show()






# sigma = 40
# zero_delay_H_convolved_array = []
# zero_delay_J_convolved_array = []
# step_con = 5.0e5

# dt = (5 * system_params.tau_L - (-5 * system_params.tau_L)) / step_con
# g2_corr_H_appended = gaussian_filter1d(h_dimer, sigma=sigma / dt, mode='reflect')
# g2_corr_J_appended = gaussian_filter1d(j_dimer, sigma=sigma / dt, mode='reflect')
# times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)


# # Main plot
# line0, = ax.plot(times_appended, g2_corr_H_appended, '-b', linewidth=3)
# line1, = ax.plot(times_appended, g2_corr_J_appended, '-r', linewidth=3)

# # Set axis labels and limits
# ax.set_xlabel(r'Time Delay $\tau$ (ns)', fontsize=55, labelpad=20)
# ax.set_ylabel(r'Photon Coincidence $g^{(2)}$', fontsize=55, labelpad=20)
# ax.set_xlim([-10, 10])
# ax.set_ylim([0.4, 1.3])
# ax.tick_params(axis='x', labelsize=40, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=40, direction='in', length=6)
# # Create the first legend for H dimer at the same position as in the first plot:
# legend_h = ax.legend([line0], ['H dimer'], loc='upper center', fontsize=40,
#                        frameon=False, handlelength=1.5, bbox_to_anchor=(0.15, 0.95))
# ax.add_artist(legend_h)  # Ensure it is added to the axes

# legend_j = ax.legend([line1], ['J dimer'], loc='upper right', fontsize=40,
#                        frameon=False, bbox_to_anchor=(0.85, 0.95), handlelength=1.5)
# ax.add_artist(legend_j) 

# # plt.plot(times_appended, g2_corr_H_appended)
# # plt.plot(times_appended, g2_corr_J_appended)
# plt.savefig('orient_g2_symm_convolved', dpi=300, bbox_inches='tight')
# plt.show()
#





# # Set a common font size for all text
# fs = 30  # change this value as needed

# # Define instrument responses
# instrument_responses = np.linspace(1, 300, 10000)

# # Load data
# data_H_nfd_slow = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_nojigg_convolved_smart.npy')
# data_H_nfd_slow_jigg_k10 = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k10.npy')
# data_H_nfd_fast_jigg_k5 = np.load('files/instrument_responses/H_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy')

# data_J_nfd_slow = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_nojigg_convolved_smart.npy')
# data_J_nfd_slow_jigg_k10 = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k10.npy')
# data_J_nfd_fast_jigg_k5 = np.load('files/instrument_responses/J_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy')

# # Colors and line styles
# colors_H = ['lightsteelblue', 'royalblue', 'mediumblue', 'midnightblue']
# colors_J = ['rosybrown', 'lightcoral', 'red', 'maroon']
# line_styles = ['-', '--', (0, (3, 1, 1, 1)), ':']

# # Create a figure with two vertically stacked subplots sharing the x-axis
# fig, (ax_H, ax_J) = plt.subplots(2, 1, figsize=(9.5, 13), sharex=True,
#                                  gridspec_kw={'hspace': 0.05, 'height_ratios': [1, 1]})

# # ---------------------------
# # Plot for H dimer
# # ---------------------------
# ax_H.plot(instrument_responses, data_H_nfd_slow, color=colors_H[0], linestyle=line_styles[0],
#           label=r'$[g^{(2)}(\infty, 0)]$', linewidth=3)
# ax_H.plot(instrument_responses, data_H_nfd_slow_jigg_k10, color=colors_H[1], linestyle=line_styles[1],
#           label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 10}$', linewidth=3)
# ax_H.plot(instrument_responses, data_H_nfd_fast_jigg_k5, color=colors_H[3], linestyle=line_styles[3],
#           label=r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', linewidth=3)

# # Set y-axis limits and ticks for H dimer plot
# ax_H.set_ylim([0.3, 1.3])
# ax_H.set_yticks([0.4, 0.7, 1.0, 1.3])

# # Apply tick parameters using the common font size
# ax_H.tick_params(axis='x', labelsize=fs, direction='in', length=6)
# ax_H.tick_params(axis='y', labelsize=fs, direction='in', length=6)

# # Add legend to H dimer plot
# ax_H.legend(loc='lower center', bbox_to_anchor=(0.6, 0.45), fontsize=fs,
#             frameon=False, handletextpad=1, borderaxespad=0.5, labelspacing=1.5, ncol=1)

# # ---------------------------
# # Plot for J dimer
# # ---------------------------
# ax_J.plot(instrument_responses, data_J_nfd_slow, color=colors_J[0], linestyle=line_styles[0],
#           label=r'$[g^{(2)}(\infty, 0)]$', linewidth=3)
# ax_J.plot(instrument_responses, data_J_nfd_slow_jigg_k10, color=colors_J[1], linestyle=line_styles[1],
#           label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 10}$', linewidth=3)
# ax_J.plot(instrument_responses, data_J_nfd_fast_jigg_k5, color=colors_J[3], linestyle=line_styles[3],
#           label=r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', linewidth=3)

# # Set y-axis limits and ticks for J dimer plot
# ax_J.set_ylim([0.48, 0.9])
# ax_J.set_yticks([0.55, 0.65, 0.75, 0.85])

# # Set x-axis label using the common font size
# ax_J.set_xlabel(r'Instrument response $\Delta \tau$ (ps)', fontsize=fs, labelpad=0)

# # Apply tick parameters to J dimer plot
# ax_J.tick_params(axis='x', labelsize=fs, direction='in', length=6)
# ax_J.tick_params(axis='y', labelsize=fs, direction='in', length=6)

# # Add legend to J dimer plot
# ax_J.legend(loc='lower center', bbox_to_anchor=(0.6, 0.45), fontsize=fs,
#             frameon=False, handletextpad=1, borderaxespad=0.5, labelspacing=1.5, ncol=1)

# # ---------------------------
# # Shared y-axis label for both plots
# # ---------------------------
# fig.text(0.0, 0.5, r'Zero delay coincidence', va='center', rotation='vertical', fontsize=fs)


# # Adjust subplot layout so that labels are not clipped
# plt.subplots_adjust(left=0.12, right=0.96, top=0.99, bottom=0.07)

# # Save the figure as an SVG (with bbox_inches='tight' to include all text)
# plt.savefig('zero_time_g2_H_J_3.svg', format='svg', dpi=300, bbox_inches='tight')
# plt.show()


















# # Plot the downsampled data
# fig, ax = plt.subplots(figsize=(14, 10))

# # # # Plot the data
# # # # line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')), '-b')
# # # # line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')), '-r')

# # line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/H_strongC_2nm_symm.npy')[1:]), np.load('files/g2/H_strongC_2nm_symm.npy')), '-b', linewidth = 3)
# # line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/J_strongC_2nm_symm.npy')[1:]), np.load('files/g2/J_strongC_2nm_symm.npy')), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2_con/ortho_strongC_2nm_symm.npy')[1:]), np.load('files/g2_con/ortho_strongC_2nm_symm.npy')), '--g', linewidth = 3)

# # # line0, = ax.plot(np.linspace(0, 1.0e-3 * tf, 100000), np.load('files/g2/H_dimer/H_strongC_2nm_symm.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(np.linspace(0, 1.0e-3 * tf, 100000), np.load('files/g2/J_dimer/J_strongC_2nm_symm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax.plot(np.linspace(0, 1.0e-3 * tf, 100000), np.load('files/g2/ortho_strongC_2nm_symm.npy'), '--g', linewidth = 3)

# line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_perp.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_perp.npy')), '-b', linewidth = 3)
# line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_par.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_par.npy')), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# line2, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_intermediate.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_intermediate.npy')), '-g', linewidth = 3)
# line3, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/ortho_strongC_2nm_symm_intermediate_neg.npy')[1:]), np.load('files/g2/ortho_strongC_2nm_symm_intermediate_neg.npy')), '-k', linewidth = 3)


# # # # line0, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_weakC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_weakC_T0_50nm_wcme.npy')), '-b')
# # # # line1, = ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/ortho_dimer/ortho_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/ortho_dimer/ortho_strongerC_T0_50nm_wcme.npy')), '-r')


# # # line0, = ax.plot(times_appended, np.load('files/g2_con/H_dimer/H_strongC_2nm_symm.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(times_appended, np.load('files/g2_con/J_dimer/J_strongC_2nm_symm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax.plot(times_appended, np.load('files/g2_con/ortho_strongC_2nm_symm.npy'), '--g', linewidth = 3)


# # # line0, = ax.plot(times_ns, np.load('files/g2/crystal/H_dimer/H_weakC_2nm_symm_unnorm.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(times_ns, np.load('files/g2/crystal/H_dimer/H_strongC_2nm_symm_fd_unnorm.npy'), '-r', linewidth = 3)


# # # line0, = ax.plot(times_appended, np.load('files/g2_con/H_strongC_2nm_symm_nfd_slow_tumb_150.npy'), '-b', linewidth = 3)
# # # line1, = ax.plot(times_appended, np.load('files/g2_con/J_strongC_2nm_symm_nfd_slow_tumb_150.npy'), '-r', linewidth = 3)

# # Set axis labels and limits
# ax.set_xlabel(r'$\tau$ (ns)', fontsize=40)
# ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=40)
# # ax.set_xlim([-1.0e-2, 1.0e-2])  
# ax.set_xlim([-40, 40]) 
# # ax.set_ylim([0.433, 0.436]) 
# ax.tick_params(axis='x', labelsize=35, direction='in', length=6)
# ax.tick_params(axis='y', labelsize=35, direction='in', length=6)

# # ax.set_xscale('log')
# # Add a legend with appropriate styling
# # ax.legend(handles=[line0], labels=[r"$\theta = 0$"], loc='best', fontsize=30, frameon=False)
# # ax.legend(handles=[line0, line1], labels=[r"$\theta = 0$", r"$\theta = \pi/2, \phi = n\pi/2$" ], loc='best', fontsize=30, frameon=False)
# # ax.legend(handles=[line0, line1, line2], labels=[r"$\theta = 0$", r"$\theta = \pi/2, \phi = n\pi/2$" , r"$\theta = \pi/2$, $\phi = \pi/4$"], loc='best', fontsize=30, frameon=False)
# ax.legend(handles=[line0, line1, line2, line3], labels=[r"$\theta = 0$", r"$\theta = \pi/2, \phi = n\pi/2$" , r"$\theta = \pi/2$, $\phi = \pi/4$", r"$\theta = \pi/2$, $\phi = 3\pi/4$"], loc='best', fontsize=30, frameon=False)
# # ax.legend(handles=[line0, line1], labels=[r"H dimer", r"J dimer"], fontsize=30, frameon=False)

# # # ax.text(-0.1, 1.05, r'(b)', size=30, transform=ax.transAxes)  # Adjusted position to upper-left

# # # Inset plot
# # ax_inset = plt.axes([0.63, 0.21, 0.24, 0.25])  # [left, bottom, width, height]

# # # # ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_weakC_T0_50nm_wcme.npy')), '-b')
# # # # ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/J_dimer/J_strongerC_T0_50nm_wcme.npy')), '-r')

# # line0, = ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/H_strongC_2nm_symm.npy')[1:]), np.load('files/g2/H_strongC_2nm_symm.npy')), '-b', linewidth = 3)
# # line1, = ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2/J_strongC_2nm_symm.npy')[1:]), np.load('files/g2/J_strongC_2nm_symm.npy')), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # line2, = ax_inset.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 100000)[1:] - np.linspace(0, 1.0e-3 * tf, 100000)[0]), np.linspace(0, 1.0e-3 * tf, 100000)), np.append(np.flip(np.load('files/g2_con/ortho_strongC_2nm_symm.npy')[1:]), np.load('files/g2_con/ortho_strongC_2nm_symm.npy')), '--g', linewidth = 3)

# # # # line0, = ax_inset.plot(times_appended, np.load('files/g2_con/H_dimer/H_strongC_2nm_symm.npy'), '-b', linewidth = 3)
# # # # line1, = ax_inset.plot(times_appended, np.load('files/g2_con/J_dimer/J_strongC_2nm_symm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')
# # # # line2, = ax_inset.plot(times_appended, np.load('files/g2_con/ortho_strongC_2nm_symm.npy'), '--g', linewidth = 3)


# # # # line0, = ax_inset.plot(times_appended, np.load('files/g2_con/H_dimer/H_strongC_2nm.npy'), '-b', linewidth = 3)
# # # # line1, = ax_inset.plot(times_appended, np.load('files/g2_con/J_dimer/J_strongC_2nm.npy'), '-r', linewidth = 3) #ax.plot(np.append(-np.flip(np.linspace(0, 1.0e-3 * tf, 20000)[1:] - np.linspace(0, 1.0e-3 * tf, 20000)[0]), np.linspace(0, 1.0e-3 * tf, 20000)), np.append(np.flip(np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')[1:]), np.load('files/g2/H_dimer/H_strongerC_T0_50nm_wcme.npy')), '-r')

 
# # ax_inset.set_xlim([-0.25, 0.25])
# # ax_inset.set_ylim([0.4, 1.25])
# # ax_inset.set_xlabel(r'$\tau$ (ns)', fontsize=27)
# # ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=27)
# # ax_inset.tick_params(axis='x', labelsize=25, direction='in', length=6)
# # ax_inset.tick_params(axis='y', labelsize=25, direction='in', length=6)
# # ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))


# plt.savefig('orient_g2_ortho_4', dpi=150, bbox_inches='tight')
# plt.show()




















# Define instrument responses
instrument_responses = np.linspace(1, 300, 10000)

# Load data

# Fast tumbling without any wiggling would give the same result as the static case because, the photon correlation doesn't change based on the direction of measurement for parallel dipoles.


data_H_nfd_slow = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_nojigg_convolved_smart.npy')
data_H_nfd_slow_jigg_k5 = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k5.npy')
data_H_nfd_slow_jigg_k10 = np.load('files/instrument_responses/H_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k10.npy')
data_H_nfd_fast_jigg_k5 = np.load('files/instrument_responses/H_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy')
# data_H_fd_slow = np.load('files/H_strongC_2nm_symm_fd_slow_tumb_convolved.npy')

data_J_nfd_slow = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_nojigg_convolved_smart.npy')
data_J_nfd_slow_jigg_k5 = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k5.npy')
data_J_nfd_slow_jigg_k10 = np.load('files/instrument_responses/J_strongC_2nm_symm_slow_tumb_jigg_convolved_smart_k10.npy')
data_J_nfd_fast_jigg_k5 = np.load('files/instrument_responses/J_strongC_2nm_symm_fast_tumb_jigg_convolved_smart_k5.npy')
# data_J_fd_slow = np.load('files/J_weakC_2nm_symm_fd_slow_tumb_convolved.npy')

# Colors for H dimer (shades of blue)
colors_H = ['lightsteelblue', 'royalblue', 'mediumblue', 'midnightblue']  # Variations of blue
colors_J = ['rosybrown', 'lightcoral', 'red', 'maroon']  # Variations of red

# Updated line styles for more differentiation
line_styles = ['-', '--', (0, (3, 1, 1, 1)), ':']  # Solid, dashed, dash-dot, dotted

# Load g2 data
g2_corr_J = np.load('files/g2/J_strongC_2nm_symm_slow_tumb_jigg_normalized_smart_k5.npy')
g2_corr_H = np.load('files/g2/H_strongC_2nm_symm_slow_tumb_jigg_normalized_smart_k5.npy')

# Inset configurations
insets_H = [
    {'ir': 50, 'pos': [0.2, 0.7, 0.21, 0.24]},  # x, y, width, height
    {'ir': 100, 'pos': [0.09, 0.14, 0.21, 0.24]},
    {'ir': 200, 'pos': [0.42, 0.08, 0.21, 0.24]}
]

insets_J = [
    {'ir': 50, 'pos': [0.18, 0.72, 0.21, 0.24]},
    {'ir': 100, 'pos': [0.09, 0.1, 0.21, 0.24]},
    {'ir': 200, 'pos': [0.42, 0.07, 0.21, 0.24]}
]


def add_inset(ax, g2_corr, instrument_response, inset_pos, x_lim, y_lim, color=None):
    step_con = 1.0e4
    dt = (5 * system_params.tau_L - (-5 * system_params.tau_L)) / step_con
    sigma = instrument_response

    g2_corr_appended = gaussian_filter1d(
        np.append(np.flip(g2_corr[1:]), g2_corr), 
        sigma=sigma / dt, 
        mode='reflect'
    )
    times_appended = np.append(-np.flip((times_ns[1:] - times_ns[0])), times_ns)

    inset_ax = ax.inset_axes(inset_pos)
    inset_ax.plot(times_appended, g2_corr_appended, color=color)
    inset_ax.set_xlim(x_lim)
    inset_ax.set_ylim(y_lim)

    inset_ax.set_xlabel(r'$\tau$ (ns)', fontsize=18)
    inset_ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=18, labelpad=-25)
    # Add y-ticks and two numerical labels
    inset_ax.set_yticks([y_lim[0], y_lim[1]])  
    inset_ax.set_yticklabels([f'{y_lim[0]:.2f}', f'{y_lim[1]:.2f}'], fontsize=18)

    # Hide x-ticks
    inset_ax.set_xticklabels([])
    inset_ax.tick_params(axis='x', which='both', length=0, labelbottom=False)

    inset_ax.text(0.5, 0.95, rf'$\Delta \tau = {instrument_response}\,\mathrm{{ps}}$', 
              ha='center', va='top', transform=inset_ax.transAxes, fontsize=18)



# Create a single figure with two subplots, stacked vertically, sharing the x-axis
fig = plt.figure(figsize=(20, 6.5))  # Adjust size as needed
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.1)

ax_H = fig.add_subplot(gs[0])
ax_J = fig.add_subplot(gs[1])

# Plot for H dimer
ax_H.plot(instrument_responses, data_H_nfd_slow, color=colors_H[0], linestyle=line_styles[0], label=r'$[g^{(2)}(\infty, 0)]$', linewidth=3)
ax_H.plot(instrument_responses, data_H_nfd_slow_jigg_k10, color=colors_H[1], linestyle=line_styles[1], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 10}$', linewidth=3)
ax_H.plot(instrument_responses, data_H_nfd_slow_jigg_k5, color=colors_H[2], linestyle=line_styles[2], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 5}$', linewidth=3)
ax_H.plot(instrument_responses, data_H_nfd_fast_jigg_k5, color=colors_H[3], linestyle=line_styles[3], label=r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', linewidth=3)

# Set title and legend for H dimer
ax_H.legend(
    loc='lower center', 
    bbox_to_anchor=(0.65, 0.45),  # Centered horizontally, slightly above the plot
    fontsize=20, 
    frameon=False, 
    handletextpad=1, 
    borderaxespad=0.5, 
    labelspacing=1.5, 
    ncol=1 
)
ax_H.tick_params(axis='x', labelsize=20, direction='in', length=6)
ax_H.tick_params(axis='y', labelsize=20, direction='in', length=6)
ax_H.set_ylim([0.3, 1.35])

# Add insets for H dimer
for inset in insets_H:
    add_inset(ax_H, g2_corr_H, inset['ir'], inset['pos'], [-2, 2], [0.43, 0.9], color=colors_H[3])


# Plot for J dimer
ax_J.plot(instrument_responses, data_J_nfd_slow, color=colors_J[0], linestyle=line_styles[0], label=r'$[g^{(2)}(\infty, 0)]$', linewidth=3)
ax_J.plot(instrument_responses, data_J_nfd_slow_jigg_k10, color=colors_J[1], linestyle=line_styles[1], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 10}$', linewidth=3)
ax_J.plot(instrument_responses, data_J_nfd_slow_jigg_k5, color=colors_J[2], linestyle=line_styles[2], label=r'$[g^{(2)}(\infty, 0)]_{O}^{\kappa = 5}$', linewidth=3)
ax_J.plot(instrument_responses, data_J_nfd_fast_jigg_k5, color=colors_J[3], linestyle=line_styles[3], label=r'$[g^{(2)}(\infty, 0)]_{(O,\textbf{q})}^{\kappa = 5}$', linewidth=3)

# Set title and legend for J dimer
ax_J.legend(
    loc='lower center', 
    bbox_to_anchor=(0.6, 0.45),  # Centered horizontally, slightly above the plot
    fontsize=20, 
    frameon=False, 
    handletextpad=1, 
    borderaxespad=0.5, 
    labelspacing=1.5, 
    ncol=1
)
fig.text(0.5, 0.01, r'Instrument response $\Delta \tau$ (ps)', ha='center', fontsize=25)

ax_J.tick_params(axis='x', labelsize=20, direction='in', length=6)
ax_J.tick_params(axis='y', labelsize=20, direction='in', length=6)
ax_J.set_ylim([0.42, 0.95])

# Add insets for J dimer
for inset in insets_J:
    add_inset(ax_J, g2_corr_J, inset['ir'], inset['pos'], [-2, 2], [0.55, 0.71], color=colors_J[3])

# Add a shared y-axis label between the plots
fig.text(0.085, 0.5, r'Zero delay coincidence', va='center', rotation='vertical', fontsize=25)
plt.subplots_adjust(left=0.12, right=0.96, top=0.99, bottom=0.1)

ax_H.set_yticks([0.4, 0.7, 1.0, 1.3])  # For H dimer plot
ax_J.set_yticks([0.4, 0.6, 0.8, 1.0])  # For J dimer plot

# Save the combined plot
plt.savefig('zero_time_g2_H_J', dpi=300, bbox_inches='tight')
# plt.show()
