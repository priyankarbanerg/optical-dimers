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

# # # Plot the intensity as a function of time
# # plt.figure(figsize=(8, 6))
# # plt.plot(system_params.gamma_q1 * times, intensity_results)
# # plt.plot(system_params.gamma_q1 * times, np.load('files/intensity/ortho_dimer/ortho_strongC.npy'))
# # # plt.plot(system_params.gamma_q1 * times, np.load('files/H_strongC.npy') - np.load('files/J_strongC.npy'))
# # # plt.plot(system_params.gamma_q1 * times * times, 2 * np.exp(-system_params.gamma_q1 * times * times), label = 'theoretical')
# # # plt.plot(gamma_q1*times, 2*np.exp(-2*gamma_q1*times)*(1 + 2*gamma_q1*times), label='theoretical')
# # plt.xlabel(r'$\gamma t$', fontsize = 14)
# # plt.ylabel(r'$I(t) / I_0$', fontsize = 14)
# # # plt.semilogx()
# # # plt.xlim([0, 5])  
# # # plt.ylim([0, 2.1])
# # plt.title('Intensity vs time', fontsize = 18)
# # plt.legend(loc='best')
# # plt.grid(True)
# # plt.show()

# # # Plot the intensity as a function of time
# # plt.figure(figsize=(8, 6))
# # # plt.plot(times, pop_results, label = 'H_noC_dark')
# # plt.plot(times, np.load('files/datafiles/H_strongC_bright.npy'), label = 'H_strongC_bright')
# # plt.plot(times, np.load('files/datafiles/H_strongC_dark.npy'), label = 'H_strongC_dark')
# # plt.plot(times, np.load('files/datafiles/J_strongC_bright.npy'), label = 'J_strongC_bright')
# # plt.plot(times, np.load('files/datafiles/J_strongC_dark.npy'), label = 'J_strongC_dark')
# # # plt.plot(gamma_q1 * times, 2 * np.exp(-gamma_q1 * times), label = 'theoretical')
# # # plt.plot(gamma_q1*times, 2*np.exp(-2*gamma_q1*times)*(1 + 2*gamma_q1*times), label='theoretical')
# # plt.xlabel(r'$\gamma t$', fontsize = 14)
# # plt.ylabel(r'$I(t) / I_0$', fontsize = 14)
# # # plt.semilogx()
# # # plt.xlim([0, 5])  
# # plt.ylim([-0.1, 1])
# # plt.title('Intensity vs time', fontsize = 18)
# # plt.legend(loc='best')
# # plt.grid(True)
# # plt.show()



# Read the CSV file
df = pd.read_csv('files/dipole_positions/dipole_pos_dir_oppo.csv')

# Remove leading spaces from column names
df.columns = df.columns.str.strip()


# # Normalize the intensity values
norm_int = norm_int / norm_int.max()

# Create figure
fig = go.Figure()

# Add initial surface plot with dynamic color scaling
# cmin = norm_int.min()
# cmax = norm_int.max()
fig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor=norm_int, colorscale='RdBu_r', opacity=0.7, cmin = 0.0, cmax = 1.0))

# Add dipoles (constant for all frames in this example)
# Plot dipoles
for index, row in df.iterrows():
    # Calculate end position of the arrow shaft (half length in each direction)
    half_length = 0.1  # Half length of the arrow shaft
    start_x = row['x'] - 0.20 * row['dir_x']
    start_y = row['y'] - 0.20 * row['dir_y']
    start_z = row['z'] - 0.20 * row['dir_z']
    end_x = row['x'] + 0.20 * row['dir_x']
    end_y = row['y'] + 0.20 * row['dir_y']
    end_z = row['z'] + 0.20 * row['dir_z']

    # Plot dipole position (spheres)
    fig.add_trace(go.Scatter3d(
        x=[row['x']],
        y=[row['y']],
        z=[row['z']],
        mode='markers',
        marker=dict(color='red', size=10),  # Color set to red
        name=''
    ))

    # Plot dipole direction shaft (arrows)
    fig.add_trace(go.Scatter3d(
        x=[start_x, end_x],
        y=[start_y, end_y],
        z=[start_z, end_z],
        mode='lines',
        line=dict(color='blue', width=10),  # Color set to blue
        name=''
    ))

    # Plot dipole direction arrowhead (cones)
    fig.add_trace(go.Cone(
        x=[end_x],
        y=[end_y],
        z=[end_z],
        u=[row['dir_x']],
        v=[row['dir_y']],
        w=[row['dir_z']],
        showscale=False,
        sizemode="absolute",
        sizeref=0.25,
        anchor="tail",
        opacity=1.0,
        colorscale=[[0, 'blue'], [1, 'blue']],  # Color set to blue
        name=''
    ))
# Customize axis labels and ticks
fig.update_layout(
    scene=dict(
        xaxis=dict(
            title=dict(text=r"x", font=dict(size=30)),  # X-axis label in LaTeX format
            tickfont=dict(size=18)  # X-axis tick font size
        ),
        yaxis=dict(
            title=dict(text=r"y", font=dict(size=30)),  # Y-axis label in LaTeX format
            tickfont=dict(size=18)  # Y-axis tick font size
        ),
        zaxis=dict(
            title=dict(text=r"z", font=dict(size=30)),  # Z-axis label in LaTeX format
            tickfont=dict(size=18)  # Z-axis tick font size
        ),
    )
)

# Customize colorbar
fig.update_traces(
    selector=dict(type='surface'),
    colorbar=dict(
        thickness=45,  # Colorbar width
        tickfont=dict(size=30),  # Colorbar tick font size
    )
)


# Display the animation
fig.show()



# # Create frames with dynamic color scaling
# frames = [go.Frame(
#     data=[go.Surface(z=z, surfacecolor=norm_int[:, :, k], colorscale='hot', opacity=0.7)], name=f'frame{k}') for k in range(200)]



# # # Normalize the g2 values
# # g2_corr = np.interp(np.array(g2_corr), (np.array(g2_corr).min(), np.array(g2_corr).max()), (0, 1))  #/ g2_corr.max()

# # Create figure
# fig = go.Figure()
# cmin = g2_corr.min() #0.0 #0.0 * 0.007
# cmax = g2_corr.max() #1.0 #1.0 * 0.0144


# # Add initial surface plot with dynamic color scaling
# fig.add_trace(go.Surface(
#     x=x, y=y, z=z, surfacecolor=g2_corr[:, :, 0], colorscale='hot', opacity=0.6, cmin=cmin, cmax=cmax
# ))

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
# #     name=''
# # ))




# # Create frames with dynamic color scaling
# frames = [go.Frame(
#     data=[go.Surface(z=z, surfacecolor=g2_corr[:, :, k], colorscale='hot', opacity=0.7, cmin=cmin, cmax=cmax)], name=f'frame{k}') for k in range(200)]


# # Add frames to the figure
# fig.frames = frames

# # Update layout with animation settings
# fig.update_layout(
#     scene=dict(
#         xaxis_title='X',
#         yaxis_title='Y',
#         zaxis_title='Z',
#         aspectmode='cube'
#     ),
#     title='Normalised Intensity Plot on a Sphere',
#     updatemenus=[{
#         'buttons': [
#             {
#                 'args': [None, {'frame': {'duration': 100, 'redraw': True}, 'fromcurrent': True}],
#                 'label': 'Play',
#                 'method': 'animate'
#             },
#             {
#                 'args': [[None], {'frame': {'duration': 0, 'redraw': True}, 'mode': 'immediate', 'transition': {'duration': 0}}],
#                 'label': 'Pause',
#                 'method': 'animate'
#             }
#         ],
#         'direction': 'left',
#         'pad': {'r': 10, 't': 87},
#         'showactive': False,
#         'type': 'buttons',
#         'x': 0.1,
#         'xanchor': 'right',
#         'y': 0,
#         'yanchor': 'top'
#     }]
# )
# # fig.write_html("g2_my.html")
# # Display the animation
# fig.show()
















# # # Plot the steady state second order correlation as a function of time
# # plt.figure(figsize=(8, 6))
# # plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(g2_corr[1:]), g2_corr), label = 'Measurement-induced')
# # plt.plot(np.append(- np.flip((times_ns[1:] - times_ns[0])), (times_ns)), np.append(np.flip(np.load('g2_mic_deph.npy')[1:]), np.load('g2_mic_deph.npy')), label = 'Measurement-induced - pure dephasing')
# # plt.xlabel(r'$\tau$ (ns)', fontsize = 14)
# # plt.ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize = 14)
# # # plt.xlim([-5, 5])
# # plt.ylim([0.0, 1.1])
# # plt.title(r'Photon coincidence with $\gamma_{\mathrm{pump}} = \gamma$', fontsize = 18)
# # plt.legend(loc='best')
# # plt.grid(True)
# # plt.savefig('g2-ortho-parallel.png')
# # # plt.show()













# # # Generate R_vals and Theta_vals for the polar plot
# # R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# # ax.plot(Theta_vals, norm_int, color='black', linewidth = 1)
# # ax.grid(True)
# # ax.set_title("Polar Heatmap of Intensity on a Sphere", va='bottom')
# # plt.savefig('2D polar J.png', dpi = 1000)
# # # plt.show()


# # # Generate R_vals and Theta_vals for the polar plot
# # R_vals, Theta_vals = np.meshgrid(np.sin(theta_vals_p), phi_vals_p)

# # # Create the polar heatmap
# # fig = plt.figure(layout='constrained')
# # ax1 = fig.add_subplot(projection='polar')
# # c = ax1.pcolormesh(Theta_vals, R_vals, g2_corr[:, :, 0], cmap='hot', shading='auto')

# # # Add colorbar
# # plt.colorbar(c, label='Intensity')
# # plt.title('Polar Heatmap of g(2) on a Sphere')


# # # Add arrows
# # # Arrow facing theta = 0 (along the x-axis)
# # ax1.annotate('', xy=(0, 0.5), xytext=(0, 0),
# #              arrowprops=dict(facecolor='blue', shrink=0, width=2, headwidth=10))

# # # Arrow facing theta = 90 degrees (along the y-axis)
# # ax1.annotate('', xy=(np.pi/2, 0.5), xytext=(0, 0),
# #              arrowprops=dict(facecolor='blue', shrink=0, width=2, headwidth=10))

# # # Save the figure
# # plt.savefig('2D polar - 2.png')











# # # Function to plot intensity distribution
# # def plot_intensity_distribution(intensity_data, label, filename):
   
# #     angles = np.linspace(0, 2 * np.pi, len(phi_vals_p)**2)
    
# #     pattern = intensity_data * np.sin(theta_vals_p)


# #     # Plot the radiation pattern
# #     plt.figure(figsize=(8, 6))
# #     ax = plt.subplot(111, polar=True)
# #     ax.plot(angles, pattern.flatten(), 'b', label=label, linewidth=1)
# #     ax.set_xlabel(r'$\theta$ (degrees)', fontsize=18, labelpad=8)  
# #     # Remove the y-axis label to prevent overlap
# #     ax.set_ylabel(r'$I(\theta, \phi)$', fontsize=18, labelpad=28)
# #     ax.set_ylim([0, 1])
# #     ax.tick_params(axis='x', labelsize=16, direction='in', length=6)
# #     ax.tick_params(axis='y', labelsize=16, direction='in', length=6)
# #     # ax.text(0, np.max(pattern), label, size=18, ha='center')
# #     # Position the legend outside the plot
# #     ax.legend(fontsize=16, frameon=False, bbox_to_anchor=(1.1, 1.1))
# #     ax.grid(True)

# #     # plt.savefig(f"{filename}.png", dpi=1000, bbox_inches='tight')
# #     plt.show()

# # # # Example usage (assuming intensity data is provided)
# # plot_intensity_distribution(norm_int, 'Anti-parallel', 'dist_int_oppo')





























# # H_dimer = np.load('files/intensity/H_dimer/H_strongerC.npy')
# # J_dimer = np.load('files/intensity/J_dimer/J_strongerC.npy')
# # ortho_dimer = np.load('files/intensity/ortho_dimer/ortho_strongerC.npy')
# # anti_parallel = np.load('files/intensity/oppo_dimer/oppo_strongerC.npy')

# # # Create main plot
# # fig, ax = plt.subplots(figsize=(8, 12))

# # # # Plot the data with appropriate labels
# # # ax.plot(system_params.gamma_q1 * times, intensity_results, '-y', label='Anti-parallel', marker='s', markersize=5, markevery=4300)

# # line1, = ax.plot(system_params.gamma_q1 * times, anti_parallel, '-k', label='Anti-parallel', marker='s', markersize=6, markevery=25000)
# # line2, = ax.plot(system_params.gamma_q1 * times, H_dimer, '-b', label='H dimer', marker='s', markersize=0, markevery=10)
# # line3, = ax.plot(system_params.gamma_q1 * times, J_dimer, '-r', label='J dimer', marker='^', markersize=0, markevery=13)
# # line4, = ax.plot(system_params.gamma_q1 * times, ortho_dimer, '--g', label='Orthogonal')

# # # Set labels and title
# # ax.set_xlabel(r'$\gamma t$', fontsize=18)
# # ax.set_ylabel(r'$I(t) / I_0$', fontsize=18)
# # ax.set_ylim([0, 2.1])
# # ax.set_xlim([system_params.gamma_q1 * times[0], system_params.gamma_q1 * times[-1]])
# # # Set tick parameters
# # ax.tick_params(axis='x', labelsize=16, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=16, direction='in', length=6)
# # ax.text(0.025, 2.035, '(a)', size=18)

# # # ax.set_title('Intensity vs time', fontsize=18)
# # # Change number of ticks on y-axis
# # ax.yaxis.set_major_locator(MaxNLocator(nbins=5))  # Change nbins to desired number of ticks

# # # Add legend and grid
# # ax.legend(handles=[line2, line3, line4, line1], 
# #           labels=['H dimer', 'J dimer', 'Orthogonal', 'Anti-parallel'],
# #           loc='best', 
# #           fontsize=16, 
# #           frameon=False)
# # # ax.grid(True)

# # diff_data = H_dimer - J_dimer
# # # # Inset plot
# # ax_inset = plt.axes([0.45, 0.3, 0.45, 0.25])  # [left, bottom, width, height]
# # ax_inset.plot(system_params.gamma_q1 * times, diff_data)
# # ax_inset.set_xlabel(r'$\gamma t$', fontsize=14)
# # ax_inset.set_ylabel(r'$\Delta I(t) / I_0$', fontsize=14)
# # ax_inset.tick_params(axis='x', labelsize=14)
# # ax_inset.tick_params(axis='y', labelsize=14)
# # # ax_inset.set_xlim([0, 4])
# # ax_inset.legend(frameon=False)
# # # ax_inset.set_ylim([-0.04, 0.0])
# # # ax_inset.grid(True)

# # # # Save the plot
# # plt.savefig("orient_intensity_strongC.png", dpi = 500)

# # # Show the plot
# # plt.show()













# # # Load g2 data
# # H_dimer_g2 = np.load('files/g2/H_dimer/H_strongC.npy')
# # J_dimer_g2 = np.load('files/g2/J_dimer/J_strongC.npy')
# # ortho_dimer_g2 = np.load('files/g2/ortho_dimer/ortho_strongC.npy')
# # anti_parallel_g2 = np.load('files/g2/oppo_dimer/oppo_strongC.npy')

# # # Concatenate times and g2 data for symmetry
# # times_sym = np.append(-np.flip(times_ns[1:] - times_ns[0]), times_ns)
# # H_dimer_g2_sym = np.append(np.flip(H_dimer_g2[1:]), H_dimer_g2)
# # J_dimer_g2_sym = np.append(np.flip(J_dimer_g2[1:]), J_dimer_g2)
# # ortho_dimer_g2_sym = np.append(np.flip(ortho_dimer_g2[1:]), ortho_dimer_g2)
# # anti_parallel_g2_sym = np.append(np.flip(anti_parallel_g2[1:]), anti_parallel_g2)

# # # Create main plot
# # fig, ax = plt.subplots(figsize=(8, 6))

# # # Plot the data with appropriate labels
# # line1, = ax.plot(times_sym, anti_parallel_g2_sym, '-k', label='Anti-parallel', marker='s', markersize=6, markevery=25000, linewidth=1)
# # line2, = ax.plot(times_sym, H_dimer_g2_sym, '-b', label='H dimer', linewidth=1)
# # line3, = ax.plot(times_sym, J_dimer_g2_sym, '--r', label='J dimer', linewidth=1)
# # line4, = ax.plot(times_sym, ortho_dimer_g2_sym, '--g', label='Orthogonal', linewidth=1)

# # # Set labels and title
# # ax.set_xlabel(r'$\tau$ (ns)', fontsize=18)
# # ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=18)
# # ax.set_xlim([-10, 10])
# # ax.set_ylim([0.0, 1.1])

# # # Set tick parameters
# # ax.tick_params(axis='x', labelsize=16, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=16, direction='in', length=6)
# # ax.text(-9.5, 1.05, '(b)', size=18)

# # # Add legend and grid
# # ax.legend(handles=[line2, line3, line4, line1], 
# #           labels=['H dimer', 'J dimer', 'Orthogonal', 'Anti-parallel'],
# #           loc='best', 
# #           fontsize=16, 
# #           frameon=False)

# # # Change number of ticks on y-axis
# # ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

# # # Inset plot
# # ax_inset = plt.axes([0.6, 0.3, 0.25, 0.25])  # [left, bottom, width, height]
# # ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), H_dimer_g2_sym, '-b', linewidth=1)
# # ax_inset.plot(np.append(-np.flip(times[1:] - times[0]), times), J_dimer_g2_sym, '-r', linewidth=1)
# # ax_inset.set_xlim([-5, 5])
# # ax_inset.set_ylim([0.0, 1.1])
# # ax_inset.set_xlabel(r'$\tau$ (ps)', fontsize=14)
# # ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=14)
# # ax_inset.tick_params(axis='x', labelsize=14, direction='in', length=6)
# # ax_inset.tick_params(axis='y', labelsize=14, direction='in', length=6)
# # ax_inset.yaxis.set_major_locator(MaxNLocator(nbins=4))


# # # Save the plot
# # # plt.savefig('orient_g2_strongC.png', dpi=500)

# # # Show the plot
# # plt.show()

# # # Load g2 data
# # H_dimer_g2 = np.load('files/g2/H_dimer/H_strongC.npy')
# # J_dimer_g2 = np.load('files/g2/J_dimer/J_strongC.npy')

# # # Concatenate times and g2 data for symmetry
# # times_sym = np.append(-np.flip(times[1:] - times[0]), times)
# # H_dimer_g2_sym = np.append(np.flip(H_dimer_g2[1:]), H_dimer_g2)
# # J_dimer_g2_sym = np.append(np.flip(J_dimer_g2[1:]), J_dimer_g2)

# # # Create new plot focusing on the inset part
# # fig, ax = plt.subplots(figsize=(8, 6))

# # # Plot the data with appropriate labels
# # line1, = ax.plot(times_sym, H_dimer_g2_sym, '-b', label='H dimer')
# # line2, = ax.plot(times_sym, J_dimer_g2_sym, '-r', label='J dimer')

# # # Set labels and title
# # ax.set_xlabel(r'$\tau$ (ps)', fontsize=18)
# # ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=18)
# # ax.set_xlim([-5, 5])
# # ax.set_ylim([-0.1, 1.2])
# # ax.text(-4.9, 1.125, '(b)', size=18)

# # # Set tick parameters
# # ax.tick_params(axis='x', labelsize=16, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=16, direction='in', length=6)

# # # Add legend and grid
# # ax.legend(loc='best', fontsize=16, frameon=False)
# # # ax.grid(True)

# # # Change number of ticks on y-axis
# # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

# # # Save the plot
# # plt.savefig('g2_strongC_inset.png', dpi=1000)

# # # Show the plot
# # plt.show()








# # times_sym = np.append(-np.flip(times_ns[1:] - times_ns[0]), times_ns)
# # g2_perp = np.append(np.flip(np.load('files/g2/ortho_dimer/g2_perp.npy')[1:]), np.load('files/g2/ortho_dimer/g2_perp.npy'))
# # g2_mic_deph = np.append(np.flip(np.load('files/g2/ortho_dimer/g2_mic_deph.npy')[1:]), np.load('files/g2/ortho_dimer/g2_mic_deph.npy'))
# # g2_par = np.append(np.flip(np.load('files/g2/ortho_dimer/g2_par.npy')[1:]), np.load('files/g2/ortho_dimer/g2_par.npy'))

# # # Create figure and axis
# # fig, ax = plt.subplots(figsize=(8, 6))

# # # Plot the data with appropriate labels
# # ax.plot(times_sym, g2_perp, '-b', label='Perpendicular', marker = 's', markersize = 3, markevery = 20000, linewidth = 1)
# # ax.plot(times_sym, g2_mic_deph, 'r', label= r'$\theta = \pi/2, \phi = \pi/4$', linewidth = 1)
# # ax.plot(times_sym, g2_par, '--g', label= 'Parallel', linewidth = 1)

# # # Set labels and title
# # ax.set_xlabel(r'$\tau$ (ns)', fontsize=18)
# # ax.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=18)

# # # Set tick parameters
# # ax.tick_params(axis='x', labelsize=16, direction='in', length=6)
# # ax.tick_params(axis='y', labelsize=16, direction='in', length=6)
# # ax.text(-14, 0.95, '(c)', size=18)

# # # Add legend and grid
# # ax.legend(loc='best', fontsize=16, frameon=False)


# # # Inset plot
# # ax_inset = plt.axes([0.625, 0.3, 0.25, 0.25])  # [left, bottom, width, height]
# # ax_inset.plot(1.0e3 * times_sym, g2_mic_deph, '-r', markersize=0)
# # ax_inset.set_xlabel(r'$t$ (ps)', fontsize=14)
# # ax_inset.set_ylabel(r'$g^{(2)}(\infty, \tau)$', fontsize=14)
# # ax_inset.tick_params(axis='x', labelsize=14)
# # ax_inset.tick_params(axis='y', labelsize=14)
# # ax_inset.set_xlim([-1, 1])
# # # ax_inset.set_ylim([-0.04, 0.0])
# # # ax_inset.grid(True)


# # # Save the plot
# # plt.savefig("orient_g2_ortho.png", dpi=500)

# # # Show the plot
# # plt.show()























