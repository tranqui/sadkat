# + ignore="True"
from droplet import *
# -

# # 4. Running the simulation

# ## 4.1. Define the Graphical User Interface (GUI) for specifying droplet simulations

# +
import ipywidgets as widgets
from IPython.display import display, HTML

from tkinter import Tk
from tkinter import filedialog

class DropletSimulationGUI:
    """This class manages a graphical user interface that allows the user to set
    droplet conditions and initiate simulations.

    Create and run this GUI by executing the following in an input cell:

    >>> gui = DropletSimulationGUI()
    >>> gui.display()
    """

    def __init__(self):
        """Constructor describes the various widgets used in the interface."""

        # The list of droplet trajectories from various simulations run by user.
        self.trajectories = []

        # Formatting options for widgets.
        dropdown_layout = widgets.Layout(width='100%')
        slide_layout = widgets.Layout(width='100%')
        label_layout = widgets.Layout(width='100%', display='flex', justify_content='center')
        centered_layout = widgets.Layout(width='100%', justify_content='center')
        coord_layout = widgets.Layout(width='100%')

        solution_grid_layout = widgets.Layout(grid_template_columns='repeat(2, 30%)')
                                             #justify_content='center')
        coord_grid_layout = widgets.Layout(grid_template_columns='repeat(3, 30%)',
                                           justify_items='center', align_items='center')
        time_grid_layout = widgets.Layout(grid_template_columns='20% 30% 50%',
                                          justify_items='center', align_items='center')
        tab_grid_layout = widgets.Layout(grid_template_columns='repeat(2, 50%)',
                                         #justify_items='center', # uncomment to horizontally align buttons in the center 
                                         align_items='center')

        # Section labels.
        self.droplet_label = widgets.Label('Droplet properties', layout=label_layout)
        self.initial_label = widgets.Label('Initial conditions', layout=label_layout)
        self.ambient_label = widgets.Label('Ambient conditions', layout=label_layout)
        self.simulation_label = widgets.Label('Simulation parameters', layout=label_layout)
        self.labels = [self.droplet_label, self.initial_label, self.ambient_label, self.simulation_label]
        # Hack to change the formatting of the section labels so they are more prominent.
        display(HTML("<style>.bolded { font-weight: bold; font-size: large }</style>"))
        for label in self.labels: label.add_class("bolded")

        ## Droplet properties.

        self.solution_dropdown = widgets.Dropdown(
            options=all_solutions.keys(),
            value=list(all_solutions.keys())[0],
            description='Solution:',
            disabled=False,
            layout=dropdown_layout
        )

        self.density_dropdown = widgets.Dropdown(
            options=['regular', 'volume additivity'],
            value='regular',
            description='Density fit:',
            disabled=False,
            layout=dropdown_layout
        )
        self.density_dropdown.observe(self.changed_density_fit, 'value')

        self.profile_dropdown = widgets.Dropdown(
            options=['uniform', 'radial'],
            value='uniform',
            description='Profile:',
            disabled=True,
            layout=dropdown_layout
        )
        self.profile_dropdown.observe(self.changed_profile, ['value', 'disabled'])

        self.npoints_select = widgets.BoundedIntText(
            value=100,
            min=10,
            max=10000,
            step=10,
            description='npoints:',
            disabled=True,
            layout=dropdown_layout
        )

        #self.solution_choices = widgets.HBox([self.solution_dropdown, self.density_dropdown])
        self.solution_choices = widgets.GridBox([self.solution_dropdown,
                                                 self.density_dropdown,
                                                 self.profile_dropdown,
                                                 self.npoints_select],
                                                 layout=solution_grid_layout)

        ## Initial conditions.

        self.initial_mfs_slider = widgets.FloatSlider(
            description='MFS', layout=slide_layout,
            value=0.05, min=0, max=1, step=0.01, readout_format='.2f')

        self.initial_R_slider = widgets.FloatSlider(
            description='R / µm', layout=slide_layout,
            value=15, min=0.1, max=500, step=0.1, readout_format='.1f')

        self.initial_T_slider = widgets.FloatSlider(
            description='T / ℃', layout=slide_layout,
            value=(body_temperature - T_freezing), min=0, max=100, step=0.1, readout_format='.1f')

        self.initial_x = widgets.FloatText(value=0, step=0.1, description='x / m', layout=coord_layout)
        self.initial_y = widgets.FloatText(value=0, step=0.1, description='y / m', layout=coord_layout)
        self.initial_z = widgets.FloatText(value=0, step=0.1, description='z / m', layout=coord_layout)
        self.initial_vx = widgets.FloatText(value=0, step=0.1, description='vx / m/s', layout=coord_layout)
        self.initial_vy = widgets.FloatText(value=0, step=0.1, description='vy / m/s', layout=coord_layout)
        self.initial_vz = widgets.FloatText(value=0, step=0.1, description='vz / m/s', layout=coord_layout)
        self.initial_coordinates = widgets.GridBox([self.initial_x, self.initial_y, self.initial_z,
                                                    self.initial_vx, self.initial_vy, self.initial_vz],
                                                    layout=coord_grid_layout)
        ## Ambient conditions.

        self.ambient_T_slider = widgets.FloatSlider(
            description='T / ℃', layout=slide_layout,
            value=(standard_temperature - T_freezing), min=0, max=100, step=0.1, readout_format='.1f')

        self.ambient_RH_slider = widgets.FloatSlider(
            description='RH', layout=slide_layout,
            value=0.5, min=0, max=1, step=0.01, readout_format='.2f')

        self.gas_vx = widgets.FloatText(value=0, step=0.1, description='vx / m/s', layout=coord_layout)
        self.gas_vy = widgets.FloatText(value=0, step=0.1, description='vy / m/s', layout=coord_layout)
        self.gas_vz = widgets.FloatText(value=0, step=0.1, description='vz / m/s', layout=coord_layout)
        self.gas_velocity = widgets.GridBox([self.gas_vx, self.gas_vy, self.gas_vz],
                                             layout=coord_grid_layout)

        self.gx = widgets.FloatText(value=0, step=0.1, description='gx / m/s', layout=coord_layout)
        self.gy = widgets.FloatText(value=0, step=0.1, description='gy / m/s', layout=coord_layout)
        self.gz = widgets.FloatText(value=-gravitational_acceleration,
                                    step=0.1, description='gz / m/s', layout=coord_layout)
        self.gravity = widgets.GridBox([self.gx, self.gy, self.gz],
                                       layout=coord_grid_layout)

        ## Simulation parameters.

        self.time_selection = widgets.BoundedFloatText(
            value=10, min=0.1, max=1e6, step=0.1,
            description='time / s',
            layout=slide_layout
        )

        self.stop_checkbox = widgets.Checkbox(
            value=True,
            description='terminate on equilibration',
            indent=True,
            layout=slide_layout
        )
        self.stop_checkbox.observe(self.clicked_stop_checkbox, 'value')

        self.stop_threshold_slider = widgets.FloatLogSlider(
            description='threshold', layout=slide_layout,
            value=1e-3, base=10, min=-8, max=-1, step=1, readout_format='.1g'
        )

        self.time_choices = widgets.GridBox([self.time_selection, self.stop_checkbox, self.stop_threshold_slider],
                                             layout=time_grid_layout)

        self.timestep_slider = widgets.FloatLogSlider(
            description='timestep / s', layout=slide_layout,
            value=1e-2, base=10, min=-5, max=-1, step=0.1, readout_format='.1g'
        )

        ## The start simulation button.

        self.run_button = widgets.Button(
            value=False,
            description='Simulate droplet',
            button_style='info',
            tooltip='Run simulation',
            icon='tint',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='50ex')
        )

        self.run_button.on_click(self.run_simulation)

        self.centered_run_button = widgets.HBox([self.run_button], layout=centered_layout)

        ## Space to show warnings etc.

        self.status_space = widgets.HBox([], layout=centered_layout)

        ## Create tabs for showing plots of the various outputs.

        tab_labels = ['radius', 'mass', 'temperature', 'position']
        self.tabs = {label: widgets.GridBox([], layout=tab_grid_layout) for label in tab_labels}

        self.output_tabs = widgets.Tab(children=list(self.tabs.values()))
        for i, label in enumerate(tab_labels):
            self.output_tabs.set_title(i, label)

        ## Putting the entire interface together.

        self.display_widgets = [
                self.droplet_label,
                self.solution_choices,
                self.initial_label,
                self.initial_mfs_slider,
                self.initial_R_slider,
                self.initial_T_slider,
                self.initial_coordinates,
                self.ambient_label,
                self.ambient_T_slider,
                self.ambient_RH_slider,
                self.gas_velocity,
                self.gravity,
                self.simulation_label,
                self.time_choices,
                self.timestep_slider,
                self.centered_run_button,
                self.status_space,
                self.output_tabs
        ]

    def changed_density_fit(self, dropdown):
        """Callback when the 'density fit' dropdown is changed.

        If volume additivity is selected, then we can enable the options to allow the concentration
        profile to vary.

        Args:
            dropdown: the dropdown widget that is changed. This is not used, but it is needed
                      to assign this function as a callback to clicking the dropdown.
        """
        self.profile_dropdown.disabled = not self.density_dropdown.value == 'volume additivity'

        if not self.density_dropdown.value == 'volume additivity':
            self.profile_dropdown.value = 'uniform'

    def changed_profile(self, dropdown):
        """Callback when the 'profile' dropdown is changed.

        Unless we are doing inhomogeneous concentration profiles, we do not need to let the user
        specify the number of points so this function disables that widget.

        Args:
            dropdown: the dropdown widget that is changed. This is not used, but it is needed
                      to assign this function as a callback to clicking the dropdown.
        """
        self.npoints_select.disabled = self.profile_dropdown.value == 'uniform' or \
                                       not self.density_dropdown.value == 'volume additivity'

    def clicked_stop_checkbox(self, checkbox):
        """Callback when the 'terminate on equilibration' checkbox is clicked.

        We disable the termination threshold slider when this box is not clicked, because that parameter
        has no purpose outside of that termination condition.

        Args:
            checkbox: the checkbox widget that is clicked. This is not used, but it is needed
                      to assign this function as a callback to clicking the checkbox.
        """
        self.stop_threshold_slider.disabled = not self.stop_checkbox.value

    @property
    def gui_state(self):
        """Get the state of *all* GUI settings describing the parameters of a simulation.

        It is useful to keep these with trajectories for reference (so e.g. old simulations can be repeated).

        Return:
            dictionary containing all parameters defining a simulation.
        """
        return dict(
            # Droplet properties
            solution = self.solution_dropdown.value,
            density_fit = self.density_dropdown.value,
            profile = self.profile_dropdown.value,
            npoints = self.npoints_select.value,

            # Initial conditions of droplet.
            initial_radius = self.initial_R_slider.value,
            initial_mfs = self.initial_mfs_slider.value,
            initial_temperature = self.initial_T_slider.value,
            initial_velocity = np.array((self.initial_vx.value, self.initial_vy.value, self.initial_vz.value)),
            initial_position = np.array((self.initial_x.value, self.initial_y.value, self.initial_z.value)),

            # Environment (i.e. the gas phase).
            ambient_temperature = self.ambient_T_slider.value,
            ambient_RH = self.ambient_RH_slider.value,
            gas_velocity = np.array((self.gas_vx.value, self.gas_vy.value, self.gas_vz.value)),
            gravity = np.array((self.gx.value, self.gy.value, self.gz.value)),

            # Simulation parameters.
            time = self.time_selection.value,
            timestep = self.timestep_slider.value,
            terminate_on_equilibration = self.stop_checkbox.value,
            terminate_threshold = self.stop_threshold_slider.value
        )

    def run(self, solution, density_fit, profile, npoints,
            initial_radius, initial_mfs, initial_temperature, initial_velocity, initial_position,
            ambient_temperature, ambient_RH, gas_velocity, gravity,
            time, timestep, terminate_on_equilibration, terminate_threshold):
        """
        Run a new simulation.

        Args:
            solution: string specifying the solution i.e. (composition of droplet).
            density_fit: string specifying the type of density fit to use (valid values are 'regular' or
                         'volume additivity').
            profile: string either 'uniform' for uniform droplets or 'radial' for droplets with a
                     radially varying concentration profile. For the latter, the 'volume additivity'
                     density fit must be selected.
            npoints: number of points to use for radially varying concentration profiles (only used if
                     profile='radial').
            initial_radius: initial droplet radius in microns.
            initial_mfs: initial mass fraction of solute in droplet, should be between 0 and 1.
            initial_temperature: initial temperature of droplet in Celsius.
            initial_velocity: a 3-dimensional vector specifying the initial droplet velocity in metres/second.
            initial_position: a 3-dimensional vector specifying the initial droplet position in metres.
            ambient_temperature: temperature of surrounding gas in Celsius.
            ambient_RH: relative humidity of surrounding gas, should be between 0 and 1.
            gas_velocity: a 3-dimensional vector specifying the velocity of background gas in metres/second.
            gravity: a 3-dimensional acceleration vector from body forces (e.g. gravity).
            time: the total time to run the simulation for in seconds.
            timestep: the timestep in the simulation in seconds. A smaller number will make the simulation
                      more accurate, but also take longer. NB: this is technically the *maximum* timestep,
                      as the integration algorithm will sometimes decide to take smaller timesteps.
            terminate_on_equilibration: a boolean; if True then the simulation will terminate early if
                      the droplet has stopped evolving (as determined by a threshold in the evaporation rate)
            terminate_threshold: the threshold for the previous termination condition. If the evaporation rate
                      divided by the initial mass falls below this threshold, then the simulation terminates.
        Returns:
            The final droplet state.
            The trajectory of (independent) droplet variables from which its complete history can be
                reconstructed.
        """

        ## Unit conversions.

        initial_radius *= 1e-6             # convert to metres
        initial_temperature += T_freezing  # convert to Kelvin
        ambient_temperature += T_freezing  # convert to Kelvin

        ## Set up droplet and environment.

        if density_fit == 'regular':
            solution = all_solutions[solution]
        elif density_fit == 'volume additivity':
            solution = all_solutions_volume_additivity[solution]
        else:
            raise ValueError('unrecognised density fit: %s' % density_fit)

        gas = Atmosphere(ambient_temperature, ambient_RH, velocity=gas_velocity)

        if profile == 'uniform':
            droplet = UniformDroplet.from_mfs(solution,
                                              gas,
                                              gravity,
                                              initial_radius,
                                              initial_mfs,
                                              initial_temperature,
                                              initial_velocity,
                                              initial_position)
        elif profile == 'radial':
            raise NotImplementedError('varying droplets not yet implemented!')
        else:
            raise ValueError('unrecognised concentration profile: %s' % profile)

        trajectory = droplet.integrate(time,
                                       timestep,
                                       terminate_on_equilibration,
                                       terminate_threshold)

        return droplet, trajectory

    def run_simulation(self, button):
        """Run a simulation with the currently selected settings.

        This function is called when the simulate button is clicked by the user.

        Args:
            button: the button widget that is clicked. This is not used, but it is needed
                    to assign this function as a callback to clicking the button.
        """

        ## Run the simulation with the current settings.
        settings = self.gui_state
        label_layout = widgets.Layout(width='100%', display='flex', justify_content='center')
        try:
            self.status_space.children = (widgets.Label("Running simulation...", layout=label_layout),)
            droplet, trajectory = self.run(**settings)
        except Exception as e:
            self.status_space.children = (widgets.Button(description=str(e),
                                                         disabled=True,
                                                         button_style='warning',
                                                         icon='exclamation-triangle',
                                                         style={'font_weight': 'bold'},
                                                         layout=label_layout),)
            return

        # On successful run, clear the status space.
        self.status_space.children = tuple()

        ## Save and output this new trajectory with some control widgets.

        # Add a new plot for this trajectory in each tab.
        plots = []
        for tab in self.tabs.values():
            plots += [widgets.Output()]
            tab.children = tab.children + (plots[-1],)
        self.plot_trajectory(droplet, trajectory, *plots)

        # Add a button so that the user can save this data to a file.
        save_button = widgets.Button(
            description='Save trajectory',
            button_style='success',
            tooltip='Save data to file',
            icon='download',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='25ex')
        )
        save_button.on_click(self.save_trajectory)

        # Add a button so that the user can get more information about this trajectory.
        info_button = widgets.Button(
            description='Info',
            button_style='info',
            tooltip='More information about this trajectory (e.g. the initial conditions)',
            icon='info',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='25ex')
        )
        info_button.on_click(self.info_trajectory)

        # Add a button so that the user can clear this data.
        remove_button = widgets.Button(
            description='Remove',
            button_style='danger',
            tooltip='Remove this trajectory',
            icon='times',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='25ex')
        )
        remove_button.on_click(self.remove_trajectory)

        # Create a space for the buttons to the right of the plot.
        buttons = widgets.HBox([widgets.VBox([save_button, info_button, remove_button])],
                                layout=widgets.Layout(align_items='center'))
        for tab in self.tabs.values():
            tab.children = tab.children + (buttons,)

        # Save the data in case the user wants to dump it to a file (or otherwise use it).
        self.trajectories += [dict(droplet=droplet, trajectory=trajectory, settings=settings,
                                   plots=plots, buttons=buttons)]

    def save_trajectory(self, button):
        """Save a simulation to a file.

        This function is called when the save button is clicked by the user.

        Args:
            button: the save button widget that is clicked. We use its identity to determine which
                    trajectory to save.
        """

        # Find the trajectory in our saved list.
        for index,trajectory in enumerate(self.trajectories):
            if button in trajectory['buttons'].children[0].children:
                break
        else: # if the trajectory was not found, we have a problem!
            raise LookupError('could not find trajectory!')

        ## Prompt the user for a filepath to save the simulation, using Tkinter

        # We do not want to create a new window, so keep the Tkinter root window from appearing
        Tk().withdraw()
        # Prompt user for save path.
        filename = filedialog.asksaveasfilename(filetypes = (("comma separated values files", "*.csv"),
                                                             ("all files", "*.*")))

        ## Save the trajectory.

        if filename:
            droplet = trajectory['droplet']
            history = trajectory['trajectory']
            settings = trajectory['settings']
            settings = '\n'.join(['# %s: %s' % (key, value) for key,value in settings.items()])
            data = droplet.complete_trajectory(history)
            with open(filename, 'w') as f:
                f.write(settings)
                f.write('\n')
                data.to_csv(f, index=False)

    def info_trajectory(self, button):
        """Show more information about a previously run simulation.

        This function is called when the info button is clicked by the user.

        Args:
            button: the info button widget that is clicked. We use its identity to determine which
                    trajectory to show more information on.
        """

        # Find the trajectory in our saved list.
        for index,trajectory in enumerate(self.trajectories):
            if button in trajectory['buttons'].children[0].children:
                break
        else: # if the trajectory was not found, we have a problem!
            raise LookupError('could not find trajectory!')

        buttons = trajectory['buttons']
        settings = trajectory['settings']

        # If we're already showing the information, remove it.
        if len(buttons.children) > 1:
            buttons.children = tuple(list(buttons.children)[:-1])
        # Otherwise, display the information:
        else:
            css_style = '<style>p{word-wrap: break-word; line-height: 1.2em; margin-left: 1em}</style>'
            formatted_settings = '\n'.join(['<p>%s: %s</p>' % (key, value) for key,value in settings.items()])
            text = widgets.HTML(css_style + formatted_settings)
            buttons.children = buttons.children + (text,)

    def remove_trajectory(self, button):
        """Remove a simulation from the saved list.

        This function is called when the close button is clicked by the user.

        Args:
            button: the close button widget that is clicked. We use its identity to determine which
                    trajectory to remove.
        """

        # Find and remove the trajectory in our saved list.
        for index,trajectory in enumerate(self.trajectories):
            if button in trajectory['buttons'].children[0].children:
                break
        else: # if the trajectory was not found, we have a problem!
            raise LookupError('could not find trajectory!')
        self.trajectories.pop(index)

        # Remove the output rows/widgets from the display.
        for tab in self.tabs.values():
            # Remove any plots.
            tab.children = tuple(child for child in tab.children if child not in trajectory['plots'])
            # Remove the buttons.
            tab.children = tuple(child for child in tab.children if child is not trajectory['buttons'])

    def plot_trajectory(self, droplet, trajectory, out_radius, out_mass, out_temperature, out_position):
        """Plot summary data for a new trajectory.
        
        Args:
            droplet: the Droplet object, containing important information about the trajectory.
            trajectory: the trajectory of independent variables itself.
            out_radius: output pane for the radius trajectory plot.
            out_mass: output pane for the mass trajectory plot.
            out_temperature: output pane for the temperature trajectory plot.
            out_position: output pane for the position trajectory plot.
        """

        # Get the equilibrium state so we can compare how far/close to it we are.
        try:
            # If RH = 1, there is no equilibrium state, so we need to catch this.
            with np.errstate(divide='raise'):
                equilibrium = droplet.equilibrium_droplet
        except:
            # This should only happen if RH = 1.
            assert droplet.environment.relative_humidity == 1
            equilibrium = droplet.copy()
            equilibrium.mass_solvent = np.inf

        data = droplet.complete_trajectory(trajectory)

        t = data['time']
        radius = data['radius']
        mass_solvent = data['mass_solvent']
        temperature = data['temperature']
        x = data['x']
        z = data['z']

        figsize = (6, 4)

        with out_radius:
            plt.figure(figsize=figsize)
            plt.plot(t, radius, label='trajectory')
            plt.axhline(y=equilibrium.radius, ls='dashed', label='equilibrium')
            plt.legend(loc='best')
            plt.xlabel('t (s)')
            plt.ylabel('R (m)')
            plt.show()

        with out_mass:
            plt.figure(figsize=figsize)
            plt.plot(t, droplet.mass_solute + mass_solvent, label='total')
            plt.plot(t, mass_solvent, label='solvent')
            plt.axhline(y=equilibrium.mass, ls='dashed', label='equilibrium')
            plt.axhline(y=droplet.mass_solute, ls='dotted', label='solute')
            plt.legend(loc='best')
            plt.xlabel('t (s)')
            plt.ylabel('m (kg)')
            plt.show()

        with out_temperature:
            plt.figure(figsize=figsize)
            plt.plot(t, temperature - T_freezing, label='trajectory')
            plt.axhline(y=(equilibrium.temperature - T_freezing), ls='dashed', label='equilibrium')
            plt.legend(loc='best')
            plt.xlabel('t (s)')
            plt.ylabel('T (℃)')
            plt.show()

        with out_position:
            plt.figure(figsize=figsize)
            plt.plot(x, z)
            plt.xlabel('x (m)')
            plt.ylabel('z (m)')
            plt.show()

    def display(self):
        """Show the GUI display."""
        display(*self.display_widgets)

    def clear(self):
        """Clear all previous simulations."""
        self.trajectories = []
        for tab in self.tabs.values():
            tab.children = []
# -

# ## 4.2. Executing the graphical program to run simulations

# +
if __name__ == '__main__':
    gui = DropletSimulationGUI()
    gui.display()
# -

# # 4.3 Iterating over input parameters

# +


def simulate(time, timestep,
             solution, ambient_temperature, ambient_RH,
             initial_radius, initial_temperature, initial_mfs,
             initial_velocity = np.zeros(3), # metres/second
             initial_position = np.zeros(3), # metres/second
             gravity = np.zeros(3)):

    gas = Atmosphere(ambient_temperature, ambient_RH)
    droplet = UniformDroplet.from_mfs(solution,
                                      gas,
                                      gravity,
                                      initial_radius,
                                      initial_mfs,
                                      initial_temperature,
                                      initial_velocity,
                                      initial_position)
    trajectory = droplet.integrate(time, timestep, terminate_on_equilibration=True, eps = 0.001)

    return droplet, trajectory


# -

# +

R0 = 25e-6 # metres
T = 293.15 # Kelvin
mfs = 0

time = 25 # seconds
timestep = 0.1 # seconds

history_list = []

RH_range = np.sqrt(np.linspace(0,100**2, 100)) / 100 #np.arange(0,1.001,0.1)

for RH in RH_range:
    
    droplet, trajectory = simulate(time, timestep, solution, T, RH, R0, T, mfs)
    
    # Obtain a table giving a history of *all* droplet parameters.
    history = droplet.complete_trajectory(trajectory)
    history_list.append(history)

RH = 0.9
specific_droplet, specific_trajectory = simulate(time, timestep, solution, T, RH, R0, T, mfs)


# -

# +


import matplotlib as mpl
cmap = mpl.cm.cool_r
norm = mpl.colors.Normalize(vmin=100 * RH_range.min(), vmax=100 * RH_range.max())

colors = plt.cm.cool_r(RH_range)

for history, RH, color in zip(history_list, RH_range, colors):
    
    
    plt.plot(history['time'], history['radius'] / 1e-6, c = color)
    

RH = 0.9
trajectory = specific_droplet.complete_trajectory(specific_trajectory)
plt.plot(trajectory['time'], trajectory['radius'] / 1e-6, '--', label = 100 * RH)
    
plt.xlabel('Time / s')
plt.ylabel('Radius / µm')
plt.hlines((0,25),0,25,'k')


plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), label = '% RH' )
plt.legend()
plt.show()

# -