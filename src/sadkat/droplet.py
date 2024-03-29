# + ignore="True"
from sadkat.gas import *
# -

# # 3. Defining the droplet and its evolution
#
# NB: See the appendix at the bottom of this notebook for the complete set of equations used to describe the droplet.

# ## 3.1. Definition

@dataclass
class UniformDroplet:
    """This class completely describes the state of the droplet during its evolution.

    The solute is assumed to be distributed uniformly throughout the droplet.
    """

    solution: object
    environment: object
    gravity: object                 # m/s^2
    mass_solute: float              # kg
    mass_solvent: float             # kg
    temperature: float              # K
    velocity: np.array=np.zeros(3)  # m/s
    position: np.array=np.zeros(3)  # m

    @staticmethod
    def from_mfs(solution, environment, gravity,
                 radius, mass_fraction_solute, temperature,
                 velocity=np.zeros(3), position=np.zeros(3)):
        """Create a droplet from experimental conditions.

        Args:
            solution: parameters describing solvent+solute
            environment: parameters of gas surrounding droplet
            body acceleration in metres/second^2 (3-dimensional vector)
            radius in metres
            mass_fraction_solute (MFS) (unitless)
            temperature in K
            velocity in metres/second (3-dimensional vector)
            position in metres (3-dimensional vector)
        """
        mass = 4*np.pi/3 * radius**3 * solution.density(mass_fraction_solute)
        mass_solvent = (1-mass_fraction_solute) * mass
        mass_solute = mass_fraction_solute * mass
        return UniformDroplet(solution, environment, gravity, mass_solute, mass_solvent, temperature, velocity, position)

    @property
    def state(self):
        """Variables that determine droplet state during its evolution, i.e. the independent variables."""
        return np.hstack((self.mass_solvent, self.temperature, self.velocity, self.position))

    @state.setter
    def state(self, new_state):
        """Update the state of the droplet by changing its independent variables.

        This is useful for reconstructing droplets at e.g. earlier points in time, or in its final
        equilibrium state.

        Args:
            new_state: the new set of independent variables.
        """
        try: self.mass_solvent, self.temperature, self.velocity, self.position = new_state
        except:
            x = new_state
            self.mass_solvent, self.temperature, self.velocity, self.position = x[0], x[1], x[2:5], x[5:]

    @property
    def complete_state(self):
        """All droplet variables, including both independent and dependent variables that completely
        determine all droplet properties.

        This form is ready for a row within a table specifying e.g. a droplet's trajectory.
        """
        return dict(mass = self.mass,
                    mass_solute = self.mass_solute,
                    mass_solvent = self.mass_solvent,
                    mass_fraction_solute = self.mass_fraction_solute,
                    mass_fraction_solvent = self.mass_fraction_solvent,
                    mole_fraction_solute = self.mole_fraction_solute,
                    mole_fraction_solvent = self.mole_fraction_solvent,
                    density = self.density,
                    radius = self.radius,
                    # This is too similar to radius so it can probably be omitted.
                    #diameter = self.diameter,
                    vapour_pressure = self.vapour_pressure,
                    temperature = self.temperature,
                    drag_coefficient = self.drag_coefficient,
                    reynolds_number = self.reynolds_number,
                    schmidt_number = self.schmidt_number,
                    prandtl_number = self.prandtl_number,
                    sherwood_number = self.sherwood_number,
                    nusselt_number = self.nusselt_number,
                    vx=self.velocity[0],
                    vy=self.velocity[1],
                    vz=self.velocity[2],
                    speed=self.speed,
                    x=self.position[0],
                    y=self.position[1],
                    z=self.position[2],
                    gx=self.gravity[0],
                    gy=self.gravity[1],
                    gz=self.gravity[2],
                    # I've assumed the user doesn't care about the derivatives, because these can be roughly inferred
                    # from the changes between frames. If more precise information is needed, then uncomment these.
                    #evaporation_rate = self.dmdt,
                    #dTdt = self.dTdt,
                    #ax = self.dvdt[0],
                    #ay = self.dvdt[1],
                    #az = self.dvdt[2]
                   )

    @property
    def mass(self):
        """Total mass of droplet (both solvent and solute components) in kg."""
        return self.mass_solute + self.mass_solvent

    @property
    def mass_fraction_solute(self):
        """Mass fraction of solute (i.e. the non-volatile component).
        NB: Should be zero for pure solvent."""
        return self.mass_solute / self.mass

    @property
    def mass_fraction_solvent(self):
        """Mass fraction of solvent."""
        return 1 - self.mass_fraction_solute

    @property
    def mole_fraction_solute(self):
        """Mole fraction of solute (i.e. the non-volatile component).
        NB: Should be zero for pure solvent."""
        return self.solution.mole_fraction_solute(self.mass_fraction_solute)

    @property
    def mole_fraction_solvent(self):
        """Mole fraction of solvent."""
        return 1 - self.mole_fraction_solute

    @property
    def density(self):
        """Droplet density in kg/m^3."""
        return self.solution.density(self.mass_fraction_solute)

    @property
    def radius(self):
        """Droplet radius in metres."""
        return (self.mass / (4*np.pi/3 * self.density))**(1/3)

    @property
    def diameter(self):
        """Droplet diameter in metres."""
        return 2*self.radius

    @property
    def vapour_pressure(self):
        """Vapour pressure at gas-liquid boundary in Pascals."""
        P = self.solution.solvent_activity(self.mass_fraction_solute) * self.solution.solvent.equilibrium_vapour_pressure(self.temperature)
        P *= kelvin_effect(self.solution.solvent.surface_tension(self.temperature),
                           self.solution.solvent.density(self.temperature),
                           self.solution.solvent.molar_mass,
                           self.temperature,
                           self.radius)
        return P

    @property
    def surface_solvent_activity(self):
        """Solvent activity at surface."""
        return self.solution.solvent_activity(self.mass_fraction_solute)

    @property
    def speed(self):
        """Magnitude of velocity vector in metres/second."""
        return np.linalg.norm(self.velocity)

    # @property
    # def jet_velocity(self):

    #     jet_initial_velocity = np.array([1,0,0]) * self.jet_initial_speed


    #     jet_dispersion_distance = 6.8 #assuming 6.8 from Xie paper
    #     jet_centreline_speed = (jet_initial_speed * jet_dispersion_distance) / (self.position[0] / self.aperture_diameter)

    #     jet_radial_velocity = 1
    #     jet_axial_velocity = 1

    #     theta = np.arctan2(self.position[2], self.position[0])
    #     r = np.linalg.norm(self.position)
    #     jet_velocity = np.array([1,
    #                              0,
    #                              3])


    #     jet_centreline_temperature = self.environment.temperature + (self.jet_initial_temperature - self.environment.temperature) *

    #                                 (5 / s_bar) * np.sqrt(self.jet_initial_temperature / self.environment.temperature)


    #     #This assumes that the closest point on the centreline wil be vertically above or below the droplet
    #     jet_centreline_position = np.array([self.position[0], 0
    #                                         np.sqrt(self.aperture_area ) * 0.0354 * self.jet_archimedes_number *
    #                                         (self.position[0] / np.sqrt(self.aperture_area) ) ** 3 *
    #                                         np.sqrt(self.jet_initial_temperature / self.environment.temperature) ])
    #     return 0

    @property
    def relative_velocity(self):
        """Velocity relative to environment in metres/second."""
        return self.velocity - self.environment.velocity
        #return self.velocity - self.jet.velocity

    @property
    def relative_speed(self):
        """Magnitude of relative velocity vector in metres/second."""
        return np.linalg.norm(self.relative_velocity)

    @property
    def drag_coefficient(self):
        """Non-dimensional number describing strength of drag forces."""
        Re = self.reynolds_number
        if Re > 1000: return 0.424
        elif Re < 1e-12: return np.inf
        else: return (24 / Re) * (1 + Re**(2/3) / 6)

    @property
    def aperture_diameter(self):
        """"meters"""
        return 0.02

    @property
    def aperture_area(self):
        """"m^2"""
        return self.aperture_diameter ** 2

    @property
    def jet_initial_temperature(self):
        return body_temperature

    @property
    def jet_initial_speed(self):
        return 1

    # @property
    # def jet_archimedes_number(self):
    #     """""""
    #     return np.linalg.norm(self.gravity) * np.sqrt(aperture_area) * volumetric_expansion_coeffcient *
    #                                 (self.jet_initial_temperature - self.environment.temperature) * / self.jet_initial_speed ** 2

    @property
    def reynolds_number(self):
        """Non-dimensional number describing the type of fluid flow."""
        return self.environment.density * self.diameter * self.speed / self.environment.dynamic_viscosity

    @property
    def schmidt_number(self):
        """Non-dimensional number describing the ratio of momentum diffusivity to mass diffusivity."""
        D_function = self.solution.solvent.vapour_binary_diffusion_coefficient
        T_inf = self.environment.temperature
        D_inf = D_function(T_inf)
        return self.environment.dynamic_viscosity / (self.environment.density * D_inf)

    @property
    def prandtl_number(self):
        """Non-dimensional number describing the ratio of momentum diffusivity to thermal diffusivity."""
        return self.environment.specific_heat_capacity * self.environment.dynamic_viscosity / self.environment.thermal_conductivity

    @property
    def sherwood_number(self):
        """Non-dimensional number describing mass transfer."""
        Re = self.reynolds_number
        Sc = self.schmidt_number
        return 1 + 0.3 * Re**(1/2) * Sc**(1/3)

    @property
    def nusselt_number(self):
        """Non-dimensional number describing conductive heat transfer."""
        Re = self.reynolds_number
        Pr = self.prandtl_number
        return 1 + 0.3 * Re**(1/2) * Pr**(1/3)

    @property
    def knudsen_number(self):
        return self.environment.mean_free_path / self.radius

    @property
    def fuchs_sutugin_correction(self):
        Kn = self.knudsen_number
        alpha = 1 # parameter in Fuchs-Sutugin theory that we can take to be one (for now).
        return (1 + Kn) / (1 + (4/3*(1 + Kn)/alpha + 0.377)*Kn)

    @property
    def dmdt(self):
        """Time derivative of mass, i.e. the rate of evaporation in kg/s."""
        Sh = self.sherwood_number

        D_function = self.solution.solvent.vapour_binary_diffusion_coefficient
        lam = D_function.lam

        T_inf = self.environment.temperature
        T = self.temperature
        D_inf = D_function(T_inf)

        # Apply temperature correction to diffusion coefficient appearing in mass flux.
        eps = 1e-8
        if np.abs(T_inf - T) < eps: C = 1 # ensure numerical stability as T -> T_inf
        else: C = (T_inf - T) / T_inf**(lam-1) * (2 - lam) / (T_inf**(2-lam) - T**(2-lam))
        D_eff = C * D_inf

        I = np.log((self.environment.pressure - self.vapour_pressure) /
                   (self.environment.pressure - self.environment.vapour_pressure))

        beta = self.fuchs_sutugin_correction

        return 4*np.pi*self.radius * self.environment.density * (self.solution.solvent.molar_mass / self.environment.molar_mass) * D_eff * Sh * beta * I

    @property
    def dTdt(self):
        """Time derivative of temperature from heat flux at the surface in K/s."""
        Nu = self.nusselt_number
        Gamma = stefan_boltzmann_constant

        r = self.radius
        m = self.mass
        rho = self.density

        T = self.temperature
        T_inf = self.environment.temperature
        K = self.environment.thermal_conductivity

        L = self.solution.solvent.specific_latent_heat_vaporisation(T)
        c = self.solution.solvent.specific_heat_capacity(T)
        r = self.radius

        return 3*K * (T_inf - T) * Nu / (c*rho*r**2) + L*self.dmdt / (c*m) - 3*Gamma * (T**4 - T_inf**4) / (c*rho*r)

    @property
    def cunningham_slip_correction(self):
        """Correction to Stokes' law for drag for small particles due to the onset
        of slip on the particle surface."""

        # Phenomenological parameters in the theory due to Davies (1945):
        A1 = 1.257
        A2 = 0.400
        A3 = 1.100

        Kn = self.knudsen_number
        return 1 + Kn * (A1 + A2*np.exp(-A3/Kn))

    @property
    def dvdt(self):
        """Time derivative of velocity, i.e. its acceleration from Newton's second law in m/s^2."""
        rho_p = self.density
        rho_g = self.environment.density
        g = self.gravity

        buoyancy = 1 - rho_g/rho_p
        acceleration = buoyancy*g

        C = self.drag_coefficient / self.cunningham_slip_correction
        if np.isfinite(C):
            acceleration -= 3*C*rho_g * self.relative_speed * self.relative_velocity / (8*rho_p*self.radius)

        return acceleration

    @property
    def drdt(self):
        """Time derivative of droplet position, i.e. its velocity in m/s."""
        return self.velocity

    @property
    def dxdt(self):
        """Time derivative of all state variables as a vector."""
        return np.hstack((self.dmdt, self.dTdt, self.dvdt, self.drdt))

    @property
    def equilibrium_state(self):
        """Final state of the droplet once it has reached equilibrium."""

        temperature = self.environment.temperature
        velocity = np.zeros(3)
        position = np.zeros(3)

        # Find the right amount of solvent to equalise the vapour pressure across the gas-liquid boundary.
        vapour_pressure = self.environment.vapour_pressure
        solvent_activity = vapour_pressure / self.solution.solvent.equilibrium_vapour_pressure(temperature)
        mass_fraction_solute = self.solution.mass_fraction_solute_from_solvent_activity(solvent_activity)
        mass = self.mass_solute / mass_fraction_solute
        mass_solvent = mass - self.mass_solute

        return mass_solvent, temperature, velocity, position

    @property
    def equilibrium_droplet(self):
        """Final droplet once it has reached equilibrium."""
        return UniformDroplet(self.solution, self.environment, self.gravity,
                              self.mass_solute, *self.equilibrium_state)

    def virtual_droplet(self, x):
        """Create droplet with new state variables.

        This is useful for creating droplets at future times.

        Args:
            x: new state variables (cf. Droplet.state function)
        """
        x = (x[0], x[1], x[2:5], x[5:])
        return UniformDroplet(self.solution, self.environment, self.gravity, self.mass_solute, *x)

    def copy(self):
        """Create an identical copy of this droplet."""
        return self.virtual_droplet(self.state.copy())

    def integrate(self, t, rtol=1e-8,
                  terminate_on_equilibration=False, equ_threshold=1e-4,
                  terminate_on_efflorescence=False, eff_threshold=0.5,
                  first_step=1e-12):
        """Integrate the droplet state forward in time.

        This solves an initial value problem with the current state as the initial conditions.
        The droplet state is updated to the final state after the integration.

        Args:
            t: total time to integrate over (s).
            rtol: relative tolerance used to set dynamic integration timestep between frames in
                trajectory (s). Smaller tolerance means a more accurate integration.
                NB: if numerical artifacts occur in the resulting trajectory, that suggests this
                parameter needs to be decreased.
            terminate_on_equilibration (default=False): if True, then the integration will stop if
                the evaporation rate falls below eps * the initial mass
            equ_threshold: threshold to use for the equilibration termination criterion.
            terminate_on_efflorescence (default=False): if True, then the integration will stop if
                the solvent activity falls below a threshold.
            eff_threshold: threshold to use for the efflorescence termination criterion.
            first_step: size of initial integration step. The subsequent timesteps are determined
                dynamically based on the rate of change and the error tolerance parameter rtol.
        Returns:
            Trajectory of historical droplets showing how it reaches the new state.
        """
        from scipy.integrate import solve_ivp

        events = []
        if terminate_on_equilibration:
            m0 = self.mass
            equilibrated = lambda t,x: np.abs(self.virtual_droplet(x).dmdt) - equ_threshold*m0
            equilibrated.terminal = True
            events += [equilibrated]

        if terminate_on_efflorescence:
            efflorescing = lambda t,x: self.virtual_droplet(x).surface_solvent_activity - eff_threshold
            efflorescing.terminal = True
            events += [efflorescing]

        dxdt = lambda t,x: self.virtual_droplet(x).dxdt
        try:
            with np.errstate(divide='raise', invalid='raise'):
                trajectory = solve_ivp(dxdt, (0, t), self.state, first_step=first_step, rtol=rtol, events=events)
        except Exception as e:
            raise RuntimeError('an error ("%r") occurred during the simulation that we don\'t know how to handle - try running again with a smaller rtol' % e) from None

        self.state = trajectory.y[:,-1]
        return trajectory

    def complete_trajectory(self, trajectory):
        """Get the trajectory of all variables (including dependent ones) from a simulation (i.e.
        the output of UniformDroplet.integrate).

        Args:
            trajectory: the output of UniformDroplet.integrate, which gives the trajectory of independent
                        variables only.
        Returns:
            A pandas dataframe detailing the complete droplet history.
        """

        variables = self.complete_state
        for label in variables:
            variables[label] = np.empty(trajectory.t.size)

        for i,state in enumerate(trajectory.y.T):
            earlier_droplet = self.virtual_droplet(state)
            earlier_state = earlier_droplet.complete_state
            for label,value in earlier_state.items():
                variables[label][i] = value

        variables['time'] = trajectory.t

        return pd.DataFrame(variables)

# ## 3.2. Example: running a droplet simulation from raw python code
# This is a minimal working example of how you might use the previously defined class to run a simulation. You can modify this example to e.g. create scripts to automatically run simulations over varying conditions.

if __name__ == '__main__':
    solution = aqueous_NaCl

    ambient_temperature = 293 # kelvin
    ambient_RH = 0.1
    gas = Atmosphere(ambient_temperature, ambient_RH)

    gravity = np.zeros(3) # ignore body forces like gravity (representing e.g. EDB setups)
    initial_radius = 30e-6 # metres
    initial_mfs = 0.2
    initial_temperature = 293 # kelvin
    initial_velocity = np.zeros(3) # metres/second
    initial_position = np.zeros(3) # metres/second
    droplet = UniformDroplet.from_mfs(solution,
                                      gas,
                                      gravity,
                                      initial_radius,
                                      initial_mfs,
                                      initial_temperature,
                                      initial_velocity,
                                      initial_position)

    time = 100 # seconds
    trajectory = droplet.integrate(time, terminate_on_equilibration=True)
    # Obtain a table giving a history of *all* droplet parameters.
    history = droplet.complete_trajectory(trajectory)

    plt.plot(history['time'], history['radius'])
    plt.xlabel('time / s')
    plt.ylabel('radius / m')
    plt.show()
