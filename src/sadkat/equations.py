# # A. Overview of model
# ## A.1. System of equations

# $\renewcommand\vec[1]{\mathbf{#1}}$
# The equations of motion in Xie et al (2007) read:
#
# \begin{align}
#   \frac{d r_p}{dt} =& \;
#   \mathrm{Sh}
#   \frac{C D_\infty}{\rho_p r_p} \frac{M_v p}{R T_\infty}
#   \ln{\left( \frac{p - p_{va}}{p - p_{v^\infty}} \right)}
#   = f_1(r_p, T_p, \vec{V}_p),
#   \\
#   \frac{d T_p}{dt} =& \;
#   4\pi r_p^2
#   \, K_g \frac{T_\infty - T_p}{c_p m_p r_p} \mathrm{Nu}
#   - \frac{L_v I}{c_p m_p}
#   - 4\pi r_p^2 \frac{\Gamma (T_p^4 - T_\infty^4)}{c_p m_p}
#   \\ =& \;
#   3 K_g \frac{T_\infty - T_p}{c_p \rho_p r_p^2} \mathrm{Nu}
#   + \frac{L_v \dot{m}_p}{c_p m_p}
#   - \frac{3 \Gamma (T_p^4 - T_\infty^4)}{c_p \rho_p r_p}
#   = f_2(r_p, T_p, \vec{V}_p),
#   \\
#   \frac{d \vec{V}_p}{dt} =& \;
#   \vec{g} \left( 1 - \frac{\rho_p}{\rho_g} \right)
#   - \frac{3 C_d \rho_g |\vec{V}_p - \vec{V}_g| (\vec{V}_p - \vec{V}_g)}{8 \rho_p r_p}
#   = f_3(r_p, T_p, \vec{V}_p),
#   \\
#   \frac{d \vec{x}_p}{dt} =& \;
#   \vec{V}_p
#   = f_4(\vec{V}_p).
# \end{align}
#
# These equations implicitly assume a further relationship to solve these.
# To see this, we divide both the numerator/denominator inside the logarithm in the first equation through by $p$ to obtain
#
# \begin{equation}
#   \frac{p - p_{va}}{p - p_{v^\infty}}
#   =
#   \frac{p (1 - \chi_w(r_p^+))}{p - p_{v^\infty}}
#   =
#   \frac{1 - \chi_w(r_p^+)}{1 - \chi_w(\infty)},
# \end{equation}
#
# which requires a further equation for $\chi_w(r_p^+)$ (the water/solvent mole fraction) for closure.
# Xie et al (2007) assume Raoult's law, i.e.
#
# \begin{equation}
#   p_{va} = \chi_w(r_p^-) \, p_{va,s}(T_p).
# \end{equation}
#
# Mechanical balance ensures that $p_{va}$ is continous across the liquid-vapour phase boundary, and so $p_{va} = \chi_w(r_p^+) \, p$ using the ideal gas law.
#
# The equilibrium vapour pressure for the solvent at the phase boundary can be calculated using the Clausius-Clapeyron relation, leading to Buck's equation for $p_{va}(T_p)$.

# ## A.2. Constant density assumption?
#
# The first equation implicitly assumes that the density is constant.
# To see this, we rewrite this equation as the following two subequations:
#
# \begin{align}
#   \frac{d r_p}{dt} =& \;
#   \frac{1}{4\pi r_p^2 \rho_p} \frac{d m_p}{dt},
#   \\
#   \frac{d m_p}{dt} =& \;
#   \mathrm{Sh} \,
#   4\pi r_p \, C D_\infty \, \frac{M_v p}{R T_\infty}
#   \ln{\left( \frac{p - p_{va}}{p - p_{v^\infty}} \right)}.
# \end{align}
#
# Next, consider the expression for total droplet mass:
#
# \begin{equation}
#   m_p = 4\pi \int_0^{r_p} \rho_p(r') r'^2 \, dr'.
# \end{equation}
#
# When the density is uniform inside the droplet $\rho_p(r') = \rho_p$ then we obtain
#
# \begin{equation}
#   m_p = \frac{4\pi r_p^3}{3} \rho_p,
# \end{equation}
#
# and so the time-derivative yields
#
# \begin{equation}
#   \frac{d m_p}{dt}
#   =
#   4\pi r_p^2 \rho_p \frac{d r_p}{dt}
#   +
#   \frac{4\pi r_p^3}{3} \frac{d \rho_p}{dt}.
# \end{equation}
#
# If we make the density constant in time then this yields the expression for $\dot{r}_p$ above.
# It is possible this is approximately valid if the density only varies slowly.
# More generally, if we have a varying density profile then we obtain
#
# \begin{equation}
#   \frac{d m_p}{dt}
#   =
#   4\pi r_p^2 \rho_p(r_p) \frac{d r_p}{dt}
#   +
#   4\pi \int_0^{r_p} \frac{\partial \rho_p(r')}{\partial t} r'^2 \, dr'.
# \end{equation}
