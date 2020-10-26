# Spherical Earth mode

The spherical Earth mode can be invoked by setting [is_sphere](parameter_list.md#is_sphere) = T. This mode assumes that input/output velocities and depths are defined with the spherical coordinate. Before forward computations, velocities, depths, and densities are adjusted for the cartesian coordinate via the Earth-flattening approximation as follows:

$$
z_{cart} = R_E\ln\left(\frac{R_E}{r}\right) ,
$$
$$
\alpha_{cart} = \alpha_{sphe} \left(\frac{r}{R_E}\right),
$$
$$
\beta_{cart} = \beta_{sphe} \left(\frac{r}{R_E}\right),
$$
and 
$$
\rho_{cart} = \rho_{sphe} \left(\frac{r}{R_E}\right), 
$$

where \\(z\\), \\(\alpha\\), \\(\beta\\), and \\(\rho\\) represent depth, P-wave velocity, S-wave velocity, and density, respectively. \\(R_E\\) is the Earth's radius and \\(r\\) is a distance to the origin.

