# Spun-up configurations

## NumericalOcean1-JRA55_v0.1.0
As built by `forced_ocean_simulation(ORCA1(GPU()))` with default arguments.

| | |
|---|---|
| **Grid** | eORCA1 tripolar (``1\degree``), 70 z-star levels to 5500 m (1.5 m at the surface) |
| **Time step** | 90 minutes, split Runge--Kutta 3 ([Silvestri et al., 2026](https://doi.org/10.22541/essoar.15002225/v1)), split-explicit free surface |
| **Equation of state** | TEOS-10 (Conservative Temperature and Absolute Salinity) |
| **Momentum advection** | WENO vector invariant, 5th order ([Silvestri et al., 2024](https://doi.org/10.1029/2023MS004130)) |
| **Tracer advection** | WENO, 7th order |
| **Vertical mixing** | CATKE ([Wagner et al., 2025](https://doi.org/10.1029/2024MS004522)) |
| **Background diffusivity** | κ = max(2×10⁻⁶, 10⁻⁵ \|sin φ\|) m² s⁻¹ ([Henyey et al., 1986](https://doi.org/10.1029/JC091iC07p08487)), ν = 3×10⁻⁵ m² s⁻¹ |
| **Mesoscale eddies** | GM ([Gent & McWilliams, 1990](https://doi.org/10.1175/1520-0485(1990)020%3C0150:IMIOCM%3E2.0.CO;2)) and Redi ([Redi, 1982](https://doi.org/10.1175/1520-0485(1982)012%3C1154:OIMBCR%3E2.0.CO;2)) on triads ([Griffies et al., 1998](https://doi.org/10.1175/1520-0485(1998)028%3C0805:IDIAZC%3E2.0.CO;2)), κ_GM = 800 m² s⁻¹, κ_Redi = 1200 m² s⁻¹ |
| **Horizontal viscosity** | biharmonic, ν₄ = (cell area)² / 50 days |
| **Sea ice** | ClimaSeaIce, 4 thickness categories, EVP dynamics ([Kimmritz et al., 2017](https://doi.org/10.1016/j.ocemod.2017.05.006)), landfast basal stress ([Lemieux et al., 2015](https://doi.org/10.1002/2014JC010678)), incremental remapping ([Lipscomb & Hunke, 2004](https://doi.org/10.1175/1520-0493(2004)132%3C1341:MSITUI%3E2.0.CO;2)) |
| **Atmosphere and runoff** | JRA55-do 1958--2018 ([Tsujino et al., 2018](https://doi.org/10.1016/j.ocemod.2018.07.002)) |
| **Air-sea fluxes** | COARE similarity theory ([Edson et al., 2013](https://doi.org/10.1175/JPO-D-12-0173.1)) |
| **Ice-ocean heat flux** | three-equation ([Holland & Jenkins, 1999](https://doi.org/10.1175/1520-0485(1999)029%3C1787:MTIOIA%3E2.0.CO;2)) |
| **Shortwave penetration** | two-band, SeaWiFS chlorophyll monthly climatology (1998–2010) |
| **Salinity restoring** | to monthly WOA, piston velocity 0.254 m day⁻¹ |
| **Initial conditions** | WOA annual temperature and salinity, ECCO4 sea ice (January 1993), ocean at rest |
| **Spin-up** | 5 cycles of 1958--2018 JRA55 |
