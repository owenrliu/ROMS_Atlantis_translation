# ROMS_Atlantis_translation
Code to translate ROMS model outputs into useable files for Atlantis models

Porting to GOA highlights one limitations of our original code: it has a hard time with large input data (e.g. NEP10K ROMS).
What changed from the “West Coast” code:
1. Grid information (coordinates etc.) for NEP is stored in a grid file – which is roms2 here.
2. Because of point 1, we have to use a custom version of angstroms::romshcoords() to extract depth at rho points.
3. In the interpolation function, Pulling vars with hyper_tibble(select_var=variable, xi_rho = between(xi_rho, min_xi, max_xi), eta_rho = between(eta_rho, min_eta, max_eta)), or else it chokes the memory.
4. We cannot subset the depth raster due to the way we index its cells to coincide with rho_idx, but we can subset the depth data frame with filter(between(xi_rho, min_xi, max_xi) & between(eta_rho, min_eta, max_eta))
5. In the interpolate_velocity function, changed the left join to inner join to make sure that there are no NAs.

(side note: angstroms::romshcoords() has code to get h at v and u points, but my attempts at making it run failed. We have been assuming that h at u and v points is the same as the closest rho point, but we might want to fix that as there is a 5 km gap between rho and u and v, and there may be substantial differences in depth in areas like WCHG and the continental slope. An evolution of this code might be a horizontal interpolation to get an estimate of depth between rho points at the location of u and v points)
