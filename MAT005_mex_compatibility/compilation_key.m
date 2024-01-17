function key = compilation_key()
    key = containers.Map;
    key("trajectory_integration_module") = {"FOR001_modules/trajectory_integration_module.f", "C__001_wrappers/trajectory_integration_module.c"};
    % key("fly_aqs")                       = {"FOR001_modules/trajectory_integration_module.f", "C__001_wrappers/fly_aqs.c"};
    key("integrate_trajectory_lite")     = {"FOR001_modules/trajectory_integration_module.f", "C__001_wrappers/integrate_trajectory_lite.c"};
    % key("iterate_laplace")               = {"FOR001_modules/potential_calculation.f", "C__001_wrappers/iterate_laplace.c"};
    % key("refined_laplace")               = {"FOR001_modules/potential_calculation.f", "C__001_wrappers/refined_laplace.c"};
    key("ray_optics_2D")                 = {"FOR001_modules/ion_optics.f", "C__001_wrappers/ray_optics_2D.c"};
    key("ray_optics_ensemble")           = {"FOR001_modules/ion_optics.f", "C__001_wrappers/ray_optics_ensemble.c"};
    key("interpolate_monotonic_1d")      = {"FOR001_modules/utils.f", "C__001_wrappers/interpolate_monotonic_1d.c"};
    key("linspace_fortran")              = {"FOR001_modules/utils.f", "C__001_wrappers/linspace_fortran.c"};
    key("fly_ensemble")                  = {"FOR001_modules/trajectory_integration_module.f", "C__001_wrappers/fly_ensemble.c"};
    key("fly_cloud")                     = {"FOR001_modules/trajectory_integration_module.f", "C__001_wrappers/fly_cloud.c"};
end