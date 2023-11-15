using RCall
using Pkg

# SETUP yours
user_libs = "C:/Users/User/AppData/Local/R/cache/R/renv/library/functional-data-regression-mip-e3349204/R-4.3/x86_64-w64-mingw32"

function set_R_lib_path(project_root)
    
    Pkg.activate(project_root)
    Pkg.instantiate()

    @rput user_libs
    R"""
    .libPaths(c(user_libs,.libPaths()))
    """
end