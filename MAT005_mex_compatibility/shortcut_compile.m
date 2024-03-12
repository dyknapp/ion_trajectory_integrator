function shortcut_compile
!gfortran -fc-prototypes -fsyntax-only -cpp E:\ion_trajectory_integrator\FOR001_modules\integrator.f > C__001_wrappers/nbody.h
!"E:\Program Files\Intel\OneAPI\setvars.bat" intel64 vs2022 & ifort /fpp /c /dll /O3 /Qparallel E:\ion_trajectory_integrator\FOR001_modules\integrator.f -o ./FOR001_modules/nbody.o
mex("C__001_wrappers/nbody.cpp", sprintf("FOR001_modules/%s.o", "nbody"), ...
            '-outdir', 'FOR001_modules', '-output', "nbody")

!gfortran -fc-prototypes -fsyntax-only -cpp E:\ion_trajectory_integrator\FOR001_modules\integrator.f > C__001_wrappers/nbody4.h
!"E:\Program Files\Intel\OneAPI\setvars.bat" intel64 vs2022 & ifort /fpp /c /dll /O3 /Qparallel E:\ion_trajectory_integrator\FOR001_modules\integrator.f -o ./FOR001_modules/nbody4.o
mex("C__001_wrappers/nbody4.cpp", sprintf("FOR001_modules/%s.o", "nbody4"), ...
            '-outdir', 'FOR001_modules', '-output', "nbody4")
end