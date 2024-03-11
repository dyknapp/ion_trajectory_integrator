function shortcut_compile
!gfortran -fc-prototypes -fsyntax-only -cpp C:\Users\h2p_l\Desktop\ion_trajectory_integrator\FOR001_modules\integrator.f > C__001_wrappers/nbody.h
!"C:\Program Files (x86)\Intel\OneAPI\setvars.bat" intel64 vs2022 & ifort /fpp /c /dll /O3 /Qparallel C:\Users\h2p_l\Desktop\ion_trajectory_integrator\FOR001_modules\integrator.f -o ./FOR001_modules/nbody.o
mex("C__001_wrappers/nbody.cpp", sprintf("FOR001_modules/%s.o", "nbody"), ...
            '-outdir', 'FOR001_modules', '-output', "nbody")

!gfortran -fc-prototypes -fsyntax-only -cpp C:\Users\h2p_l\Desktop\ion_trajectory_integrator\FOR001_modules\integrator.f > C__001_wrappers/nbody4.h
!"C:\Program Files (x86)\Intel\OneAPI\setvars.bat" intel64 vs2022 & ifort /fpp /c /dll /O3 /Qparallel C:\Users\h2p_l\Desktop\ion_trajectory_integrator\FOR001_modules\integrator.f -o ./FOR001_modules/nbody4.o
mex("C__001_wrappers/nbody4.cpp", sprintf("FOR001_modules/%s.o", "nbody4"), ...
            '-outdir', 'FOR001_modules', '-output', "nbody4")
end