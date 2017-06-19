%--------------------------------------------------------------------------
% Script for mex compilation of C files. Use with care. 
%
% With great power comes great responsibility.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
             
%% Alternative build, with command lines (use -v option for verbose mode)
% mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ispolycw.c
%%
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_cr3bp_rel.c  lib/libgsl.a lib/libgslcblas.a  C/ode78_rel.c  C/cr3bp_derivatives_rel.c  C/custom.c C/custom_ode.c
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_bcp_event_rel.c   lib/libgsl.a lib/libgslcblas.a  C/ode78_rel.c C/custom_odezero.c C/cr3bp_derivatives_rel.c C/custom.c C/custom_ode.c
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_cr3bp_event_rel.c lib/libgsl.a lib/libgslcblas.a  C/ode78_rel.c C/custom_odezero.c C/cr3bp_derivatives_rel.c C/custom.c C/custom_ode.c
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_bcp_rel.c    lib/libgsl.a lib/libgslcblas.a  C/ode78_rel.c  C/custom_odezero.c  C/cr3bp_derivatives_rel.c  C/custom.c  C/custom_ode.c

mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_cr3bp.c  lib/libgsl.a lib/libgslcblas.a  C/ode78.c  C/cr3bp_derivatives.c  C/custom.c C/custom_ode.c 
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_bcp_event.c   lib/libgsl.a lib/libgslcblas.a  C/ode78.c C/custom_odezero.c C/cr3bp_derivatives.c C/custom.c C/custom_ode.c
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_cr3bp_event.c lib/libgsl.a lib/libgslcblas.a  C/ode78.c C/custom_odezero.c C/cr3bp_derivatives.c C/custom.c C/custom_ode.c
mex  -I/lib/gsl -I/lib CFLAGS="\$CFLAGS -std=c99" -I/C  C/ode78_bcp.c    lib/libgsl.a lib/libgslcblas.a  C/ode78.c  C/custom_odezero.c  C/cr3bp_derivatives.c  C/custom.c  C/custom_ode.c


%% Build C files via mex()
% Define some paths
gslpath         = '-I/lib/gsl';
incpath         = '-I/C';
libgslpath      = 'lib/libgsl.a';
libgslcblaspath = 'lib/libgslcblas.a';

% Build ode78_cr3bp.c
mex(gslpath, 'ode78_cr3bp.c', libgslpath,... 
                              libgslcblaspath,... 
                              'C/ode78.c',... 
                              'C/cr3bp_derivatives.c',... 
                              'C/custom.c',... 
                              'C/custom_ode.c');
%
% Build ode78_cr3bp_event.c
mex(gslpath, 'ode78_cr3bp_event.c', libgslpath,...
                                    libgslcblaspath, ...
                                    'C/ode78.c',... 
                                    'C/custom_odezero.c', ...
                                    'C/cr3bp_derivatives.c', ...
                                    'C/custom.c', ...
                                    'C/custom_ode.c');
                                
% Build ode78_bcp.c                               
mex(gslpath, 'ode78_bcp.c', libgslpath,... 
                              libgslcblaspath,... 
                              'C/ode78.c',... 
                              'C/cr3bp_derivatives.c',... 
                              'C/custom.c',... 
                              'C/custom_ode.c');
                          
% Build ode78_bcp_event.c                          
mex(gslpath, 'ode78_bcp_event.c', libgslpath,...
                                    libgslcblaspath, ...
                                    'C/ode78.c',... 
                                    'C/custom_odezero.c', ...
                                    'C/cr3bp_derivatives.c', ...
                                    'C/custom.c', ...
                                    'C/custom_ode.c');