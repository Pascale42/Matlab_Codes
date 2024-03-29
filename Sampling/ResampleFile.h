/*
 * MATLAB Compiler: 3.0
 * Date: Mon Jul 12 20:54:12 2004
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "-v" "ResampleFile" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __ResampleFile_h
#define __ResampleFile_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_ResampleFile(void);
extern void TerminateModule_ResampleFile(void);
extern _mexLocalFunctionTable _local_function_table_ResampleFile;

extern void mlfResampleFile(mxArray * inname,
                            mxArray * outname,
                            mxArray * numchannel,
                            mxArray * resampldown,
                            mxArray * resamplup);
extern void mlxResampleFile(int nlhs,
                            mxArray * plhs[],
                            int nrhs,
                            mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
