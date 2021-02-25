# How to Build Applications with Intel&reg; Software to Get Best Performance on Intel&reg; Hardware
## Intel Programming Tools Are Free and Easy to Install
It makes sense that Intel software products like compilers and libraries will deliver the best performance on Intel hardware, but did you know that they are available for free? The [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html) contains compilers, libraries, and analysis tools to improve the performance of your applications and the open source applications that you use.

You can [download](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html) the toolkit, select your operating system (Linux, Windows, or MacOS) and distribution preference (web or local installation, APT Package Manager, Docker Hub, YUM Package Manager, Zypper Package Manager) from the pulldown menus, and follow the installation instructions. It’s that easy, and did I mention the toolkit is free? (See [Free Intel® Software Development Tools](https://software.intel.com/content/www/us/en/develop/articles/free-intel-software-developer-tools.html) if you don’t believe me.) Add-on toolkits for HPC, IoT, advanced rendering, and AI analytics can also be [downloaded](https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html) for free.
## Intel Programming Tools Are Easy to Use
In addition to being free, Intel software is easy to incorporate into your computing environment. For example, many Linux-based, open-source applications are built using familiar configure, make, and cmake commands. Configuration scripts and makefiles often define the C and C++ compilers using CC and CXX variables, respectively. Replacing the default compilers with the Intel compilers is as easy as setting these variables on the configure or make command-lines, e.g.:
```
$ ./configure CC=icx CXX=icpx
$ make CC=icx CXX=icx
```
Similarly, CMake uses variables to change compilers, e.g.:

    $ CC=icx CXX=icpx cmake /path_to_source

If your project uses a cloud continuous integration (CI) system like GitHub Actions Intel makes that easy too: [Intel® oneAPI CI Samples](https://github.com/oneapi-src/oneapi-ci).
## Intel Software Delivers Great Performance on Intel Hardware

A real open-source application demonstrates how the Intel programming tools improve performance on Intel hardware. The compute-intensive [RELION](https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page) application from the MRC Laboratory of Molecular Biology at the University of Cambridge is used to solve biomolecular structures, [like the COVID-19 virus](https://www.nature.com/articles/s41586-020-2665-2), from cryo-electron microscopy images. The default RELION configuration (which uses the GNU compilers and the FFTW library) is compared to the configuration that uses the Intel oneAPI HPC Toolkit to target Intel hardware.

The RELION developers publish [benchmarks](https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Benchmarks_%26_computer_hardware#Intel_.22Skylake.22_systems_-_RELION_built_for_Intel.C2.AE_AVX-2_with_GCC_7.3) for a standard data set, along with system details and build options for each set of performance results. The results shown in the table below were collected on a cluster of dual-socket Intel® Xeon® Gold 6148 processors (@2.40 GHz with Intel® Hyper-Threading Technology enabled).

|Compiler|Vectorization|# Nodes|# MPI Ranks|Time (hr:min)|Intel Compiler Speedup Over GNU Compiler
|--|--|--|--|--|--
|GNU GCC 7.3|AVX-2|1|8|13:34||
|ICC 2018.3|AVX-2|1|8|8:21|1.6|
|ICC 2018.3|AVX-512|1|8|4:05|3.3|
|GNU GCC 7.3|AVX-2|4|16|3:36||
|ICC 2018.3|AVX-2|4|16|2.10|1.6|
|ICC 2018.3|AVX-512|4|16|1.08|3.3|

The Intel compiler gives better performance than the GNU compiler, especially when the application is compiled to take advantage of the AVX-512 instruction set. RELION is more than 3x faster when built with the Intel compiler.

If you want to unlock the performance potential of Intel hardware, use Intel software products.
