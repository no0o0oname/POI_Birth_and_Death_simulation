# Understanding the original code
- compile with `cygwin`, but debug with `mingw`
  - since `mingw` lacks some c files in the code, better to build with `cygwin`
  - but there appears to be some weird problems in `cygwin` debugging relating to _the path of source c files_
- Since it is hard to know the length of C array pointers, the array `pointer` along with its `max_length` are usually used together.
- kernels are categorized into two types:
  - `ckernel`s: connectivity kernels; 
  - `okernel`s: other kernels, including global kernels
- activity
  - `pactivity`: process activity, the rate parameter of exponential distribution, controls how soon (the time interval before) the next process event will occur, equal to the _multiplication_ of all `cactivity`s involved in this process
  - `cactivity`: kernel activity, related to each kernel, equal to the _sum_ of the intensity of every kernel function at this location
  
# Modifications to code
- remove torus topology, (perhaps)