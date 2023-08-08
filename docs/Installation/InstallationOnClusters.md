\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Installation on Clusters {#installation_on_clusters}

\tableofcontents

The installation instructions are the same for most systems because we use shell
scripts to set up the environment for each supercomputer. We describe the
generic installation instructions once, and only note special instructions if
necessary for the particular system. If you have already built SpECTRE and just
want to load the modules, source the shell file for your system and run
`spectre_load_modules`.

\note Sample submit scripts for some systems are available in
`support/SubmitScripts`.

## General Instructions

1. Run `export SPECTRE_HOME=/path/to/where/you/want/to/clone`
2. Clone SpECTRE using `git clone SPECTRE_URL $SPECTRE_HOME`
3. Run `cd $SPECTRE_HOME && mkdir build && cd build`
4. Run `. $SPECTRE_HOME/support/Environments/SYSTEM_TO_RUN_ON_gcc.sh`, where
   `SYSTEM_TO_RUN_ON` is replaced by the name of the system as described in the
   relevant section below.
5. If you haven't already installed the dependencies, run
   `export SPECTRE_DEPS=/path/to/where/you/want/the/deps`
   Then run `spectre_setup_modules $SPECTRE_DEPS`. This
   will take a while to finish. Near the end the command will tell you how to
   make the modules available by providing a `module use` command. Make
   sure you are providing an absolute path to `spectre_setup_modules`.
6. Run `module use $SPECTRE_DEPS/modules`
7. Run `spectre_run_cmake`, if you get module loading errors run
   `spectre_unload_modules` and try running `spectre_run_cmake` again. CMake
   should set up successfully.
8. Build the targets you are interested in by running, e.g.
   `make -j4 test-executables`

## Anvil at Purdue University

You should build and run tests on a compute node. You can get a compute node by
running
```
sinteractive -N1 -n 20 -p debug -t 60:00
```
Avoid running `module purge` because this also removes various default modules
that are necessary for proper operation. Instead, use `module
restore`. Currently the tests can only be run in serial, e.g. `ctest -j1`
because all the MPI jobs end up being launched on the same core.

## Cedar, Graham, and Niagara at ComputeCanada

Use `compute_canada` as the `SYSTEM_TO_RUN_ON` in the general instructions.

## Frontera at TACC

Follow the general instructions, using `frontera` for `SYSTEM_TO_RUN_ON`.

Processes running on the head nodes have restrictions on memory use
that will prevent linking the main executables.  It is better to
compile on an interactive node.  Interactive nodes can be requested
with the `idev` command.

For unknown reasons, incremental builds work poorly on frontera.
Running `make` will often unnecessarily recompile SpECTRE libraries.

## Wheeler at Caltech

Follow the general instructions using `wheeler` for `SYSTEM_TO_RUN_ON`, except
you do not need to install any dependencies, so you can skip steps 5 and 6. You
can optionally compile using LLVM/Clang by sourcing `wheeler_clang.sh` instead
of `wheeler_gcc.sh`

If you are running jobs on a Wheeler interactive compute
node, make sure that when you allocate the interactive node using
`srun`, use the `-c <CPUS_PER_TASK>` option to `srun`, and not the `-n
<NUMBER_OF_TASKS>` option.  If you use the `-n <NUMBER_OF_TASKS>`
option and pass the number of cores for NUMBER_OF_TASKS, then you will
get multiple MPI ranks on your node and the run will hang.

## CaltechHPC at Caltech

Follow the general instructions, using `caltech_hpc` for `SYSTEM_TO_RUN_ON`.
When you go to build, you will need to get an interactive node (login nodes
limit the amount of memory accessible to individual users, to below the amount
necessary to build SpECTRE).
To ensure you get an entire node to build on, use the command:
```
srun -t 02:00:00 -N 1 -c 32 --mem=192000 -A sxs -D . --pty /bin/bash
```
being sure to re-source environment files once you get the interactive node
shell.

## Ocean at Fullerton

Follow the general instructions, using `ocean` for `SYSTEM_TO_RUN_ON`,
you do not need to install any dependencies, so you can skip steps 5 and 6.
