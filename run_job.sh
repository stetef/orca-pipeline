#!/bin/bash

shopt -s extglob nullglob

# Examples:
#   ./run_job.sh ZnHis4.in
#     Generate and submit the PBS script.
#   ./run_job.sh --dry-run ZnHis4.in
#     Generate the PBS script but skip qsub submission.
#   ./run_job.sh -n ZnHis4.in
#     Same as --dry-run.

dry_run=0
if [[ "$1" == "--dry-run" || "$1" == "-n" ]]; then
  dry_run=1
  shift
fi

if [[ -z "$1" ]]; then
  echo "Usage: $0 [--dry-run|-n] <input.in>" >&2
  exit 2
fi

if [[ ! -f "$1" ]]; then
  echo "Input file not found: $1" >&2
  exit 2
fi

basename="${1%.*}"

# get number of procs from ORCA input
nprocs=$(awk 'BEGIN{IGNORECASE=1}
  /^!/{
    if (match($0, /PAL[[:space:]]*([0-9]+)/, m)) {print m[1]; exit}
  }
  /^%pal/{inpal=1}
  inpal && /nprocs/{
    if (match($0, /nprocs[[:space:]]*([0-9]+)/, m)) {print m[1]; exit}
  }
  inpal && /^end/{inpal=0}
' "$1")
nprocs=${nprocs:-1}
nodes=$nprocs

cat > generated-${basename}-orca.script <<EOF
#!/bin/bash
#PBS -l nodes=1:ppn=${nodes:=1}
#PBS -S /bin/bash
#PBS -N orca-${basename}
#PBS -q workq
#PBS -V

# Load modules
if [[ -f /etc/profile.d/modules.sh ]]; then
  # Initialize environment-modules for non-interactive shells
  source /etc/profile.d/modules.sh
fi
module purge
module load intel/2021.2.0
module load openmpi/4.0.5-intel-2021.2.0
module list

shopt -s extglob nullglob

export PBS_O_WORKDIR=\${PBS_O_WORKDIR:-\$PWD}
export PBS_O_HOST=\${PBS_O_HOST:-\$(hostname)}
export PBS_O_PATH=\${PBS_O_PATH:-\$PATH}

# Reduce UCX verbosity
export UCX_LOG_LEVEL=error

export ORCA_HOME=/opt/bin
export PATH="\$ORCA_HOME:\$PATH"
ORCA_EXEC="\$ORCA_HOME/orca"
if [[ ! -x "\$ORCA_EXEC" ]]; then
  echo "ORCA not found at \$ORCA_EXEC" >&2
  exit 1
fi

logfile=\$PBS_O_WORKDIR/${basename}-orca.log
tdir=\$(mktemp -d /home/stetef/working-station/scratch/orca-run-${basename}__XXXXXX)

trap '
echo "Job terminated from outer space!" >> \$logfile
rm -rf \$tdir
exit
' TERM 

cp \$PBS_O_WORKDIR/$1 \$tdir
for f in "\$PBS_O_WORKDIR"/*.gbw "\$PBS_O_WORKDIR"/*.pot; do
  [[ -e "\$f" ]] || continue
  cp "\$f" "\$tdir"
done
cd \$tdir

# Ensure mpirun exists if module only provides mpiexec
MPIEXEC=\$(command -v mpirun || command -v mpiexec || true)
if [[ -z "\$MPIEXEC" ]]; then
  echo "mpirun/mpiexec not found in PATH" >> \$logfile
  exit 1
fi
MPI_DIR=\$(dirname "\$MPIEXEC")
if [[ ! -x "\$MPI_DIR/mpirun" && -x "\$MPI_DIR/mpiexec" ]]; then
  mkdir -p "\$tdir/mpi-bin"
  ln -s "\$MPI_DIR/mpiexec" "\$tdir/mpi-bin/mpirun"
  export PATH="\$tdir/mpi-bin:\$MPI_DIR:\$PATH"
  echo "Created mpirun shim -> mpiexec in \$tdir/mpi-bin" >> \$logfile
fi

echo "Job started from \${PBS_O_HOST}, running on \$(hostname) in \$tdir using \$ORCA_EXEC" > \$logfile
echo "--- MODULES ---" >> \$logfile
module list 2>> \$logfile
echo "--- MODULE SHOW (openmpi/4.0.5-intel-2021.2.0) ---" >> \$logfile
module show openmpi/4.0.5-intel-2021.2.0 2>> \$logfile
echo "--- MPI ---" >> \$logfile
echo "PATH=\$PATH" >> \$logfile
which mpirun >> \$logfile 2>&1
which mpiexec >> \$logfile 2>&1

"\$ORCA_EXEC" "$1" 1>>"\$logfile" 2>&1

cp !(*.inp|*.tmp*) "\$PBS_O_WORKDIR"/
rm -rf \$tdir

EOF

if [[ $dry_run -eq 1 ]]; then
  echo "Dry run: generated ${basename}-orca.script (qsub skipped)"
else
  qsub -j oe -o generated-${basename}-orca.script.out generated-${basename}-orca.script
fi