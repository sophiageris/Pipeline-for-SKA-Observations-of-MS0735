#!/bin/bash
#SBATCH --job-name=sequential_job
#SBATCH --output=output_%A_%j.out
#SBATCH --error=error_%A_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --time=00:30:00
#SBATCH --partition=quicktest

# Enable nullglob to handle patterns with no matches
shopt -s nullglob

FILES=()

# Search for files with _thermal and .sh, and add them to the array
for file in *_thermal*.sh; do
    FILES+=("$file")
done

# Search for files with _nonthermal and .sh, and add them to the array
for file in *_nonthermal*.sh; do
    FILES+=("$file")
done

# Disable nullglob
shopt -u nullglob

# Check if there are any files to process
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No eligible files found."
    exit 1
fi

for file in "${FILES[@]}"; do
    echo "$file"
    current_file="$(realpath "$file")"  # Store the full file path
    echo $current_file
    base_name="${current_file%.*}"

    source /nfs/scratch/projects/icm/AMI_git/profile/sourceme.sh
    /nfs/scratch/projects/icm/AMI_git/profile/profile.linux < "$current_file" > "${base_name}_profile_2.out"
done
