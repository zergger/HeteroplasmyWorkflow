# icHET Enhanced Version

This is a modified version of icHET: Interactive Visualization of Cytoplasmic Heteroplasmy. This repository optimizes the original icHET method (Phan et al. 2019) through threshold-based enhancements, including the optimization of the `percentage_threshold` and the introduction of new parameters such as `rm_alter_align` and `count_threshold`.

## Table of Contents

- [icHET Enhanced Version](#ichet-enhanced-version)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Step-by-Step Installation](#step-by-step-installation)
  - [Usage](#usage)
  - [Configuration](#configuration)
  - [References](#references)

## Features

- **Threshold-based Enhancements**: Improved the `percentage_threshold` parameter for better accuracy.
- **New Parameters**: Introduced `rm_alter_align` and `count_threshold` for more control over data processing.
- **Optimized Performance**: Enhanced performance and accuracy of the original method.

## Installation

### Dependencies

Ensure the following third-party software dependencies are installed:

- [samtools](http://www.htslib.org/download/) For manipulating high-throughput sequencing data.
- [bwa-meme](https://github.com/kaist-ina/BWA-MEME) An optimized version of BWA for aligning sequences. (Note: Using the original BWA requires code modification)
- [sambamba](https://lomereiter.github.io/sambamba/) A high-performance tool for working with SAM/BAM files.
- **Python 3**: For running the scripts.

### Step-by-Step Installation

1. **Install samtools**:
    \`\`\`bash
    sudo apt-get install samtools
    \`\`\`

2. **Install bwa-meme**:
    Follow the instructions on the [bwa-meme GitHub repository](https://github.com/kaist-ina/BWA-MEME).

3. **Install sambamba**:
    \`\`\`bash
    sudo apt-get install sambamba
    \`\`\`

4. **Install Python 3**:
    \`\`\`bash
    sudo apt-get install python3
    \`\`\`

## Usage

1. Clone this repository:
    \`\`\`bash
    git clone https://github.com/zergger/HeteroplasmyWorkflow.git
    cd HeteroplasmyWorkflow
    \`\`\`

2. Run the main script:
    \`\`\`bash
    python run_hpc.py config.txt read_ids.txt
    \`\`\`

## Configuration

Ensure the `config.txt` file includes the following parameters:

- `PE = 1` for paired-end reads, or `PE = 0` for single-end reads.
- `rm_alter_align = 1` to remove alter alignments (Numts). Set to `0` to disable, but enabling is recommended.

## References

- Phan et al. 2019, *Interactive Visualization of Cytoplasmic Heteroplasmy*.
- [Original icHET repository](https://github.com/vtphan/HeteroplasmyWorkflow)


