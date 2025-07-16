#!/bin/bash

#############################################################################
## Identify the candidate drivers using OncodriveFML
## https://oncodrivefml.readthedocs.io/en/latest/configuration.html
#############################################################################

## example files
# Option 1: Create directory and use curl
mkdir -p /path/to/download/folder
curl -L https://github.com/bbglab/oncodrivefml/archive/refs/tags/2.5.0.tar.gz -o /path/to/download/folder/2.5.0.tar.gz
cd /path/to/download/folder

# Option 2: Use wget with output directory option
mkdir -p /path/to/download/folder
wget -P /path/to/download/folder https://github.com/bbglab/oncodrivefml/archive/refs/tags/2.5.0.tar.gz
cd /path/to/download/folder

# Option 3: Direct download with output option
mkdir -p /path/to/download/folder
wget https://github.com/bbglab/oncodrivefml/archive/refs/tags/2.5.0.tar.gz -O /path/to/download/folder/2.5.0.tar.gz
cd /path/to/download/folder
tar xvzf 2.5.0.tar.gz