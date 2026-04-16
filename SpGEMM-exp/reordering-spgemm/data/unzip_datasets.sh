#!/bin/bash

# Use the following Linux command to unzip all .tar.gz files in a directory:
for file in *.tar.gz; do tar -xzf "$file"; done

# Use the following Linux command to remove the original .tar.gz files after extraction, use:
for file in *.tar.gz; do tar -xzf "$file" && rm "$file"; done
