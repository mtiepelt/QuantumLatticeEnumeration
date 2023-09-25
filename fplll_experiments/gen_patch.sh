#!/usr/bin/env bash

# https://stackoverflow.com/a/45353002

# Enumerate the files modified on the desired commit
for file in $(git diff --name-only HEAD~1); do
    # Generate the name of the patch file: replac '/' with '_'
    # in the paths of the modified files and add the .patch termination
    patch=${file//\//_}.patch
    # Create one patch for each modified file
    git diff HEAD~1 -- $file > $patch
done
