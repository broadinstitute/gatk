# Tracking the provenance of Work performed on AoU Callsets

## Overview/Background
This document describes how we will track the provenance of WDLs used during the generation of the AoU Foxtrot and later releases. We wish to be able to track the provenance of the WDLs used at all stages of the generation of the Foxtrot release in order to be able to determine the state of the code used at a particular point in time.

As we have discussed, we had wanted to rely on git tagging to track the provenance, but due to issues with dockstore, we have decided to use long-lived one-off branches as a substitute for tags. Note that in the Echo release we used a specific branch for all Echo-related work but this proved to be very cumbersome and hard to use in practice as we were managing making changes to the Echo call set branch and the ah_var_store itself.

## Steps:

For any WDL that is to be run for the FoxTrot release:

1. Create a branch off of ah_var_store at the point you are ready to run the WDL. 
   - Foxtrot release branches should be of the form: `FoxtrotRelease_20250528_1` (where 20250528 is the calendar date in the format YYYYMMDD and `_1` uniqueifies the branch if you branch more than once per date (e.g. `_2`, `_3`, …)
1. Note the githash of the branch IMMEDIATELY as it was created. This is the hash that also exists on ah_var_store and this will be used to keep track of exactly what version of the WDL was used.
2. Add the branch name to the .dockstore.yml in the branch you have just created / are using.
1. Commit and push the branch to GitHub
1. In Terra, run that branched version of the WDL (by selecting the branch from the workflow drop down)
1. Document your work!
   - Include the branch name and the original githash of the branch in the Terra workflow comment (when launching the WDL, or after the fact by editing it)
   - Indicate the branch name and the original githash of the branch in the Jira for the task you ran
   - Update the Foxtrot tracking document to inclue the branch name and the original githash of the branch.
1. After your workflow has completed, you can delete the branch from GitHub.

