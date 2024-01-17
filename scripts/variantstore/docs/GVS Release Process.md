Steps for creating a new version of GVS:

1. Create a *pre*-release branch off the main GVS branch. This branch can be called anything, it will soon be deleted after being merged back to the main GVS branch. This branch is created as a separate step from the release branch to make cumulative changes to the `CHANGELOG.md`.
   - Update `CHANGELOG.md` as appropriate for the new version. Decide here what the new semantic version number will be and describe added / changed / fixed / removed functionality ([here's a good guide](https://common-changelog.org/)).
   - PR, get thumbs, merge back to main GVS branch.
1. Based at the pre-release commit on the main GVS branch, create a release branch named like `gvs_<major>.<minor>.<patch>` corresponding to the new semantic version just added to the `CHANGELOG.md`. The name of this branch does matter because it will be long-lived and will be referenced by the `dockstore.yml`. On this branch:
   1. Update the `version` declaration in `GetToolVersions` within `GvsUtils.wdl` to the name of the new release branch.
   1. Update `.dockstore.yml` for `GvsBeta` and `GvsExtractCohortFromSampleNames` to add the new release branch.
   1. Commit these changes and push the new branch to origin e.g. `git push --set-upstream origin <new branch name>`
1. Update the `GvsBeta` and `GvsExtractCohortFromSampleNames` workflows in Dockstore to make the new branch the default version.

Done!
