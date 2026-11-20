<!-- Delete the lines that do not apply, including these comments. The description becomes the squash-merge commit message on master. -->

## What this changes

<!-- One or two sentences: what changes and why. Jira keys such as [VS-NNNN] go in the PR title. -->

- Issue: <!-- #NNN if there is one. Write "Fixes #NNN" only if merging should close it. -->
- User-visible changes (output format, VCF header, argument defaults): none

## Testing

- [ ] Tests added or updated (bug fixes need a regression test)
- [ ] Ran the affected tests locally after `git lfs pull`, for example `./gradlew test --tests '*MyToolIntegrationTest'`

<!-- If expected-output files changed: which tool, which fields, and why the new values are right.
     Cloud tests need repository secrets, so on a PR from a fork the cloud job passes without running them.
     Every other partition, including the Docker ones, runs on every PR. -->

## Checklist

- [ ] Code follows the [developer guidelines](https://github.com/broadinstitute/gatk#general-guidelines-for-gatk4-developers)
- [ ] Tool and argument docs are updated (class javadoc, `@CommandLineProgramProperties` and `@Argument` doc strings generate the tool documentation)
- [ ] Test files over 100 KB are under `src/test/resources/large/` (Git LFS picks them up from .gitattributes; do not run `git lfs track`)
- [ ] Changes to `Dockerfile`, `scripts/docker/`, `scripts/gatkcondaenv.yml.template` or dependencies in `build.gradle` are described above
