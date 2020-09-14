# Contributing

1. Before making a substantial pull request, consider first [filing an issue](https://github.com/broadinstitute/ml/issues) describing the feature addition or change you wish to make.
1. [Get setup](#setup-for-code-contributions)
1. [Follow the coding style](#python-coding-style)
1. [Test your code](#testing)
1. Send a [pull request](https://github.com/broadinstitute/ml/pulls)

## Setup for code contributions

### Get setup for GitHub

Small typos in code or documentation may be edited directly using the GitHub web interface. Otherwise:

1. If you are new to GitHub, don't start here. Instead, work through a GitHub tutorial such as https://guides.github.com/activities/hello-world/.
1. Create a fork of https://github.com/broadinstitute/ml
1. Clone your fork.
1. Work from a feature branch. See the [Appendix](#appendix) for detailed `git` commands.

### Install precommit

[`pre-commit`](https://pre-commit.com/) is a framework for managing and maintaining multi-language pre-commit hooks.

```
# Install pre-commit
pip3 install pre-commit
# Install the git hook scripts by running this within the git clone directory
cd ${HOME}/ml
pre-commit install
```

See [.pre-commit-config.yaml](https://github.com/broadinstitute/ml/blob/master/.pre-commit-config.yaml) for the currently configured pre-commit hooks for ml4cvd.

### Install git-secrets

```git-secrets``` helps us avoid committing secrets (e.g. private keys) and other critical data (e.g. PHI) to our
repositories. ```git-secrets``` can be obtained via [github](https://github.com/awslabs/git-secrets) or on MacOS can be
installed with Homebrew by running ```brew install git-secrets```.

To add hooks to all repositories that you initialize or clone in the future:

```git secrets --install --global```

To add hooks to all local repositories:

```
git secrets --install ~/.git-templates/git-secrets
git config --global init.templateDir ~/.git-templates/git-secrets
```

We maintain our own custom "provider" to cover any private keys or other critical data that we would like to avoid
committing to our repositories. Feel free to add ```egrep```-compatible regular expressions to
```git_secrets_provider_ml4cvd.txt``` to match types of critical data that are not currently covered by the patterns in that
file. To register the patterns in this file with ```git-secrets```:

```
git secrets --add-provider -- cat ${HOME}/ml/git_secrets_provider_ml4cvd.txt
```

### Install pylint

[`pylint`](https://www.pylint.org/) is a Python static code analysis tool which looks for programming errors, helps enforcing a coding standard, sniffs for code smells and offers simple refactoring suggestions.

```
# Install pylint
pip3 install pylint
```

See [pylintrc](https://github.com/broadinstitute/ml/blob/master/pylintrc) for the current lint configuration for ml4cvd.

# Python coding style

Changes to ml4cvd should conform to [PEP 8 -- Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/). See also [Google Python Style Guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md) as another decription of this coding style.

Use `pylint` to check your Python changes:

```bash
pylint --rcfile=${HOME}/ml/pylintrc myfile.py
```

Any messages returned by `pylint` are intended to be self-explanatory, but that isn't always the case.

* Search for `pylint <alphanumeric code>` or `pylint <keyword-code>` for more details on the recommended code change to resolve the lint issue.
* Or add comment `# pylint: disable=<keyword-code>` to the end of the line of code.

# Testing

## Testing of `recipes`

Unit tests can be run in Docker with
```
${HOME}/ml/scripts/tf.sh -T ${HOME}/ml/tests
```
Unit tests can be run locally in a conda environment with
```
python -m pytest ${HOME}/ml/tests
```
Some of the unit tests are slow due to creating, saving and loading `tensorflow` models.
To skip those tests to move quickly, run
```
python -m pytest ${HOME}/ml/tests -m "not slow"
```
pytest can also run specific tests using `::`. For example

```
python -m pytest ${HOME}/ml/tests/test_models.py::TestMakeMultimodalMultitaskModel::test_u_connect_segment
```

For more pytest usage information, checkout the [usage guide](https://docs.pytest.org/en/latest/usage.html).

## Testing of `visualization_tools`

The code in [ml4cvd/visualization_tools](https://github.com/broadinstitute/ml/tree/master/ml4cvd/visualization_tools) is primarily interactive so we add test cases to notebook [test_error_handling_for_notebook_visualizations.ipynb](https://github.com/broadinstitute/ml/blob/master/notebooks/review_results/test_error_handling_for_notebook_visualizations.ipynb) and visually inspect the output of `Cells -> Run all`.

# Appendix

For the ml4cvd GitHub repository, we are doing ‘merge and squash’ of pull requests. So that means your fork does not match upstream after your pull request has been merged. The easiest way to manage this is to always work in a feature branch, instead of checking changes into your fork’s master branch.


## How to work on a new feature

(1) Get the latest version of the upstream repo

```
git fetch upstream
```

Note: If you get an error saying that upstream is unknown, run the following remote add command and then re-run the fetch command. You only need to do this once per git clone.

```
git remote add upstream https://github.com/broadinstitute/ml.git
```

(2) Make sure your master branch is “even” with upstream.

```
git checkout master
git merge --ff-only upstream/master
git push
```

Now the master branch of your fork on GitHub should say *"This branch is even with broadinstitute:master."*.


(3) Create a feature branch for your change.

```
git checkout -b my-feature-branch-name
```

Because you created this feature branch from your master branch that was up to date with upstream (step 2), your feature branch is also up to date with upstream. Commit your changes to this branch until you are happy with them.

(4) Push your changes to GitHub and send a pull request.

```
git push --set-upstream origin my-feature-branch-name
```

After your pull request is merged, its safe to delete your branch!

## I accidentally checked a new change to my master branch instead of a feature branch. How to fix this?

(1) Soft undo your change(s). This leaves the changes in the files on disk but undoes the commit.

```
git checkout master
# Moves pointer back to previous HEAD
git reset --soft HEAD@{1}
```

Or if you need to move back several commits to the most recent one in common with upstream, you can change ‘1’ to be however many commits back you need to go.

(2) “stash” your now-unchecked-in changes so that you can get them back later.

```
git stash
```

(3) Now do the [How to work on a new feature](#how-to-work-on-a-new-feature) step to bring master up to date and create your new feature branch that is “even” with upstream. Here are those commands again:

```
git fetch upstream
git merge --ff-only upstream/master
git checkout -b my-feature-branch-name
```

(4) “unstash” your changes.

```
git stash pop
```
Now you can proceed with your work!
