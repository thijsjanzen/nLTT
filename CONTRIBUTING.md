# Contributing

Awesome that you are reading this!

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

## Branching policy

 * The `master` branch should always build successfully
 * The `development` branch is for developers

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage above 95%
 * Coding style should follow the default style by `lintr`

These are all checked by Travis CI when submitting
a Pull Request. 

## git usage

To get started working on `nLTT` do:

```
git clone https://github.com/thijsjanzen/nLTT.git
```

Development is done on the `develop` branch. 
To download and checkout the `develop` branch, 
first go into the `nLTT` folder (`cd nLTT`), then do:

```
git checkout develop
```

Then the workflow is the common `git` workflow:

```
git pull
git add --all :/
git commit -m "Did something awesome"
git push
```
