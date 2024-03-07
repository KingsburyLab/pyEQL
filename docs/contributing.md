(contributing)=

# Contributing to pyEQL

## Reporting Issues

You can help the project simply by using pyEQL and comparing the output to experimental data and/or other models and tools.
If you encounter any bugs, packaging issues, feature requests, comments, or questions, please report them
using the [issue tracker](https://github.com/KingsburyLab/pyEQL/issues) on [github](https://github.com/KingsburyLab/pyeql).

```{tip}
Please don't forget to include the closed issues in your search. Sometimes a solution was already reported, and the problem is considered solved.
```

New issue reports should include information about your programming environment (e.g., operating system, Python version) and steps to reproduce the problem. Please try also to simplify the reproduction steps to a very minimal example that still illustrates the problem you are facing. By removing other factors, you help us to identify the root cause of the issue.

## Documentation Improvements

You can help improve `pyEQL` docs by making them more readable and coherent, or
by adding missing information and correcting mistakes.

`pyEQL` documentation uses [Sphinx](https://www.sphinx-doc.org/en/master/) as its main documentation compiler.
This means that the docs are kept in the same repository as the project code, and
that any documentation update is done in the same way was a code contribution.

```{tip}
Please notice that the [GitHub web interface] provides a quick way of
propose changes in `pyEQL`'s files. While this mechanism can
be tricky for normal code contributions, it works perfectly fine for
contributing to the docs, and can be quite handy.

If you are interested in trying this method out, please navigate to
the `docs` folder in the source [repository], find which file you
would like to propose changes and click in the little pencil icon at the
top, to open [GitHub's code editor]. Once you finish editing the file,
please write a message in the form at the bottom of the page describing
which changes have you made and what are the motivations behind them and
submit your proposal.
```

When working on documentation changes in your local machine, you can
compile them using [tox] :

```
tox -e docs
```

and use Python's built-in web server for a preview in your web browser
(`http://localhost:8000`):

```
python3 -m http.server --directory 'docs/_build/html'
```

## Contributing Code

To contribute bug fixes, documentation enhancements, or new code, please fork pyEQL and send us a pull request. It's not as hard as it sounds! Beginning with version 0.6.0, we follow the [GitHub flow](https://docs.github.com/en/get-started/quickstart/github-flow) workflow model.

The [Scientific Python Guide](https://learn.scientific-python.org/development/guides/) is also an excellent technical reference for new and longtime developers.

### Submit an issue

Before you work on any non-trivial code contribution it's best to first create
a report in the [issue tracker](https://github.com/KingsburyLab/pyEQL/issues) to start a discussion on the subject. This often provides additional considerations and avoids unnecessary work.

### Hacking pyEQL, step by step

1. [Fork the pyEQL repository](https://help.github.com/articles/fork-a-repo/) on Github

2. Clone your repository to a directory of your choice:

   ```
   git clone https://github.com/<username>/pyEQL
   ```

3. Install the package and the test dependencies by running the following command from the repository directory:

   ```
   pip install -e '.[testing]``
   ```

4. Create a branch for your work. Preferably, start your branch name with "feature-", "fix-", or "doc-" depending on whether you are contributing **bug fixes**, **documentation** or a **new feature**, e.g.
   prefix your branch with "fix-" or "doc-" as appropriate:

   ```
   git checkout -b mybranch
   ```

   or

   ```
   git checkout -b doc-mydoc
   ```

   or

   ```
   git checkout -b feature-myfeature
   ```

5. Make changes to the code until you're satisfied.

6. Push your work back to Github:

   ```
   git push origin feature-myfeature
   ```

7. Create a pull request with your changes. See [this tutorial](https://yangsu.github.io/pull-request-tutorial) for instructions.

## Guidelines

Please abide by the following guidelines when contributing code to `pyEQL`:

- All changes you make to quacc should be accompanied by unit tests and should not break existing tests. To run the full test suite, run `pytest tests/` from the repository directory.

- Code coverage should be maintained or increase. Each PR will report code coverage after the tests pass, but you can check locally using [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/), by running `pytest --cov tests/`

- All code should include type hints and have internally consistent documentation for the inputs and outputs.

- Use Google style docstrings

- Lint your code with [`ruff`](https://github.com/astral-sh/ruff) by running `ruff check --fix src/` from the repo directory. Alternatively, you can install the `pre-commit` hooks by running `pre-commit install` from the repository directory. This will prevent committing new changes until all linting errors are fixed.

- Update the `CHANGELOG.md` file.

- Ask questions and be open to feedback!

## Documentation

Improvements to the documentation are most welcome! Our documentation system uses `sphinx` with the [Materials for Sphinx](https://bashtage.github.io/sphinx-material/) theme. To edit the documentation locally, run `tox -e autodocs` from the repository root directory. This will serve the documents to http://localhost:8000/ so you can view them in your web browser. When you make changes to the files in the `docs/` directory, the documentation will automatically rebuild and update in your browser (you might have to refresh the page to see changes).

## Changelog

We keep a `CHANGELOG.md` file in the base directory of the repository. Before submitting your PR, be sure to update the `CHANGELOG.md` file under the "Unreleased" section with a brief description of your changes. Our `CHANGELOG.md` file lossely follows the [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format, beginning with `v0.6.0`.
