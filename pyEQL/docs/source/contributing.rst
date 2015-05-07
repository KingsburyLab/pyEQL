.. _contributing:

Contributing to pyEQL
*********************

Reporting Issues
================

You can report any bugs, packaging issues, feature requests, comments, or questions
using the `issue tracker <https://github.com/rkingsbury/pyEQL/issues>`_ on `github <https://github.com/rsking84/pyeql>`_.

Contributing Code
=================

To contribute bug fixes, documentation enhancements, or new code, please 
fork pyEQL and send us a pull request. It's not as hard as it sounds!

It is **strongly** recommended that you read the following short articles
before starting your work, especially if you are new to the open source community.

* `Open Source Contribution Etiquette <http://tirania.org/blog/archive/2010/Dec-31.html>`_
* `Don't "Push" Your Pull Requests <https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/>`_
* `A Successful Git Branching Model <http://nvie.com/posts/a-successful-git-branching-model>`_

Hacking pyEQL in Six Easy Steps:
---------------------------------

1. `Fork the pyEQL repository <https://help.github.com/articles/fork-a-repo/>`_ on Github

2. Clone your repository to a directory of your choice::

    git clone https://github.com/<username>/pyEQL

3. Create a branch for your work. We loosely follow the branching guidelines
   outlined at http://nvie.com/posts/a-successful-git-branching-model.

   If you are adding **documentation** or **bug fixes**, start with the **master** branch and
   prefix your branch with "fix-" or "doc-" as appropriate::

    git checkout -b fix-myfix master

    git checkout -b doc-mydoc master

   If you are adding a **new feature**, start with the **develop** branch and prefix your
   branch with "feature-"::

    git checkout -b feature-myfeature develop

4. Hack away until you're satisfied.

5. Push your work back to Github::

    git push origin feature-myfeature

6. Create a pull request with your changes. See `this tutorial <https://yangsu.github.io/pull-request-tutorial>`_ for instructions.

Generating Test Cases
=====================

pyEQL has many capabilities that have not been tested thoroughly. You can help
the project simply by using pyEQL and comparing the output to experimental data
and/or more established models. Report back your results on the 
`issue tracker <https://github.com/rkingsbury/pyEQL/issues>`_.

Even better, write up an automated test case (see the tests/ directory for examples).

Making a Donation
=================

If you'd like to leave a 'tip' for the project maintainer to support the time and effort
required to develop pyEQL, simply send it via Paypal to RyanSKingsbury@alumni.unc.edu