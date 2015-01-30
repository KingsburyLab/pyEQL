.. _contributing:

Contributing to pyEQL
*********************

Reporting Issues
================

You can report any bugs, packaging issues, feature requests, comments, or questions
using the `issue tracker <URL>`_ on `github <https://github.com/rsking84/pyeql>`_.

Contributing Code
=================

To contribute bug fixes, documentation enhancements, or new code, please 
fork pyEQL and send us a pull request. It's not as hard as it sounds!

Hacking pyEQL in Five Easy Steps:
---------------------------------

1. Clone the repository::

    git clone https://github.com/rsking84/pyeql

2. Create a branch for your work. We loosely follow the branching guidelines
   outlined at http://nvie.com/posts/a-successful-git-branching-model.

   If you are adding **documentation** or **bug fixes**, start with the **master** branch and
   prefix your branch with "fix-" or "doc-" as appropriate::

    git checkout -b fix-myfix master

    git checkout -b doc-mydoc master

   If you are adding a **new feature**, start with the **develop** branch and prefix your
   branch with "feature-"::

    git checkout -b feature-myfeature develop

3. Hack away until you're satisfied.

4. Push your work to the pyEQL repository::

    git push origin feature-myfeature

5. Create a pull request with your changes. See the tutorial at https://yangsu.github.io/pull-request-tutorial for instructions.

