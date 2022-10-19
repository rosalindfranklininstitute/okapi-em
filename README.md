# okapi-em

https://github.com/rosalindfranklininstitute/okapi-em

<!--
[![License](https://img.shields.io/pypi/l/okapi-em.svg?color=green)](https://github.com/rosalindfranklininstitute/okapi-em/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/okapi-em.svg?color=green)](https://pypi.org/project/okapi-em)
[![Python Version](https://img.shields.io/pypi/pyversions/okapi-em.svg?color=green)](https://python.org)
[![tests](https://github.com/perdigao1/okapi-em/workflows/tests/badge.svg)](https://github.com/rosalindfranklininstitute/okapi-em/actions)
[![codecov](https://codecov.io/gh/perdigao1/okapi-em/branch/main/graph/badge.svg)](https://codecov.io/gh/rosalindfranklininstitute/okapi-em)
[![napari hub](https://img.shields.io/endpoint?url=https://api.napari-hub.org/shields/okapi-em)](https://napari-hub.org/plugins/okapi-em)
-->

A napari plugin for processing serial-FIB-SEM data.

Powered by [chafer] and [quoll].

This [napari] plugin contains the following tools:

    - Two charge artifact suppression filter
        - directional fourier bandapass filter
        - line-by-line filter function optimiser and subtraction (requires charge artifact labels)
    - slice alignment using constrained SIFT
    - FRC estimation

----------------------------------

This [napari] plugin was generated with [Cookiecutter] using [@napari]'s [cookiecutter-napari-plugin] template.

<!--
Don't miss the full getting started guide to set up your new package:
https://github.com/napari/cookiecutter-napari-plugin#getting-started

and review the napari docs for plugin developers:
https://napari.org/plugins/stable/index.html
-->

## Installation

You can install `okapi-em` via [pip]:

    `pip install okapi-em`

For development mode it can be installed by navigating to the cloned `okapi-em` folder and run:

    `pip install -e .`

This should install in any machine, however ...

Currently the FRC calculation provided by the [quoll] package which is optional because
of its stringent environmemt requirements from miplib package. These currently are:
    - python 3.7
    - linux OS

This issue will be addressed in future version.


To install okapi-em with quoll included:
    
    `pip install okapi-em[all]`

Note that to run napari in python 3.7 you will need to use the command:

    `python -m napari`



## Contributing

Contributions are very welcome. Tests can be run with [tox], please ensure
the coverage at least stays the same before you submit a pull request.

## License

Distributed under the terms of the [Apache Software License 2.0] license,
"okapi-em" is free and open source software

## Issues

If you encounter any problems, please [file an issue] along with a detailed description.

[quoll]: https://github.com/rosalindfranklininstitute/quoll
[chafer]: https://github.com/rosalindfranklininstitute/chafer
[napari]: https://github.com/napari/napari
[Cookiecutter]: https://github.com/audreyr/cookiecutter
[@napari]: https://github.com/napari
[MIT]: http://opensource.org/licenses/MIT
[BSD-3]: http://opensource.org/licenses/BSD-3-Clause
[GNU GPL v3.0]: http://www.gnu.org/licenses/gpl-3.0.txt
[GNU LGPL v3.0]: http://www.gnu.org/licenses/lgpl-3.0.txt
[Apache Software License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[Mozilla Public License 2.0]: https://www.mozilla.org/media/MPL/2.0/index.txt
[cookiecutter-napari-plugin]: https://github.com/napari/cookiecutter-napari-plugin


[tox]: https://tox.readthedocs.io/en/latest/
[pip]: https://pypi.org/project/pip/
[PyPI]: https://pypi.org/
