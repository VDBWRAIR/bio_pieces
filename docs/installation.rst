Installation
============

It is recommended to install into a virtualenv. If you know what you are doing
and don't want to install into virtualenv, then you can skip right to step 3

1. Setup Virtualenv

    It is assumed you have virtualenv already installed. If not see
    https://virtualenv.pypa.io/en/latest/installation.html

    .. code-block:: bash

        virtualenv env

2. Activate virtualenv

    .. code-block:: bash

        . env/bin/activate

3. Install dependencies

    .. code-block:: bash

        pip install -r requirements.txt

   For python 2.6 you will need to also install some additional packages

        .. code-block:: bash

            pip install -r requirements-py26.txt

4. Install bio_bits

    .. code-block:: bash

        python setup.py install
