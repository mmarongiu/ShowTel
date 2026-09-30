# ShowTel
Tool for the real-time monitoring of the operations of a radio telescope during an observing session, in particular: (1) the status of the antenna control system, (2) general information about the observing session, (3) the angular distance between the Sun (and the Moon) and the radio telescope pointing, (4) the sky position of the Sun, Moon and the astrophysical object under observation, and (5) the Radio Frequency Interference (RFI) of the sky region pointed by the radio telescope, at the running observing frequency.



## Installation



### Preparation and dependencies



#### Pyenv and virtual environment (recommended)
We strongly suggest to install the
[Pyenv](https://github.com/pyenv/pyenv) Python distribution.
Once the installation is complete, you should create a new Python environment:

    $ pyenv install 3.13.4



##### Download of the de441.bsp ephemeris
You must to download the de442s.bsp ephemeris from this link:

    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/

Once you downloaded this file (31 MB), you must put this file in the directory "utilities".



### Cloning and execution

Clone the repository:

    $ cd /my/software/directory/
    $ git clone https://github.com/mmarongiu/ShowTel.git

or if you have deployed your SSH key to Github:

    $ git clone git@github.com:mmarongiu/ShowTel.git

Then go to the main directory and type:

    $ pip install .


### Updating

To update the code, simply run `git pull` and reinstall:

    $ git pull


### Contribution guidelines

See the file CONTRIBUTING.md for more details.

This code is written in Python 3.13+. Tests run at each commit during Pull Requests, so it is easy to single out points in the code that break this compatibility.


### If you use this code

First of all... **This code is under development!**... so, it might well be that something does not work as expected. For any inquiries, bug reports, or suggestions, please use the [Issues](https://github.com/mmarongiu/ShowTel/issues) page.

If you used this software package to reduce data for a publication, please write in the acknowledgements something along these lines:

    This work makes use of the ShowTel tool (https://github.com/mmarongiu/ShowTel)