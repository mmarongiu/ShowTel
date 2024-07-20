# ShowTel
Tool for the real-time monitoring of the operations of a radio telescope during an observing session, in particular: (1) the status of the antenna control system, (2) general information about the observing session, (3) the angular distance between the Sun (and the Moon) and the radio telescope pointing, (4) the sky position of the Sun, Moon and the astrophysical object under observation, and (5) the Radio Frequency Interference (RFI) of the sky region pointed by the radio telescope, at the running observing frequency.



## Installation



### Preparation and dependencies



#### Anaconda and virtual environment (recommended)
We strongly suggest to install the
[Anaconda](https://www.anaconda.com) Python distribution.
Once the installation has finished, you should have a working `conda`
command in your shell. First of all, create a new environment:

    $ conda create -n py3 python=3

load the new environment:

    $ source activate py3

and install the dependencies (including a few optional but recommended):

    (py3) $ conda install numpy, astropy, datetime, pytz, time, sched, threading, socket, os, tkinter, tkmacosx, PIL, gtts

and
    (py3) $ sudo apt-get install ffmpeg libavcodec-extra
    (py3) $ pip install pydub



##### Download of the de441.bsp ephemeris
You must to download the de441.bsp ephemeris following this procedure in ipython environment:
    In [1]: from astropy.coordinates import solar_system_ephemeris

and
    In [2]: solar_system_ephemeris.set("ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de441.bsp")

Once you downloaded this file (3.3 GB), you must put this file in the directory "utilities".



### Cloning and execution

Clone the repository:

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://github.com/mmarongiu/ShowTel.git

or if you have deployed your SSH key to Github:

    (py3) $ git clone git@github.com:mmarongiu/ShowTel.git

Then:

    (py3) $ cd Showtel
    (py3) $ python showtel_v03.py


### Updating

To update the code, simply run `git pull` and reinstall:

    (py3) $ git pull


### Contribution guidelines

See the file CONTRIBUTING.md for more details.

This code is written in Python 3.11+. Tests run at each commit during Pull Requests, so it is easy to single out points in the code that break this compatibility.


### If you use this code

First of all... **This code is under development!**... so, it might well be that something does not work as expected. For any inquiries, bug reports, or suggestions, please use the [Issues](https://github.com/mmarongiu/ShowTel/issues) page.

If you used this software package to reduce data for a publication, please write in the acknowledgements something along these lines:

    This work makes use of the ShowTel tool (https://github.com/mmarongiu/ShowTel)
