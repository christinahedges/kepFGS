kepFGS
======

Tools and tutorials for using the Kepler and K2 FGS data. The FGS.py file included in this repo contains tools to download, clean and stitch FGS quarters and campaigns. Currently this is only compatible with python 2.7.

Dependencies
------------

Currently FGS.py requires

>
* pandas 
* numpy
* urllib
* glob
* progressbar
* astropy
>

This list of dependencies may get smaller.

Example use
-----------

A full tutorial of how to use FGS data and the tools provided in FGS.py is given in the notebook **FGS_tutorial.ipynb**. A simple example is below.

	$ import FGS
	$ FGS.get_data() #Download all Kepler FGS data
	$ #generate light curve from all quarters for Kepler ID 2564891
	$ time,flux,column,row=FGS.gen_lc(datadir,ID=2564891)

