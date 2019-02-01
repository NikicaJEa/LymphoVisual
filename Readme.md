
## LymphoVisual - Clonal Rearrangement Visualisation Tool##
----------
This bioinformatic tool helps us to  visualise the generated LymphoTrack data (invivoscribe, San Diego, USA). The clone rearrangement frequency graph helps in detecting the proliferative activity of lymphocytes. 
Furthermore, the comparison graph compares the clonal rearrangements of two patients. Those visualisation tools could be helpful in a lymphoma diagnosis and used in monitoring the lymphoma patients.

## Getting Started ##
----------
**Prerequisites**
Programs: GNUPLOT (5.2 patchlevel 4)
Libraries for the environment : see requirements.txt

**Installing GNUPLOT**

apt-get install gnuplot

** Setting up the environment **

Note: make sure that none of the files have a space in their name!

1. install virtualenv
   sudo apt-get install virtualenv
   sudo apt-get update
2. install pip3
   sudo apt-get install python3-pip

3. create your environment folder   
   virtualenv -p python3 --no-site-packages (target folder)
   
4. activate the environment
   source (target folder)bin/activate
   
5. install the libraries
   pip3 install -r requirements.txt
   
6. pip3 instal matplotlib==2.2.2

7. pip3 install pandas

8. pip3 install plotly==2.7.0

 

**Running the test**

Unpack the whole file, install the necessary programs and setup the environment. Finally, run the runme.sh script which has to be in the same folder as the Test_Data

## Build With ##
----------
 
 - Python 3

## Contributing to IGH FR1, IGH FR2, IGK and TRG Clonal Rearrangement Visualisation Tool ##
----------
LymphoVisual is open-source software under the term of an [GNU GPL v3 licence.](https://www.gnu.org/licenses/gpl-3.0.html) The source code is hosted on GitHub.

LymphoVisual
 was created by Nikola Vinko with a help of  Karl Kashofer and Stefan Sauer.