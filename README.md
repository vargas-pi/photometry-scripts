# Alhadeff Lab Photometry Scripts

This repository contains some of the code we use in lab for online viewing of photometry signals and post-hoc analyses. The code for viewing the data is specific to TDT systems but all code for post-hoc analyses should be generalizable so long as you can load the data in the appropriate format. The goal of this code is to provide a user-friendly interface for performing fiber photometry analysis to remove the burden from researchers of writing custom scripts for general pre-processing steps, while allowing those with some knowledge of Python the flexibility to make the analyses their own.

## Installation 
This code base is not currently pip or conda installable (hopefully in the future it will be!). As a result, before using the pipeline it is important that you first have the repository downloaded locally. Feel free to either clone into this repository through git or simply download and extract a zip file. If you are familiar with using git, this may be the preferable approach a it will allow easy access to future updates. 

Once  you have the code downloaded, we recommend setting up an environment for running the code. This makes it so that you have an isolated setting for running your analyses with all the correct versions of the necessary packages for running things. The requirements are not particularly strict in this case so this step is not entirely necessary so long as you have already have the appropriate versions of the dependencies listed in the file **environment.yml**. This is however good programming practice so I'll briefly detail the steps for doing this below. If you use CPython and are more familiar with virtualenv feel free to use that instead. We will go over the steps for creating an environment from a .yml file using the [Anaconda](https://www.anaconda.com/products/distribution) distribution of Python below. 

1. Open a conda enabled terminal or command prompt. You can launch this from the Anaconda Navigator on Windows. On Mac, your terminal should already be Anaconda enabled.

2. Navigate to the photometry scripts folder. This can be done with the command `cd path\to\photometry-scripts`

3. Run the command `conda env create --file environment.yml`

4. Whenever you're ready to run an analysis just be sure to activate this environment first by running `conda activate photometry`

Optionally, if you will be using Jupyter Notebooks to run any post-hoc analyses, it may be useful to create a designated kernel for any notebooks you create for these analyses. Kernels are basically the programs behind the scenes that actually runs the code you type in the Noteboks. When you install Jupyter Notebooks you'll have a default kernel that runs code in your base conda environment. Creating a kernel specifically for this new environmnet just makes it easier to ensure your notebook is running inside the photometry environment you've created above. To do this, you can run the following commands 

```
conda activate photometry
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=photometry
```



## Workflow
The overall workflow for how we use this code in lab and the corresponding files in the repository are as follows:

1. `import_plot.py` - visualize data during data collection
2. `import_plot.py` - export the data to an npy file which contains a list of the recordings stored in our custom data structure called a mouse_data object
3. `analysis_pipeline_cli.py` | `analysis_pipeline.py`  - run post-hoc analyses on the exported files (normalization and visualization) either using the provided command-line interface or importing the necessary functions into a Jupyter Notebook

With the exception of the Jupyter Notebooks method of runnin post-hoc analyses, everything that follows is run from the command-line. Before performing any of these tasks, be sure to activate the photometry environment created above, and navigate to the folder with the photometry code. The workflow we'll go through below assumes the user is collecting data with a TDT system. In the next section we'll discuss how you may prep any data to be fed through the pos-hoc analaysis pipeline.

### Online Data Visualization
In order to use the **import_plot.py** script for anything, you need to first create a file in your copy of the repository called **sys_parameters.py**. An example of what should go inside of this file is provided in **example_sys_parameters.py**. The idea is we need to set some important system specific parameters so we know where the data is being streamed to in realtime and what field the data is stored under in the output files. Once this file is prepared, and after you've started a recording, you can stream and plot the data by running the following command.

```
python import_plot.py -mouse1 ntn123 -mouse2 ntn234
```

For newcomers to Python, `--arg value` or `-arg value` are just ways of specifying and argument for a given script. To see all the possible arguments for import_plot.py simply run `python import_plot.py --help`. The important ones here are `-mouse1` and `-mouse2`. At least one of these needs to be specified. 1 and 2 refer to the line that the mouse is being run on, corresponding to however you've labeled the channels in the `sys_parameters.py` file.

### Exporting the data
Exporting data uses the same **import_plot.py** script but with different arguments. The key difference here is we also need to tell the script the location of the raw data files using the `-path` argument. Notably, this is not exactly the file path you define in Synapse before starting the recording. What we need here is a few folders deeper. It should be the last subfolder of this folder. You will also need to a time for the stimulus or each recording your are exporting using the arguments `t_stim1` and `t_stim2`.

### Post-Hoc Analyses through the Command-Line Interface
We provide a command-line interface for running post-hoc analyses. To launch the CLI, simply run the command `python analysis_pipeline_cli.py`. From here, the CLI will prompt you to provide all the information it needs to run the analysis. 

### Post-Hoc Analyses through Jupyter Notebooks
See the notebooks folder for some sample notebooks that walk through the process of doing the same analyses inside of a Jupyter Notebook. For any of these analyses you will need to import all of the functions from the file **analysis_pipeline.py**. These are the same functions that are imported and run by the CLI. For those that plan to submit pull requests to this repository, We ask that you keep any custom notebooks,if you want to keep your notebooks inside this repository locally, we ask that you store your notebooks under **notebooks/custom-scripts** where it will not be tracked by git.

## Prepping External Data For Post-Hoc Analysis
As aforementioned, the pre-processing pipeline is specific to TDT systems. Future iterations of the code will hopefully be able to handle data from other sources implicitly. For now we can provide a blueprint for what would be needed to use our post-hoc analysis. The post-hoc analysis expects as input an npy file that contains a list of objects which are instances of the `mouse_data` object defined in **utilities.py**. As a result, in order to feed any external data into the pipeline, one will need to create such instances. The main things you will need to specify when doing this are as follows

* an identifier for the mouse
* a numpy array of the 465nm channel fluorescence values
* a numpy array of the 405nm channel fluorescence values
* the sampling rate
* a stimulus time for the recording
* *Optionally:* a condition name for this recording

Below is an example taking into account the order of the positional arguments to the class's init function
```
m = mouse_data("ntn123", F465, F405, fs, t_stim = 500, cond = "control" ):
```

## Contact
For any assistance feel free to reach out to [Nathaniel Nyema](mailto:nnyema@gmail.com) or [Alexandra Vargas](mailto:alexandragve@gmail.com)
