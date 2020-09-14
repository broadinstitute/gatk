# ml4h
`ml4h` is a project aimed at using machine learning to model multi-modal cardiovascular
time series and imaging data. `ml4h` began as a set of tools to make it easy to work
with the UK Biobank on Google Cloud Platform and has since expanded to include other data sources
and functionality.   


Getting Started
* [Setting up your local environment](#setting-up-your-local-environment)
* [Setting up a remote VM](#setting-up-a-remote-vm)
* Modeling/Data Sources/Tests [(`ml4h/DATA_MODELING_TESTS.md`)](ml4h/DATA_MODELING_TESTS.md)
* [Contributing Code](#contributing-code)

Advanced Topics:
* Tensorizing Data (going from raw data to arrays suitable for modeling, in `ml4h/tensorize/README.md, TENSORIZE.md` )

## Setting up your local environment

Clone the repo
```
git clone git@github.com:broadinstitute/ml.git
```
Make sure you have installed the [Google Cloud SDK (gcloud)](https://cloud.google.com/sdk/docs/downloads-interactive). With [Homebrew](https://brew.sh/), you can use
```
brew cask install google-cloud-sdk
```

If you don't have your gcloud already configured -- set the project to broad-ml4cvd (This step is optional, scripts will also set projects directly when you run them.)

```gcloud config set project broad-ml4cvd```


### Conda (Python package manager)
* Download onto your laptop the Miniconda `bash` or `.pkg` installer for `Python 3.7` and `Mac OS X`
from [here](https://conda.io/en/latest/miniconda.html), and run it. If you installed Python via a package manager
such as `Homebrew`, you may want to uninstall that first, to avoid potential conflicts.
* On your laptop, at the root directory of your `ml4h` GitHub clone, load the `ml4h` environment via
    ```
    conda env create -f env/ml4h_osx64.yml
    ```
    If you get an error, try updating your `Conda` via
    ```
    sudo conda update -n base -c defaults conda
    ```
    If you have get an error while installing gmpy, try installing gmp:
    ```
    brew install gmp
    ```
    The version used at the time of this writing was `4.6.1`.

    If you plan to run jupyter locally, you should also (after you have `conda activate ml4h`, run `pip install ~/ml` (or wherever you have stored the repo)
* Activate the environment:
    ```
    source activate ml4h
    ```
You may now run code on your `Terminal`, like so
```
python recipes.py --mode ...
```
**Note** that *recipe*s require having the right input files in place and running them without proper inputs will not
yield meaningful results.   

### PyCharm (Python IDE if interested)
* Install PyCharm either directly from [here](https://www.jetbrains.com/pycharm/download/#section=mac), or download
the [Toolbox App](https://www.jetbrains.com/toolbox/app/) and have the app install PyCharm. The latter makes
PyCharm upgrades easier. It also allows you to manage your JetBrains IDEs from a single place if you have multiple
(e.g. IntelliJ for Java/Scala).
* Launch PyCharm.
* (Optional) Import the custom [settings](https://drive.google.com/open?id=1YvNVgVEH-rzsCJtrJ0mCi1nyAxG8Xync) as
described [here](https://www.jetbrains.com/help/pycharm/exporting-and-importing-settings.html).
* Open the project on PyCharm from the `File` menu by pointing to where you have your GitHub repo.
* Next, configure your Python interpreter to use the Conda environment you set up previously:
    * Open `Preferences` from `PyCharm -> Preferences...`.
    * On the upcoming `Preferences` window's left-hand side, expand `Project: ml4h` if it isn't already.
    * Highlight `Project Interpreter`.
    * On the right-hand side of the window, where it says `Project Interpreter`, find and select your `python`
    binary installed by `Conda`. It should be a path like `~/conda/miniconda3/envs/ml4h/bin/python` where `conda`
    is the directory you may have selected when installing `Conda`.
    * For a test run:
        * Open `recipes.py` (shortcut `Shift+Cmd+N` if you imported the custom settings).
        * Right-click on `if __name__=='__main__'` and select `Run recipes`.
        * You can specify input arguments by expanding the `Parameters` text box on the window
         that can be opened using the menu `Run -> Edit Configurations...`.    

## Setting up a remote VM
To create a VM without a GPU run:
```
./scripts/vm_launch/launch_instance.sh ${USER}-cpu
```
With GPU (not recommended unless you need something beefy and expensive)
```
./scripts/vm_launch/launch_dl_instance.sh ${USER}-gpu
```
This will take a few moments to run, after which you will have a VM in the cloud.  Remember to shut it off from the command line or [console](https://console.cloud.google.com/compute/instances?project=broad-ml4cvd) when you are not using it!  

Now ssh onto your instance (replace with proper machine name, note that you can also use regular old ssh if you have the external IP provided by the script or if you login from the GCP console)
```
gcloud --project broad-ml4cvd compute ssh ${USER}-gpu --zone us-central1-a
```

Next, clone this repo on your instance (you should copy your github key over to the VM, and/or if you have Two-Factor authentication setup you need to generate an SSH key on your VM and add it to your github settings as described [here](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/#platform-linux)):
```
git clone git@github.com:broadinstitute/ml.git
```

Because we don't know everyone's username, you need to run one more script to make sure that you are added as a docker user and that you have permission to pull down our docker instances from GCP's gcr.io. Run this while you're logged into your VM:
```
./ml/scripts/vm_launch/run_once.sh
```

Note that you may see warnings like below, but these are expected:
```
WARNING: Unable to execute `docker version`: exit status 1
This is expected if `docker` is not installed, or if `dockerd` cannot be reached...
Configuring docker-credential-gcr as a registry-specific credential helper. This is only supported by Docker client versions 1.13+
/home/username/.docker/config.json configured to use this credential helper for GCR registries
```

You need to log out after that (`exit`) then ssh back in so everything takes effect.


### Finish setting up docker, test out a jupyter notebook
Now let's run a Jupyter notebook.  On your VM run:

```
${HOME}/ml/scripts/jupyter.sh -p 8889
```
Add a ```-c``` if you want a CPU version.

This will start a notebook server on your VM. If you a Docker error like
```
docker: Error response from daemon: driver failed programming external connectivity on endpoint agitated_joliot (1fa914cb1fe9530f6599092c655b7036c2f9c5b362aa0438711cb2c405f3f354): Bind for 0.0.0.0:8888 failed: port is already allocated.
```
overwrite the default port (8888) like so
```
${HOME}/ml/scripts/dl-jupyter.sh 8889
```
The command also outputs two command lines in red.
Copy the line that looks like this:
```
ssh -i ~/.ssh/google_compute_engine -nNT -L 8888:localhost:8888 <YOUR VM's IP ADDRESS>
```
Open a terminal on your local machine and paste that command.  

If you get a public key error run: `gcloud compute config-ssh`

Now open a browser on your laptop and go to the URL `http://localhost:8888`

## Contributing code

Want to contribute code to this project? Please see [CONTRIBUTING](./CONTRIBUTING.md) for developer setup and other details.
