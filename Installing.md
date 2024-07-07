# Installing Tutorials

This Tutorial have been put together by Jean-Christophe LONE.

In this tutorial we will provide help to install Python, R and Rstudio. We will cover quickly the creation of Python env. 

## Downloading the material from GitHub 

The code for 'HitchhikerGuide' is available online at GitHub.
<details>
  <summary> GitHub for Beginners </summary>
  Git is a version control system that helps in storing the source code of a project and tracking all changes made to that code. It enables developers to collaborate more effectively by managing changes from multiple contributors. GitHub enhances this experience by offering a hosting service and web interface for Git repositories, along with tools for collaboration. Think of GitHub as a social networking site for software developers, fostering community and collaboration.
</details>

To download the student code, visit the following link: `https://github.com/iMM-Workshops/2024_Hitchhiker_Guide_scRNA-seq`

1. Click the "Clone or download" button to get the lab tools in a zip format: `2024_Hitchhiker_Guide_scRNA-seq-main.zip`
2. Unzip the file in your desired working directory (eg Desktop). This folder will serve as your workspace for the project.

**_NOTE: Click on the arrow to open collapsible section_** 

<details>
  <summary> Workspace Setup for Beginners </summary>
  Unzip the file to your chosen working directory (e.g., your desktop). This workspace is like a physical desk where you keep everything you need for a task. Having all necessary files and tools in one place allows for easy access and use. If something is missing, it's like having to retrieve a pen from a drawer—additional effort is needed to continue your work.
</details>




## Installing Python from Scratch with miniconda

If you are new to Python, we recommend using miniconda. You can find installers for Windows, Linux, and OSX platforms here : `https://docs.anaconda.com/miniconda/`

<details>
  <summary> Miniconda for Noobs </summary>
Installing Python with Miniconda is like setting up a versatile workbench in your digital workshop. Miniconda serves as the sturdy foundation, providing the essential tools: conda and Python, much like a well-built bench with a few basic, yet crucial, instruments. From this solid starting point, you can easily add more specialized tools as your projects require, just like hanging additional gadgets and equipment on a pegboard above your workbench. It's a lightweight, flexible, and powerful setup, transforming your workspace into a hub of endless possibilities, ready to support your creative and technical endeavors.
</details>


### Windows

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe
start /wait "" miniconda.exe /S
del miniconda.exe
```

### MacOsX using curl

```bash
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

```bash
# Initiate the bash shells
~/miniconda3/bin/conda init bash
```

### MacOsX using homebrew

You can also do it with `homebrew`you can do it as well
```bash
# Install homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install --cask miniconda
brew install virtualenv
```

<details>
  <summary> Homebrew for Beginners </summary>
  Imagine Homebrew as a magical toolbox that helps you install software programs with just a few simple commands. It's like having a friendly wizard who knows all the best tools and can quickly fetch them for you. Homebrew works on macOS and Linux, making it easy to set up new programs without any hassle. The name "Homebrew" is like brewing a custom potion just for your computer, tailored to your needs and preferences.
</details>


## Installing with Pip (create environment)

If you are familiar with Python, you might already know about the pip package installer. Using pip and a virtual environment can make it easier for you to manage your Python packages without causing conflicts with existing Python installations.

Steps to Install and Create a Virtual Environment:

1. Install virtualenv: This tool helps you create a separate space for your Python projects, just like having a special playroom where you can keep all your toys organized without mixing them up with the rest of the house.

```bash
# Install the virtualenv
sudo pip install virtualenv
```

2.  Create a Virtual Environment: Once you have virtualenv installed, you can create a new environment for your project. This is like setting up your own secret clubhouse where you have everything you need for your project.

```bash
# Create a new virtual environment
virtualenv myenv
```

3. Activate the Virtual Environment: Before you start installing packages, you need to enter your clubhouse. This step ensures that all the packages you install will stay inside this special environment.

```bash
# Activate the virtual environment
source myenv/bin/activate  # On Windows use `myenv\Scripts\activate`
```

4. Install Packages with Pip: Now that you are in your virtual environment, you can use pip to install the packages you need for your project. It's like stocking your clubhouse with all the cool gadgets and tools you need to get started.

```bash
# Install packages using pip
pip install package_name
```

By following these steps, you can easily manage your Python projects without worrying about messing up other Python installations on your computer.

<details>
  <summary> Virtual environment for Noobs </summary>
Think of a virtual environment like a special playroom just for your toys and games. When you create a virtual environment in Python, it's like setting up this playroom on your computer. This special space keeps all your project tools and packages organized and separate from everything else on your computer, so nothing gets mixed up. Just like you keep your toys in your playroom, a virtual environment keeps your Python projects neat and tidy, making it easier and safer to work on them without messing up other projects.
</details>



## Installing R and Rstudio

### R from CRAN Download from CRAN 

You will need to install `R`, the required `R` libraries, and download the codebase. It is crucial to complete this before the Workshop starts to avoid unnecessary delays. Please follow the installation instructions provided below.

If you are new to `R`, the best option is to install it from the CRAN Website: 
```
https://cran.r-project.org/
```

There, you can find installers for your operating system (Windows, Linux, and macOS).


<details>
  <summary> R for Noobs </summary>
R is a language and environment for statistical computing and graphics. R is similar to a magic box of tools for solving puzzles. These puzzles could be anything from math problems to understanding data from a science experiment. R is a special computer language that helps you use these tools to work on different projects. Just like you might have different kinds of toys for different games, R has different tools and packages for different kinds of data puzzles.
</details>


### Rstudio as IDE from Posit

An Integrated Development Environment (IDE) includes a text editor and various tools to debug and interpret complex code. Think of `R` as the engine of a car, and the IDE as the components that make using the engine easier (like the wheels, seat, or radio).

We recommend using RStudio for this workshop:
```
https://posit.co/download/rstudio-desktop/
```

However, feel free to use the software you are most comfortable with. `Visual Studio Code` (`VSCode`) or `Jupyter Notebook` are also good options.

<details>
  <summary> Rstudio for Noobs </summary>
Now, imagine RStudio as your super organized toy room where you can easily find and use all the tools from your R magic box. RStudio helps you keep everything in one place, making it easier to write code, see results, and keep track of your projects. It's like having a special desk where you can draw, build, and create, with all your favorite supplies right at your fingertips. With RStudio, using R becomes much more fun and easy, just like playing in your perfect toy room!
</details>





