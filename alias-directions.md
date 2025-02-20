# Tired of typing long directory paths?

We are always navigating to the same directory,let's create an alias so you can jump there with a short command.

## Create a temporary alias (Works Only In The Current Session)
``` bash
alias myocean="cd /ocean/projects/agr250001p/your-psc-username"
```
- Now you can go to your personal folder in ocean by typing `myocean`. This alias will work only until you log out.

## Create a permanent alias (Works Every Time You Log In!)
These two command lines will create your alias forever on this HPC!
``` bash
echo 'alias myocean="cd /ocean/projects/agr250001p/your-psc-username"' >> ~/.bashrc
source ~/.bashrc
```
- NNow, every time you log onto the HPC, you can access your folder by typing `myocean`
