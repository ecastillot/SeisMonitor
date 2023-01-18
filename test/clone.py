import os
from git import Repo

def git_clone_aipicker(name,repo_dir):
    """
    Params:
    -------
    name: str
        EQTransformer or PhaseNet
    repo_dir: str
        Directory path to place the repository
    """
    if name == "PhaseNet":
        git_url = "https://github.com/ecastillot/PhaseNet.git"
    elif name == "EQTransformer":
        git_url = "https://github.com/ecastillot/EQTransformer.git"
    else:
        return False

    repo_dir = os.path.join(repo_dir,name)
    if os.path.isdir(repo_dir):
        print("There is alaready")
        return True
    else:
        Repo.clone_from(git_url, repo_dir)
        return True

git_clone_aipicker(name="EQTransformer",repo_dir="/home/emmanuel")