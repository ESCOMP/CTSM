Any netCDF files in this directory (and elsewhere in the repository) should be
stored via Git LFS (Large File Support). To retrieve any netCDF files, or to add
or modify any netCDF files, you will need to install the Git LFS tool.

Installing Git LFS on your machine is a two-step process; step (1) needs to be
done once per machine, and step (2) needs to be done once per user:
1. Install the Git LFS tool: Follow the instructions on the [Git LFS
page](https://git-lfs.github.com/) for installing Git LFS on your platform.
  - On derecho the system default version of git already has Git LFS installed.
  - On cheyenne and casper, Git LFS is already available as long as you are using a git
    module rather than the default system-level git. So just make sure that you
    are always using git via a git module (`module load git`).
  - On a Mac using homebrew, this can be done with `brew install git-lfs`.
2. Set up your git configuration to use Git LFS by running: `git lfs install`
   (this will add a few lines to your `.gitconfig` file in your home directory).
   (If you're not sure whether you have already done this, it is safe to rerun
   that command to be sure.)

If any netCDF file appears to be a text file like this:

```
version https://git-lfs.github.com/spec/v1
oid sha256:df03c76138bbdcdd12497a4bb0a86591283d1d8ef90cba80d823660f549381c0
size 896
```

then that is a sign that you may not have Git LFS installed. You can first try
running "git lfs pull". If that doesn't work, then install Git LFS as documented
above.

When you have new large files to commit, use the usual syntax of adding and
committing:
git add <filename>
git commit

For more information on using Git LFS with CTSM, search for lfs here:
<https://github.com/ESCOMP/CTSM/wiki/Directions-for-editing-CLM-documentation-on-github-and-sphinx>.
Most of the notes about using Git LFS for the documentation's image files apply
here as well, although (in contrast to image files) netCDF files *are* pulled
down by default, so you usually should *not* need to explicitly run `git lfs
pull` to get new or updated netCDF files.
