# omics-toolkit
Collection of scripts and utilities to analyze various types of 'omics data



## Installation

#### Github

To install directly from github, first clone this repository.

```bash
git clone git@github.com:BakelLab/omics-toolkit.git
```

Next, add the subdirectories within the `omics-toolkit` repository to your `PATH` environment variable. You can easily collect the ':' delimited string of all subdirectories by executing the following command from within the repository folder:

```bash
realpath */ | xargs | perl -pe 's/ /:/g'
```



