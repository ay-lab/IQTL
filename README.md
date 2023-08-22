# IQTL (Interaction QTL)

## Description
Quantitative trait loci (QTL) associated with chromatin interaction data.

	Although the current workflow is based on HiChIP data, in practice, the workflow can be adapted to process other types of chromatin interaction data, such as HiC, PCHi-C, Micro-C, ChIA-PET, etc.


## ********
## Installation
## ********

Following packages (and associated libraries) need to be installed to run IQTL pipeline and other processing scripts:

	1. R (version 3.6.1 or higher) along with the following libraries / packages:

		data.table, edgeR

	2. HiC-pro (https://github.com/nservant/HiC-Pro). A package to align the HiChIP datasets with respect to the reference genome. We recommend installing the latest version.

	3. FitHiChIP (https://github.com/ay-lab/FitHiChIP) and its associated dependencies, to call significant HiChIP loops. Check its documentation for the detailed installation. The significant HiChIP loops will be used as an input to the IQTL derivation pipeline.

	4. RASQUAL (https://github.com/natsuhiko/rasqual). Download the source code and install. The directory containing the downloaded package will be used as an input of the IQTL pipeline.

	5. GATK (https://gatk.broadinstitute.org/hc/en-us). Download and install. The path of the GATK executable needs to be provided as a configuration parameter.

	6. bedtools (https://bedtools.readthedocs.io/en/latest/)


## ********
## Step-by-step execution
## ********

## ===========
## A. Preprocessing
## ===========

## A.0. Creating the sample list and the metadata

Check the file *Data/DonorList_Annotated.txt* provided along with this repository. User should maintain the format of this file, specifying the list of samples and the associated metadata.







## ===========
## B. Running IQTL pipeline (SnakeMake)
## ===========















## A.4. Processing HiChIP alignment files (generated from HiCPro)

Check *Preprocessing/HiChIP_Alignment* for details on how to 


## ===========
## B. Derivation of IQTLs
## ===========

Check the section *IQTL_derive* for a step-by-step description on how to use the donor / sample specific FitHiChIP loops to derive the IQTLs.





## ********
## Support
## ********

For any queries, please e-mail: 

Sourya Bhattacharyya <sourya@lji.org>

Ferhat Ay <ferhatay@lji.org>







































## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.lji.org/sourya/iqtl.git
git branch -M master
git push -uf origin master
```

## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.lji.org/sourya/iqtl/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Automatically merge when pipeline succeeds](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thank you to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.



## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.


## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.


## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.



