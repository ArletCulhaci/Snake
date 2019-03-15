# RAD-Seq analysis pipeline using a reference genome
What does this software do in simple words?
What are the expected outputs?

### Getting Started/Requirements/Prerequisites/Dependencies
This Snakemake pipeline make use of the [conda package manager](https://docs.conda.io/en/latest/) to install softwares and dependencies.
1. First, make sure you have conda installed on your system. Use [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) and follow the [installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).  
2. Using `conda`, create a virtual environment called `snakemake` to install Snakemake (version 5.4.3 or higher) by executing the following code in a Shell window: `conda create --name snakemake -c bioconda snakemake=5.4.3`. This will install `snakemake version 5.4.3` in a new environment called __snakemake__.
3. You can now run the pipeline since snakemake will use conda to install softwares and packages for each rule.

### Dry and real runs
For a dry run, type: `snakemake --use-conda -np`  
For a real run type: `snakemake --use-conda`  


# For improvement: a template for your README file
Taken from [Zalando](https://github.com/zalando/zalando-howto-open-source/edit/master/READMEtemplate.md)

Clear documentation is critical to the success of your project. This checklist is meant to help you cover all your bases. Not every section/subsection will be relevant to your project; pick and choose what is. Inspired by READMEs of very successful projects like etcd.

Please copy-paste this into a new document and save as you build your READMEs. For alternative formats, you might create a [Structured README](https://github.com/shaloo/structuredreadme), which offers a thorough breakdown of optional README ingredients for you to consider. You might also take a look at [this similar checklist](https://github.com/cfpb/open-source-project-template); or check out [art-of-readme](https://github.com/noffle/art-of-readme).

### Project Name/Intro

- Describe very briefly but clearly what the project does.
- State if it is out-of-the-box user-friendly, so it’s clear to the user.
- List its most useful/innovative/noteworthy features.
- State its goals/what problem(s) it solves.
- Note and briefly describe any key concepts (technical, philosophical, or both) important to the user’s understanding.
- Link to any supplementary blog posts or project main pages.
- Note its development status.
- Include badges.
- If possible, include screenshots and demo videos.

### Core Technical Concepts/Inspiration

- Why does it exist?
- Frame your project for the potential user.
- Compare/contrast your project with other, similar projects so the user knows how it is different from those projects.
- Highlight the technical concepts that your project demonstrates or supports. Keep it very brief.
- Keep it useful.

### Getting Started/Requirements/Prerequisites/Dependencies
Include any essential instructions for:
- Getting it
- Installing It
- Configuring It
- Running it

### More Specific Topics (+ sample sub-categories)
- Versioning: Services, APIs, Systems
- Common Error Messages/related details
- Tests
- Is it a Swift project? Please take a look at Mattt Thompson & Nate Cook's [Swift documentation](http://nshipster.com/swift-documentation/) guide

### Contributing
- Contributor Guidelines
- Code Style/Requirements
- Format for commit messages
- Thank you (name contributors)

### TODO
- Next steps
- Features planned
- Known bugs (shortlist)

### Contact
- Email address
- Google Group/mailing list (if applicable)
- IRC or Slack (if applicable)

### License
