# Security Policy

## Supported versions

goSPL uses CalVer (`YYYY.M.D`) and only the latest release on PyPI and conda is supported with security fixes. Please upgrade before reporting an issue.

## Reporting a vulnerability

Please do not open a public issue for security problems.

Report privately using GitHub's [private vulnerability reporting](https://github.com/Geodels/gospl/security/advisories/new).

Include the goSPL version, a description of the issue, and steps or an input file to reproduce it.

You can expect an acknowledgement within 7 days. If the report is confirmed, we will aim to release a fix within 90 days and credit the reporter unless they prefer otherwise.

## Scope

goSPL is a research code run locally or on HPC systems. Relevant reports include unsafe handling of input files (YAML configuration, meshes, forcing data) and vulnerabilities in the build, release or container images. Numerical or scientific bugs should go through the regular [issue tracker](https://github.com/Geodels/gospl/issues).
