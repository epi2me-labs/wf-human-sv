## Introduction

The pipeline performs the following steps:
* Maps reads using [minimap2](https://github.com/lh3/minimap2)
* Calls variants using [sniffles2](https://github.com/fritzsedlazeck/Sniffles)
* Filters variants by minimum/maximum length, read support, or type (e.g. insertion, deletion, etc.)
* Optionally evaluates yours calls against a truthset using [truvari](https://github.com/spiralgenetics/truvari)
