# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.0]
### Changed
- Improved `--help` command
- Improved robustness of `fastq_ingress` utility
- `--bai` no longer needs to be provided with `--bam` if in the standard location (`<bam_path>.bai`)
- Replaced `lra/cutesv` workflow with `minimap2/sniffles2`
### Added
- `--tr_bedfile` option to provide a tandem repeat BED for your reference to improve the accuracy of SV calling
### Removed
- Test data is now hosted in AWS S3
### Fixed
- Set `out_dir` option format to ensure output is written to correct directory on Windows
- Prevent pipeline from stalling on Launcher when not providing a `--target_bedfile`

## [v0.0.6]
### added
- cumulative distribution plots

## [v0.0.5]
### Changed
- Report supports human references where chromosome names start with 'chr'

## [v0.0.4]
### Added
- Support for inputting single FASTQ files
- `--bam` and `--bai` options for starting from BAM files instead of FASTQ

### Fixed
- Now correctly overlap the target bedfile or generated bedfile with the truthset
- Remove benchmarkUseTruthsetBed option

## [v0.0.3]
### Added
- Multi-sample analysis is now supported

### Changed
- Standardised the CLI and reporting with other epi2me-labs workflows
- Benchmarking now uses an auto-generated bedfile from ref chromosomes by default
- Added benchmarkUseTruthsetBed to permit using truth-set bed for benchmarking instead

## [v0.0.2]
### Fixed
- Mosdepth limited to the target region if provided
- Correct conda profile environment file path
- Workflow internal reference naming
### Changed
- Pinned lra to version 1.1.2
- Make sample name optional, include in output names
### Added
- Additional summary stats and plots to report

## [v0.0.1]
First release.
