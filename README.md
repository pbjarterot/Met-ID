<img src="https://github.com/pbjarterot/Met-ID/assets/46728406/115bcc2d-3c16-42ab-8f50-484b2dd5d253" width="200" height="200">

![example workflow](https://github.com/pbjarterot/Met-ID/actions/workflows/main.yml/badge.svg)    [![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    

# Met-ID
Metabolite identification in Mass Spectrometry Imaging.
Met-ID has been developed to automate metabolite identification in Mass Spectrometry Imaging, at the moment most of this is done manually by experts which in the world of high throughput studies is not feasable.

Met-ID has a particular focus on derivatizing matrices leading to other ions than the common [M+H]+ in positive mode and [M-H]- in negative mode. As [FMP-10](https://www.nature.com/articles/s41592-019-0551-3) was developed in house, it features heavily in the software, however, this is mostly to show the point at which to start as Met-ID is extendable to use any derivatizing matrix with the tools to do local version changes right from inside the software.

Met-ID also includes tools to study MS2 spectra, as a critical part of metabolite identification is to study MS2 spectra of unknown compounds and compare them to those generated by standard. The base version of Met-ID comes with a number of MS2 spectra collected from chemical standards with an FT-ICR using the FMP-10 chemical matrix. Met-ID allows the user to input MS2 spectra from unknown compounds to compare to standards and also ways to add user generated spectra from chemical standards to the database. It should be noted that all changes to the database only affect that specific install of Met-ID and not any other.


## Met-ID is under development
To help with debugging, a metid_log.log file on the desktop (it will be moved elsewhere in future versions to not clutter the computer), which can be sent together with the steps to recreate the error to patrik.bjarterot@uu.se for some technical support. For some reason, reinstalling newer versions of Met-ID does not update the database files it comes bundled with and as such the database files will have to removed manually as the Met-ID installer will add the correct files if there is no other files already there, this means removing the "com.farmbio.metid" folder found in "C:/Users/"username"/AppData/Roaming/" in windows, "/Users/"username"/Library/Application Support/" on macos, and a similar file on linux, "username" needs to be changed for the personal username used on your own computer, on macos, there may be an extra button to check to view the Library folder. Met-ID then needs to be reinstalled to have the intended database files.


## How to Install
Files to install Met-ID are located on the right hand side under "Releases" or can be found by clicking [here](https://github.com/pbjarterot/Met-ID/releases).

### Windows
To install Met-ID on Windows, download the .msi file found in releases and follow the instructions.

### MacOS
There are currently some issues with adding functional groups on MacOS. I am working on it as fast as I can.
To install Met-ID on MacOS, download the .dmg file found in releases and follow the instructions.

### Linux (Debian systems, e.g. Ubuntu)
To install Met-ID on Debian systems, download the .deb file found in releases and follow the instructions.

## Issues
To report any issue with the software, please use the "Issues" tab above.

## Discussions
To discuss anything regarding the software, please use the "Discussions" tab above.

## Preprint / Paper
Preprint:

https://www.biorxiv.org/content/10.1101/2025.01.24.634674v1


Paper:

https://pubs.acs.org/doi/full/10.1021/acs.analchem.5c00633

## Platforms 
[![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)](https://www.linux.org/)
[![macOS](https://img.shields.io/badge/mac%20os-000000?style=for-the-badge&logo=macos&logoColor=F0F0F0)](https://www.apple.com/se/macos)
[![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white)](https://www.microsoft.com/)

## Framework
[![Tauri](https://img.shields.io/badge/tauri-%2324C8DB.svg?style=for-the-badge&logo=tauri&logoColor=%23FFFFFF)](https://tauri.app/)

## CI/CD
[![GitHub Actions](https://img.shields.io/badge/github%20actions-%232671E5.svg?style=for-the-badge&logo=githubactions&logoColor=white)](https://github.com/pbjarterot/Met-ID/actions)

## Written in
[![Rust](https://img.shields.io/badge/rust-%23000000.svg?style=for-the-badge&logo=rust&logoColor=white)](https://www.rust-lang.org/)
[![TypeScript](https://img.shields.io/badge/typescript-%23007ACC.svg?style=for-the-badge&logo=typescript&logoColor=white)](https://www.typescriptlang.org/)
[![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)](https://www.python.org/)
[![HTML5](https://img.shields.io/badge/html5-%23E34F26.svg?style=for-the-badge&logo=html5&logoColor=white)](https://en.wikipedia.org/wiki/HTML5)
[![CSS3](https://img.shields.io/badge/css3-%231572B6.svg?style=for-the-badge&logo=css3&logoColor=white)](https://en.wikipedia.org/wiki/CSS)

## Contact
[![Linktree](https://img.shields.io/badge/linktree-1de9b6?style=for-the-badge&logo=linktree&logoColor=white)](https://linktr.ee/patrikbja)












